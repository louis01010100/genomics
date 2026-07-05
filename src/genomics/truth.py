"""Backbone-driven exact-match SNV truth calling.

`genomics snv-truth` consumes an `axiom backbone` output VCF as `--backbone-file`
and builds a per-sample truth VCF (+ TSV) by exact `(chrom, pos, ref, ALT-set)`
matching against a normalized single-sample VCF. Backbone records that share an `ID`
form a SNP family: if one member matches the sample it is emitted with the sample's
genotype and its siblings are dropped; if none matches, every member is emitted with a
depth- and gender-aware homozygous-reference or no-call genotype.

Input VCFs (bgzip-compressed and indexed) may be split by chromosome — a sample's
calls spread across several per-chromosome files — or be whole-genome multi-sample, or
any mix. Work is partitioned into (sample x chromosome) units run in parallel, then the
per-chromosome truth slices are assembled into each sample's truth VCF and one combined
`truth.tsv.gz`.
"""
import gzip
import math
from collections import OrderedDict
from importlib import resources
from pathlib import Path

from pathos.multiprocessing import ProcessPool

from .utils import (
    is_gzip,
    load_dict,
    load_list,
    init_logging,
    log_start,
    log_info,
    execute,
)
from .vcf import Vcf, list_samples, concat

AUTOSOMES = {f'chr{i}' for i in range(1, 23)}
MITO = {'chrM', 'chrMT'}
SEX = {'chrX', 'chrY'}

PLOIDY_ASSEMBLIES = ('hg38', 'hg19')


def load_ploidy(assembly):
    """Read the chr-prefixed `bcftools +fixploidy` ploidy table (CHROM FROM TO SEX
    PLOIDY) for a genome assembly from package resources. Unlisted regions default
    to ploidy 2, so autosomes and the PAR stay diploid; only male non-PAR chrX /
    chrY (MSY) and chrM are haploid, and female chrY is absent (ploidy 0)."""
    if assembly not in PLOIDY_ASSEMBLIES:
        raise ValueError(
            f'unsupported assembly {assembly!r}; expected one of {PLOIDY_ASSEMBLIES}')
    return (resources.files('genomics') / 'resources' / f'ploidy-{assembly}.txt').read_text()

def export_snv_truth(
    coordinates_file: Path,
    vcf_files: list,
    autosomes_depths_file: Path,
    sex_depths_file: Path,
    samples_file: Path,
    genders_file: Path,
    genome_file: Path,
    output_dir: Path,
    n_threads: int = 1,
    min_depth: int = 2,
    assembly: str = 'hg38',
    prod: bool = True,
):
    output_dir.mkdir(parents=True, exist_ok=True)
    init_logging(output_dir / 'snv_truth.log')

    info = OrderedDict()
    info['backbone-file'] = coordinates_file
    info['n_vcf_files'] = len(vcf_files)
    info['samples-file'] = samples_file
    info['genders-file'] = genders_file
    info['autosomes-depths-file'] = autosomes_depths_file
    info['sex-depths-file'] = sex_depths_file
    info['genome-file'] = genome_file
    info['output-dir'] = output_dir
    info['min-depth'] = min_depth
    info['assembly'] = assembly
    info['n-threads'] = n_threads
    log_start(banner='SNV Truth Creation', info=info)

    tmp_dir = output_dir / 'tmp'
    tmp_dir.mkdir(parents=True, exist_ok=True)

    log_info('load backbone coordinates')
    coordinates = load_coordinates(coordinates_file)
    families = group_families(coordinates)
    # Chromosomes to process come from the backbone; families are single-chromosome
    # (left-alignment never crosses contigs) and a cross-contig family is an error.
    families_by_chrom = partition_families_by_chrom(families)
    chroms = list(families_by_chrom.keys())

    log_info('load samples and genders')
    samples = load_list(samples_file)
    sample2gender = load_dict(genders_file)

    log_info('validate and inventory input vcfs')
    # Inputs must be bgzip-compressed AND indexed (fail fast; never mutate them).
    inventory = build_inventory(vcf_files)
    sources = resolve_sources(inventory, samples, chroms)

    # Depth is consulted only for the few sites that fall through to the depth
    # decision, so don't pre-fetch the panel; validate the tables up front (present +
    # tabix-indexed, no region reads) and let each unit query on demand.
    _require_indexed_depth(autosomes_depths_file)
    _require_indexed_depth(sex_depths_file)

    ploidy_file = output_dir / f'{assembly}.ploidy'
    ploidy_file.write_text(load_ploidy(assembly))

    # One header per sample (from the reference .fai) shared across that sample's
    # per-chromosome slices, so slices concat/sort cleanly and a no-source
    # chromosome can still be emitted with a valid header.
    headers = {}
    for sample in samples:
        sample_dir = tmp_dir / sample
        sample_dir.mkdir(parents=True, exist_ok=True)
        headers[sample] = build_sample_header(sample, genome_file, sample_dir / 'header.vcf')

    log_info('build (sample x chrom) truth')
    jobs = [
        {
            'sample': sample,
            'chrom': chrom,
            'sources': sorted(sources[(sample, chrom)], key=str),
            'families_chrom': families_by_chrom[chrom],
            'gender': sample2gender[sample],
            'genome_file': genome_file,
            'autosomes_depths_file': autosomes_depths_file,
            'sex_depths_file': sex_depths_file,
            'min_depth': min_depth,
            'ploidy_file': ploidy_file,
            'header_file': headers[sample],
            'tmp_dir': tmp_dir,
        }
        for sample in samples
        for chrom in chroms
    ]

    unit_out = {}
    if n_threads > 1:
        with ProcessPool(n_threads) as pool:
            for sample, chrom, vcf_path, body_path in pool.uimap(process_unit, jobs):
                unit_out[(sample, chrom)] = (Path(vcf_path), Path(body_path))
    else:
        for job in jobs:
            sample, chrom, vcf_path, body_path = process_unit(job)
            unit_out[(sample, chrom)] = (Path(vcf_path), Path(body_path))

    log_info('assemble per-sample truth vcfs')
    for sample in samples:
        parts = [unit_out[(sample, chrom)][0] for chrom in chroms]
        concat(parts, output_dir / f'{sample}.truth.vcf.bgz', tmp_dir / sample, preprocess=False)

    log_info('combine truth into one tsv')
    combine_tsv(samples, chroms, unit_out, output_dir / 'truth.tsv.gz')

    log_info('done')


def _index_path(vcf_file):
    for ext in ('.csi', '.tbi'):
        candidate = Path(f'{vcf_file}{ext}')
        if candidate.exists():
            return candidate
    return None


def _require_indexed_bgz(vcf_file):
    """Inputs must be bgzip-compressed AND indexed; fail fast otherwise. Inputs are
    never auto-indexed and never mutated."""
    vcf_file = Path(vcf_file)
    if not vcf_file.exists():
        raise ValueError(f'input vcf not found: {vcf_file}')
    if not is_gzip(vcf_file):
        raise ValueError(f'input vcf is not bgzip-compressed: {vcf_file}')
    if _index_path(vcf_file) is None:
        raise ValueError(
            f'input vcf is not indexed: {vcf_file}; index it with '
            f'`bcftools index {vcf_file}`')


def _vcf_contigs(vcf_file):
    """Contigs present in the index (works for both .csi and .tbi)."""
    out = execute(f'bcftools index -s {vcf_file}', pipe=True)
    return {line.split('\t')[0] for line in out if line.strip()}


def build_inventory(vcf_files):
    """For each input VCF (validated bgz+indexed), its sample set and contig set."""
    inventory = []
    for vcf_file in sorted(vcf_files, key=str):
        _require_indexed_bgz(vcf_file)
        inventory.append({
            'path': Path(vcf_file),
            'samples': set(list_samples(vcf_file)),
            'contigs': _vcf_contigs(vcf_file),
        })
    return inventory


def resolve_sources(inventory, samples, chroms):
    """Map each (sample, chrom) to the input paths holding both. A target sample in
    no input is an error; an empty source list for a (sample, chrom) is legal (that
    unit emits depth-based fills)."""
    absent = [s for s in samples if not any(s in e['samples'] for e in inventory)]
    if absent:
        raise ValueError(f'target samples not found in any input vcf: {sorted(absent)}')
    sources = {}
    for sample in samples:
        for chrom in chroms:
            sources[(sample, chrom)] = [
                e['path'] for e in inventory
                if sample in e['samples'] and chrom in e['contigs']
            ]
    return sources


def partition_families_by_chrom(families):
    """Split families by chromosome, preserving backbone (first-seen) order. A family
    whose members span more than one contig is an error (never silently truncated)."""
    by_chrom = OrderedDict()
    for family_id, members in families.items():
        member_chroms = {member['chrom'] for member in members}
        if len(member_chroms) != 1:
            raise ValueError(
                f'backbone family {family_id!r} spans multiple contigs: '
                f'{sorted(member_chroms)}')
        chrom = next(iter(member_chroms))
        by_chrom.setdefault(chrom, OrderedDict())[family_id] = members
    return by_chrom


def _require_indexed_depth(depths_file):
    """A depth table must be bgzip-compressed AND tabix-indexed. Validated up front
    (presence + index only; no region reads), so a missing/unindexed table is
    reported before processing rather than deep inside a worker."""
    depths_file = Path(depths_file)
    if not depths_file.exists():
        raise ValueError(f'depth table not found: {depths_file}')
    if not is_gzip(depths_file):
        raise ValueError(f'depth table is not bgzip-compressed: {depths_file}')
    if _index_path(depths_file) is None:
        raise ValueError(
            f'depth table is not indexed: {depths_file}; index it with '
            f'`tabix -s 1 -b 2 -e 2 -S 1 {depths_file}`')


class DepthProvider:
    """Per-unit, memoized, on-demand depth lookup. A site's depth is fetched from the
    chromosome-appropriate tabix-indexed table only when a genotype decision needs it
    (unmatched family member with no observed 0/0 / ./.), and cached by (chrom, pos)
    so the same position is never queried twice. Process-local; never shared across
    parallel units. `n_queries` counts the tabix region queries issued (for tests)."""

    def __init__(self, autosomes_depths_file, sex_depths_file):
        self.autosomes_depths_file = autosomes_depths_file
        self.sex_depths_file = sex_depths_file
        self._autosomes = {}   # (chrom, pos) -> float | None
        self._sex = {}         # (chrom, pos) -> {'male', 'female'} | None
        self.n_queries = 0

    def autosomes(self, chrom, pos):
        """depth_mean for (chrom, pos), or None if the site is absent."""
        key = (chrom, pos)
        if key not in self._autosomes:
            row = self._query(self.autosomes_depths_file, chrom, pos)
            self._autosomes[key] = _to_float(row[2]) if row else None
        return self._autosomes[key]

    def sex(self, chrom, pos):
        """{'male', 'female'} means for (chrom, pos), or None if the site is absent."""
        key = (chrom, pos)
        if key not in self._sex:
            row = self._query(self.sex_depths_file, chrom, pos)
            self._sex[key] = (
                {'male': _to_float(row[2]), 'female': _to_float(row[4])} if row else None)
        return self._sex[key]

    def _query(self, depths_file, chrom, pos):
        """One single-position tabix region query; returns the split row or None."""
        self.n_queries += 1
        for line in execute(f'tabix {depths_file} {chrom}:{pos}-{pos}', pipe=True):
            if line.strip():
                return line.split('\t')
        return None


def build_sample_header(sample, genome_file, out_file):
    """A minimal, complete single-sample VCF header from the reference .fai: every
    contig (so per-chromosome slices concat/sort cleanly), a GT FORMAT line, and the
    sample column. Reused as the header of every per-chromosome truth slice for the
    sample so no-source chromosomes still get a valid header."""
    fai = Path(f'{genome_file}.fai')
    if not fai.exists():
        raise ValueError(
            f'reference {genome_file} is not indexed: {fai} not found; '
            f'index it with `samtools faidx {genome_file}`')
    lines = ['##fileformat=VCFv4.2']
    with fai.open('rt') as fh:
        for line in fh:
            name, length = line.split('\t')[:2]
            lines.append(f'##contig=<ID={name},length={length}>')
    lines.append('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">')
    lines.append('\t'.join(
        ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', sample]))
    out_file.write_text('\n'.join(lines) + '\n')
    return out_file


def _extract_slice(src, sample, chrom, out_file):
    """bcftools view of one sample restricted to one chromosome -> bgz slice."""
    execute(f'bcftools view -r {chrom} -s {sample} -O z -o {out_file} {src}')
    return out_file


def process_unit(job):
    """Build the truth for one (sample, chromosome): region-extract the sample's
    slice for the chromosome (concat if several sources; none -> all depth fills),
    run the unchanged prepare -> match -> fill -> ploidy pipeline over that
    chromosome's families, and emit a per-chromosome truth VCF slice + TSV body."""
    sample = job['sample']
    chrom = job['chrom']
    sources = job['sources']
    genome_file = job['genome_file']

    unit_tmp = job['tmp_dir'] / sample / chrom
    unit_tmp.mkdir(parents=True, exist_ok=True)

    if sources:
        if len(sources) == 1:
            slice_vcf = _extract_slice(sources[0], sample, chrom, unit_tmp / 'slice.vcf.bgz')
        else:
            parts = [
                _extract_slice(src, sample, chrom, unit_tmp / f'part{i}.vcf.bgz')
                for i, src in enumerate(sources)
            ]
            slice_vcf = concat(parts, unit_tmp / 'slice.vcf.bgz', unit_tmp).filepath
        prepared = prepare_sample(slice_vcf, sample, genome_file, unit_tmp)
        lookup, nonvariant = build_lookup(prepared)
    else:
        lookup, nonvariant = {}, {}

    matched_lines, fill_lines = call_families(
        families=job['families_chrom'],
        lookup=lookup,
        nonvariant=nonvariant,
        gender=job['gender'],
        depth_provider=DepthProvider(job['autosomes_depths_file'], job['sex_depths_file']),
        min_depth=job['min_depth'],
    )

    fixed_fill_lines = fixploidy_lines(
        fill_lines=fill_lines,
        header_file=job['header_file'],
        sample_name=sample,
        gender=job['gender'],
        ploidy_file=job['ploidy_file'],
        tmp_base=unit_tmp,
    )

    truth_txt = unit_tmp / 'truth.vcf'
    with open(job['header_file'], 'rt') as hfh, truth_txt.open('wt') as bfh:
        for line in hfh:
            bfh.write(line)
        for line in matched_lines:
            bfh.write(line + '\n')
        for line in fixed_fill_lines:
            bfh.write(line + '\n')

    vcf = Vcf(truth_txt, unit_tmp).bgzip().sort().index()
    body_file = unit_tmp / 'body.tsv'
    write_body_tsv(vcf.filepath, body_file)

    return sample, chrom, str(vcf.filepath), str(body_file)


def write_body_tsv(vcf_file, output_file):
    """bcftools query -> headerless TSV body: chrom pos id ref alt tgt. The header
    and sample_name column are added once when per-unit bodies are combined."""
    cmd = (''
           f"bcftools query"
           f"      -f '%CHROM\\t%POS\\t%ID\\t%REF\\t%ALT\\t[%TGT]\\n'"
           f"      {vcf_file}"
           f"      > {output_file}"
           '')
    execute(cmd)
    return output_file


def combine_tsv(samples, chroms, unit_out, output_file):
    """Concatenate the per-(sample, chrom) TSV bodies into one gzip TSV, writing the
    header once and prepending each row with its sample_name. Samples are emitted in
    sorted order and chromosomes in backbone contig order, so the combined output is
    deterministic regardless of processing order."""
    with gzip.open(output_file, 'wt') as out:
        out.write('sample_name\tchrom\tpos\tid\tref\talt\ttgt\n')
        for sample in sorted(samples):
            for chrom in chroms:
                body_file = unit_out[(sample, chrom)][1]
                with body_file.open('rt') as fh:
                    for line in fh:
                        out.write(f'{sample}\t{line}')
    return output_file


def prepare_sample(vcf_file, sample, genome_file, tmp_base):
    """Extract one sample, trim unused ALTs, left-align + parsimonious normalize
    (no multiallelic split). Records with no called ALT (observed 0/0 / ./.) are
    kept — with ALT='.' — so their genotype can be honored during filling."""
    vcf = (
        Vcf(vcf_file, tmp_base)
        .subset_samples({sample})
        .keep_format(fields=['GT'])
        .drop_info()
        .trim_alts()
        .normalize(genome_file)
        .uppercase()
        .index()
    )
    return vcf.filepath


def build_lookup(prepared_vcf):
    """From the single-sample VCF: the exact-match variant lookup
    (chrom, pos, ref, ALT-set) -> GT, plus a non-variant lookup (chrom, pos) ->
    'homref'|'nocall' built from the retained ALT='.' records (observed 0/0 / ./.)."""
    lookup = dict()
    nonvariant = dict()
    with gzip.open(prepared_vcf, 'rt') as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            items = line.rstrip('\n').split('\t')
            chrom, pos, _id, ref, alt = items[0], items[1], items[2], items[3], items[4]
            gt = items[9].split(':')[0]
            if alt == '.':
                alleles = gt.replace('|', '/').split('/')
                nonvariant[(chrom, int(pos))] = 'nocall' if '.' in alleles else 'homref'
            else:
                lookup[match_key(chrom, pos, ref, alt)] = gt
    return lookup, nonvariant


def call_families(families, lookup, nonvariant, gender, depth_provider, min_depth):
    """Per family: matched member(s) with the sample GT (siblings dropped), else one
    fill record per member. A fill honors the sample's observed 0/0 / ./. at the
    member's position; only a site with no sample record falls back to the depth
    decision. All fill genotypes are ploidy-corrected downstream."""
    matched_lines = list()
    fill_lines = list()

    for members in families.values():
        matched = [
            (member, lookup[match_key(member['chrom'], member['pos'], member['ref'], member['alt'])])
            for member in members
            if match_key(member['chrom'], member['pos'], member['ref'], member['alt']) in lookup
        ]

        if matched:
            for member, gt in matched:
                matched_lines.append(vcf_line(member, gt))
        else:
            # Members of a family can sit at different (chrom, pos) after axiom
            # backbone re-normalizes/left-aligns the expanded records, so each
            # member is resolved at its own position.
            for member in members:
                observed = nonvariant.get((member['chrom'], member['pos']))
                if observed == 'homref':
                    gt = '0/0'
                elif observed == 'nocall':
                    gt = './.'
                else:
                    depth = depth_for(member['chrom'], member['pos'], gender, depth_provider)
                    gt = '0/0' if is_homref(depth, min_depth) else './.'
                fill_lines.append(vcf_line(member, gt))

    return matched_lines, fill_lines


def vcf_line(member, gt):
    return '\t'.join([
        member['chrom'],
        str(member['pos']),
        member['id'],
        member['ref'],
        member['alt'],
        '.', '.', '.', 'GT', gt,
    ])


def fixploidy_lines(fill_lines, header_file, sample_name, gender, ploidy_file, tmp_base):
    """Apply GRCh38 gender-aware ploidy to the fill records ONLY (never to matched
    real genotypes). Returns ploidy-corrected data lines."""
    if not fill_lines:
        return []

    fills_file = tmp_base / 'fills.vcf'
    with header_file.open('rt') as hfh, fills_file.open('wt') as ffh:
        for line in hfh:
            ffh.write(line)
        for line in fill_lines:
            ffh.write(line + '\n')

    sex_file = tmp_base / 'sex.txt'
    sex = 'M' if gender == 'male' else 'F'
    sex_file.write_text(f'{sample_name}\t{sex}\n')

    cmd = (''
           f'bcftools +fixploidy {fills_file} --'
           f'      -p {ploidy_file}'
           f'      -s {sex_file}'
           f'      -t GT'
           '')
    out = execute(cmd, pipe=True)
    return [line for line in out if not line.startswith('#')]


def load_coordinates(coordinates_file):
    """Read the backbone VCF (bgzip or plain) into coordinate records."""
    opener = gzip.open if is_gzip(coordinates_file) else open
    rows = list()
    with opener(coordinates_file, 'rt') as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            if not line.strip():
                continue
            items = line.rstrip('\n').split('\t')
            rows.append({
                'chrom': items[0],
                'pos': int(items[1]),
                'id': items[2],
                'ref': items[3].upper(),
                'alt': items[4].upper(),
            })
    return rows


def group_families(coordinates):
    """Group coordinate records into SNP families by the exact ID string."""
    families = OrderedDict()
    for row in coordinates:
        families.setdefault(row['id'], []).append(row)
    return families


def depth_for(chrom, pos, sex, depth_provider):
    """Select the depth for a site via the on-demand provider: contig class picks the
    table; sex picks the sex-chromosome column. Missing site -> None (treated as
    below threshold)."""
    pos = int(pos)
    if chrom in AUTOSOMES or chrom in MITO:
        return depth_provider.autosomes(chrom, pos)
    if chrom in SEX:
        record = depth_provider.sex(chrom, pos)
        if record is None:
            return None
        return record['male'] if sex == 'male' else record['female']
    return None


def match_key(chrom, pos, ref, alt):
    """Exact-match key: (chrom, pos, ref, set-of-ALT-alleles), case-insensitive."""
    return (chrom, int(pos), ref.upper(), frozenset(a.upper() for a in alt.split(',')))


def is_homref(depth, min_depth):
    if depth is None:
        return False
    if isinstance(depth, float) and math.isnan(depth):
        return False
    return depth >= min_depth


def _to_float(value):
    value = value.strip()
    if value == '':
        return None
    try:
        return float(value)
    except ValueError:
        return math.nan
