"""Backbone-driven exact-match SNV truth calling.

`genomics snv-truth` consumes an `axiom backbone` output VCF as `--backbone-file`
and builds a per-sample truth VCF (+ TSV) by exact `(chrom, pos, ref, ALT-set)`
matching against a normalized single-sample VCF. Backbone records that share an `ID`
form a SNP family: if one member matches the sample it is emitted with the sample's
genotype and its siblings are dropped; if none matches, every member is emitted with a
depth- and gender-aware homozygous-reference or no-call genotype.

Each input VCF (bgzip-compressed and indexed) must contain exactly one chromosome. Each
chromosome is decoded once — the target samples present are extracted in a single pass
and split per sample — then each sample's per-chromosome pieces are concatenated into
one per-sample VCF and processed in parallel across samples into that sample's truth VCF
and one combined `truth.tsv.gz`.
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

    log_info('load samples and genders')
    samples = load_list(samples_file)
    sample2gender = load_dict(genders_file)

    log_info('validate inputs (one chromosome per vcf, bgzip + indexed)')
    # Each input VCF must hold exactly one chromosome; fail fast and never mutate.
    inputs = build_inputs(vcf_files)
    target = set(samples)
    all_input_samples = set().union(*(inp['samples'] for inp in inputs)) if inputs else set()
    absent = [s for s in samples if s not in all_input_samples]
    if absent:
        raise ValueError(f'target samples not found in any input vcf: {sorted(absent)}')

    # Depth is consulted only for the few sites that fall through to the depth
    # decision, so don't pre-fetch the panel; validate the tables up front (present +
    # tabix-indexed, no region reads) and let each sample query on demand.
    _require_indexed_depth(autosomes_depths_file)
    _require_indexed_depth(sex_depths_file)

    ploidy_file = output_dir / f'{assembly}.ploidy'
    ploidy_file.write_text(load_ploidy(assembly))

    # One header per sample (from the reference .fai) — every contig + GT + the
    # sample column — used to wrap that sample's truth records.
    headers = {}
    for sample in samples:
        sample_dir = tmp_dir / sample
        sample_dir.mkdir(parents=True, exist_ok=True)
        headers[sample] = build_sample_header(sample, genome_file, sample_dir / 'header.vcf')

    # Extract each chromosome ONCE for all target samples present, then split per
    # sample. Each input file is decoded exactly once (independent of sample count).
    log_info('extract each chromosome once and split by sample')
    extract_dir = tmp_dir / 'extract'
    extract_dir.mkdir(parents=True, exist_ok=True)
    extract_jobs = [
        {
            'path': inp['path'],
            'chrom': inp['chrom'],
            'targets_present': sorted(target & inp['samples']),
            'tmp_dir': extract_dir,
        }
        for inp in inputs
    ]
    sample_pieces = {sample: {} for sample in samples}
    if n_threads > 1:
        with ProcessPool(n_threads) as pool:
            for chrom, pieces in pool.uimap(extract_and_split, extract_jobs):
                for sample, piece in pieces.items():
                    sample_pieces[sample][chrom] = piece
    else:
        for job in extract_jobs:
            chrom, pieces = extract_and_split(job)
            for sample, piece in pieces.items():
                sample_pieces[sample][chrom] = piece

    # Combine each sample's per-chromosome pieces and build its truth, in parallel
    # across samples.
    log_info('build per-sample truth (parallel over samples)')
    process_jobs = [
        {
            'sample': sample,
            'pieces': [sample_pieces[sample][c] for c in sorted(sample_pieces[sample])],
            'gender': sample2gender[sample],
            'families': families,
            'genome_file': genome_file,
            'autosomes_depths_file': autosomes_depths_file,
            'sex_depths_file': sex_depths_file,
            'min_depth': min_depth,
            'ploidy_file': ploidy_file,
            'header_file': headers[sample],
            'output_dir': output_dir,
            'tmp_dir': tmp_dir,
        }
        for sample in samples
    ]
    unit_out = {}
    if n_threads > 1:
        with ProcessPool(n_threads) as pool:
            for sample, body_path in pool.uimap(process_sample, process_jobs):
                unit_out[sample] = Path(body_path)
    else:
        for job in process_jobs:
            sample, body_path = process_sample(job)
            unit_out[sample] = Path(body_path)

    log_info('combine truth into one tsv')
    combine_tsv(samples, unit_out, output_dir / 'truth.tsv.gz')

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


def _require_single_chrom(vcf_file):
    """Return the single chromosome an input VCF contains; error if it spans more
    than one (or none). Enforces the one-chromosome-per-input contract from the
    index alone (no record scan)."""
    contigs = sorted(_vcf_contigs(vcf_file))
    if len(contigs) != 1:
        raise ValueError(
            f'input vcf must contain exactly one chromosome: {vcf_file} spans {contigs}')
    return contigs[0]


def build_inputs(vcf_files):
    """Validate each input (bgzip + indexed, exactly one chromosome) and record its
    chromosome and sample set. Inputs are read-only and never mutated."""
    inputs = []
    for vcf_file in sorted(vcf_files, key=str):
        _require_indexed_bgz(vcf_file)
        chrom = _require_single_chrom(vcf_file)
        inputs.append({
            'path': Path(vcf_file),
            'chrom': chrom,
            'samples': set(list_samples(vcf_file)),
        })
    return inputs


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


def extract_and_split(job):
    """Extract the target samples present in one chromosome's input in a single pass
    and split them into per-sample pieces. The input file is decoded exactly once.
    Returns (chrom, {sample: piece_path})."""
    chrom = job['chrom']
    targets_present = job['targets_present']
    work = job['tmp_dir'] / chrom
    work.mkdir(parents=True, exist_ok=True)
    if not targets_present:
        return chrom, {}
    # One streamed pass: strip INFO + all non-GT FORMAT, then select the target
    # samples, then split per sample. bcftools +split copies INFO/FORMAT into every
    # per-sample file (~88x write amplification) and the truth pipeline discards
    # both downstream, so stripping before selection keeps the read single-pass and
    # is output-neutral.
    pieces = Vcf(job['path'], work).strip_select_split(set(targets_present))
    return chrom, {sample: Path(path) for sample, path in pieces.items()}


def process_sample(job):
    """Build the truth for one sample: concatenate its per-chromosome pieces into one
    whole-genome VCF, then run the unchanged prepare -> match -> fill -> ploidy
    pipeline over the full backbone. Emits {sample}.truth.vcf.bgz + a TSV body."""
    sample = job['sample']
    genome_file = job['genome_file']

    tmp_base = job['tmp_dir'] / sample
    tmp_base.mkdir(parents=True, exist_ok=True)

    # Pieces are one-per-chromosome and disjoint, so a plain concat (no
    # --allow-overlaps) suffices and needs no per-piece index; the trailing
    # sort inside concat() still coordinate-sorts the combined result.
    combined = concat(
        job['pieces'], tmp_base / 'combined.vcf.bgz', tmp_base,
        preprocess=False, allow_overlaps=False).filepath
    prepared = prepare_sample(combined, genome_file, tmp_base)
    lookup, nonvariant = build_lookup(prepared)

    matched_lines, fill_lines = call_families(
        families=job['families'],
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
        tmp_base=tmp_base,
    )

    truth_txt = tmp_base / 'truth.vcf'
    with open(job['header_file'], 'rt') as hfh, truth_txt.open('wt') as bfh:
        for line in hfh:
            bfh.write(line)
        for line in matched_lines:
            bfh.write(line + '\n')
        for line in fixed_fill_lines:
            bfh.write(line + '\n')

    vcf = Vcf(truth_txt, tmp_base).sort().index()
    out_vcf = job['output_dir'] / f'{sample}.truth.vcf.bgz'
    vcf.move_to(out_vcf)

    body_file = tmp_base / 'body.tsv'
    write_body_tsv(out_vcf, body_file)

    return sample, str(body_file)


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


def combine_tsv(samples, unit_out, output_file):
    """Concatenate the per-sample TSV bodies into one gzip TSV, writing the header
    once and prepending each row with its sample_name. Samples are emitted in sorted
    order; rows within a sample are already contig+position sorted (the per-sample
    truth VCF is sorted), so the combined output is deterministic."""
    with gzip.open(output_file, 'wt') as out:
        out.write('sample_name\tchrom\tpos\tid\tref\talt\ttgt\n')
        for sample in sorted(samples):
            body_file = unit_out[sample]
            with body_file.open('rt') as fh:
                for line in fh:
                    out.write(f'{sample}\t{line}')
    return output_file


def prepare_sample(vcf_file, genome_file, tmp_base):
    """Trim unused ALTs, left-align + parsimonious normalize (no multiallelic split).
    Records with no called ALT (observed 0/0 / ./.) are kept — with ALT='.' — so their
    genotype can be honored during filling.

    Preparation is trim+normalize only: no uppercasing and no indexing. The prepared
    VCF's case and index are never used — match_key is case-insensitive, and the truth
    VCF's REF/ALT come from the backbone, so a case-folding/indexing pass here would be
    dead work.

    The input is already single-sample, GT-only, and INFO-stripped (each piece comes
    from `bcftools +split` after the lean extraction, then a same-sample concat), so
    sample subsetting, FORMAT-keeping, and INFO-dropping are unnecessary here. This
    relies on that extraction contract; if it changes, restore those steps."""
    vcf = Vcf(vcf_file, tmp_base).trim_alts_normalize(genome_file)
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
