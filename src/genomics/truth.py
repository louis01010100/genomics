"""SNV-family-driven exact-match SNV truth calling.

`genomics snv-truth` consumes an `axiom snv-family` output TSV as `--snv-family-file`
and builds a per-sample truth VCF (+ TSV) by exact `(chrom, pos, ref, ALT-set)`
matching against a normalized single-sample VCF. snv-family records that share a family
id (`fmid`, `FM-<n>`) form a SNP family: if one member matches the sample it is emitted
with the sample's genotype and its siblings are dropped; if none matches, every member is
emitted with a depth- and gender-aware homozygous-reference or no-call genotype. The
family id (`fmid`) is the truth output's `id`.

Each input VCF (bgzip-compressed and indexed) is a single-sample, whole-genome VCF
(e.g. the output of `genomics samples`); its sole sample column names the sample. Each
per-sample VCF is processed in parallel across samples into that sample's truth VCF and
gzipped TSV, written together under `<output-dir>/samples/` (`<sample>.truth.vcf.bgz`,
its `.csi`, and `<sample>.tsv.gz`).
"""
import gzip
import math
import pickle
from collections import OrderedDict
from importlib import resources
from pathlib import Path

from pathos.multiprocessing import ProcessPool

from .utils import (
    is_gzip,
    load_dict,
    init_logging,
    log_start,
    log_info,
    execute,
)
from .vcf import Vcf, list_samples

AUTOSOMES = {f'chr{i}' for i in range(1, 23)}
MITO = {'chrM', 'chrMT'}
SEX = {'chrX', 'chrY'}

PLOIDY_ASSEMBLIES = ('hg38', 'hg19')

# The snv-family depths are pre-loaded once per run and pickled to a file; each worker
# loads that file at most once and caches it here, so the (large, whole-genome) depth
# tables are read exactly once and the ~family-sized maps are never pickled per job — only
# the small file path travels in the job. Robust across fork/spawn worker start methods.
_DEPTHS_CACHE = {}   # depths_file path -> (autosomes_depths, sex_depths)


def _shared_depths(depths_file):
    """Return the once-per-run pre-loaded depth maps, loading the pickled file at most once
    per worker process and caching it (via the live module, so the cache persists across a
    worker's samples regardless of how the pool serializes `process_sample`)."""
    import genomics.truth as _self
    key = str(depths_file)
    if key not in _self._DEPTHS_CACHE:
        with open(key, 'rb') as fh:
            _self._DEPTHS_CACHE[key] = pickle.load(fh)
    return _self._DEPTHS_CACHE[key]


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
    snv_family_file: Path,
    vcf_files: list,
    autosomes_depths_file: Path,
    sex_depths_file: Path,
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
    info['snv-family-file'] = snv_family_file
    info['n_vcf_files'] = len(vcf_files)
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

    log_info('load snv-family coordinates')
    coordinates = load_snv_family(snv_family_file)
    families = group_families(coordinates)

    log_info('load genders')
    sample2gender = load_dict(genders_file)

    log_info('validate inputs (one sample per vcf, bgzip + indexed)')
    # Each input VCF is a per-sample, whole-genome VCF (e.g. `genomics samples`
    # output); validate one sample per file and map sample -> path. Never mutate.
    sample2vcf = build_sample_inputs(vcf_files)
    samples = sorted(sample2vcf)

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

    # The set of (chrom, pos) matching ever consults — every snv-family member position.
    # Derived once and shared across workers; each per-sample lookup is scoped to it so it
    # is bounded by the snv-family, not the genome.
    family_positions = frozenset(
        (member['chrom'], member['pos'])
        for members in families.values()
        for member in members
    )

    # Pre-load the depths for the snv-family positions ONCE — the tables are panel-level
    # and sample-independent, so a single streaming pass per table replaces the per-site
    # tabix queries that would otherwise run tens of thousands of times per sample. Publish
    # the maps to a module global BEFORE the pool is created so forked workers inherit them
    # read-only (copy-on-write); they are never carried in the per-job payload.
    log_info('pre-load snv-family depths (once, shared across workers)')
    targets_file = write_depth_targets(family_positions, tmp_dir / 'depth.targets.tsv')
    depths = (
        load_depth_map(autosomes_depths_file, targets_file, is_sex=False),
        load_depth_map(sex_depths_file, targets_file, is_sex=True),
    )
    depths_file = tmp_dir / 'depths.pkl'
    with open(depths_file, 'wb') as fh:
        pickle.dump(depths, fh)

    # Each input is already a single-sample, whole-genome VCF (e.g. `genomics samples`
    # output), so process each one directly into its truth, in parallel across samples.
    # Per-sample outputs (truth VCF + gzipped TSV) are written together under samples/.
    (output_dir / 'samples').mkdir(parents=True, exist_ok=True)

    log_info('build per-sample truth (parallel over samples)')
    process_jobs = [
        {
            'sample': sample,
            'vcf': sample2vcf[sample],
            'gender': sample2gender[sample],
            'families': families,
            'family_positions': family_positions,
            'genome_file': genome_file,
            'depths_file': depths_file,
            'min_depth': min_depth,
            'ploidy_file': ploidy_file,
            'header_file': headers[sample],
            'output_dir': output_dir,
            'tmp_dir': tmp_dir,
        }
        for sample in samples
    ]
    if n_threads > 1:
        with ProcessPool(n_threads) as pool:
            for _sample in pool.uimap(process_sample, process_jobs):
                pass
    else:
        for job in process_jobs:
            process_sample(job)

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


def _require_single_sample(vcf_file):
    """Return the single sample an input VCF contains; error if it holds != 1.
    Each input is a per-sample VCF, so its sole sample column names the sample."""
    samples = list(list_samples(vcf_file))
    if len(samples) != 1:
        raise ValueError(
            f'input vcf must contain exactly one sample: {vcf_file} holds {sorted(samples)}')
    return samples[0]


def build_sample_inputs(vcf_files):
    """Validate each input (bgzip + indexed, exactly one sample) and map its sample
    to its path. Inputs are read-only and never mutated."""
    sample2vcf = {}
    for vcf_file in sorted(vcf_files, key=str):
        _require_indexed_bgz(vcf_file)
        sample = _require_single_sample(vcf_file)
        if sample in sample2vcf:
            raise ValueError(
                f'sample {sample!r} appears in more than one input vcf: '
                f'{sample2vcf[sample]} and {vcf_file}')
        sample2vcf[sample] = Path(vcf_file)
    return sample2vcf


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
    """Read-only depth lookup over pre-loaded maps. The depths for the snv-family
    positions are read once per run (see `load_depth_map`) into `(chrom, pos)`-keyed
    maps that are the same for every sample, so a lookup is an in-memory dict hit with
    NO per-site subprocess. A position absent from a map resolves to None (treated as
    below threshold), identical to a single-position tabix query that finds no row."""

    def __init__(self, autosomes_depths, sex_depths):
        self._autosomes = autosomes_depths   # (chrom, pos) -> float | None
        self._sex = sex_depths               # (chrom, pos) -> {'male', 'female'} | None

    def autosomes(self, chrom, pos):
        """depth_mean for (chrom, pos), or None if the site is absent."""
        return self._autosomes.get((chrom, pos))

    def sex(self, chrom, pos):
        """{'male', 'female'} means for (chrom, pos), or None if the site is absent."""
        return self._sex.get((chrom, pos))


def write_depth_targets(family_positions, regions_file):
    """Write the snv-family positions as a sorted BED regions file — one single-base
    interval per position (`chrom`, 0-based start `pos-1`, 1-based end `pos`) — the format
    `tabix -R` consumes to fetch only those sites from a depth table."""
    with open(regions_file, 'wt') as fh:
        for chrom, pos in sorted(family_positions):
            fh.write(f'{chrom}\t{pos - 1}\t{pos}\n')
    return regions_file

def load_depth_map(depths_file, regions_file, is_sex):
    """Read the depths for exactly the family positions from one depth table via a tabix
    index range read, returning a `(chrom, pos)` -> value map. `is_sex` selects the row
    shape: autosomes -> `depth_mean` (col 3) as float|None; sex -> `{'male','female'}`
    (cols 3 and 5) as float|None.

    `tabix -R <regions_file>` uses the depth table's index to fetch only the blocks holding
    the snv-family sites — the whole table is never scanned. Depth rows are single positions,
    so the single-base regions return exactly the rows at family positions; a POS-membership
    guard (against the regions' 1-based ends) drops any header/overlap row. Parsing matches
    the previous per-site query exactly (same `_to_float(col)` on the same columns), so the
    map value at every position equals what a single-position `tabix` query returned."""
    family = set()
    with open(regions_file, 'rt') as fh:
        for line in fh:
            if not line.strip():
                continue
            chrom, _start, end = line.rstrip('\n').split('\t')
            family.add((chrom, int(end)))

    depths = dict()
    for line in execute(f'tabix -R {regions_file} {depths_file}', pipe=True):
        if not line.strip():
            continue
        items = line.split('\t')
        key = (items[0], int(items[1]))
        if key not in family:
            continue
        if is_sex:
            depths[key] = {'male': _to_float(items[2]), 'female': _to_float(items[4])}
        else:
            depths[key] = _to_float(items[2])
    return depths


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


def process_sample(job):
    """Build the truth for one sample from its per-sample, whole-genome input VCF:
    run the unchanged prepare -> match -> fill -> ploidy pipeline over the full
    snv-family. Emits {sample}.truth.vcf.bgz + a TSV body."""
    sample = job['sample']
    genome_file = job['genome_file']

    tmp_base = job['tmp_dir'] / sample
    tmp_base.mkdir(parents=True, exist_ok=True)

    prepared = prepare_sample(job['vcf'], genome_file, tmp_base)
    lookup, nonvariant = build_lookup(prepared, job['family_positions'])

    matched_lines, fill_lines = call_families(
        families=job['families'],
        lookup=lookup,
        nonvariant=nonvariant,
        gender=job['gender'],
        depth_provider=DepthProvider(*_shared_depths(job['depths_file'])),
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

    samples_dir = job['output_dir'] / 'samples'
    vcf = Vcf(truth_txt, tmp_base).sort().index()
    out_vcf = samples_dir / f'{sample}.truth.vcf.bgz'
    vcf.move_to(out_vcf)

    # Per-sample gzipped TSV: header + the query body (chrom pos id ref alt tgt).
    body_file = tmp_base / 'body.tsv'
    write_body_tsv(out_vcf, body_file)
    with gzip.open(samples_dir / f'{sample}.tsv.gz', 'wt') as out, body_file.open('rt') as fh:
        out.write('chrom\tpos\tid\tref\talt\ttgt\n')
        for line in fh:
            out.write(line)

    return sample


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


def build_lookup(prepared_vcf, family_positions):
    """From the single-sample VCF, restricted to the snv-family positions: the exact-match
    variant lookup (chrom, pos, ref, ALT-set) -> GT, plus a non-variant lookup
    (chrom, pos) -> 'homref'|'nocall' from the retained ALT='.' records (observed 0/0 / ./.).

    Only the family positions are retained — a single streaming pass (`bcftools view -T`)
    filters the prepared VCF in C, and only records whose own POS is a family position are
    kept — so both dicts are bounded by the snv-family, not the genome. This is byte-identical,
    for every key matching can consult, to a whole-genome build: `call_families` only ever
    queries the dicts at snv-family member (chrom, pos), so every dropped key is one matching
    never reads.

    A streaming `-T` filter is used rather than an indexed `-R` range read: the family spans
    ~10^5-10^6 scattered positions, for which one sequential pass is far faster than that many
    random index seeks, and it needs no index (so `prepare_sample` stays index-free). The
    prepared VCF is read exactly once, so an index would never be amortized.

    `family_positions` is a set of (chrom, pos) covering every snv-family member."""
    prepared_vcf = Path(prepared_vcf)

    # Materialize the family positions as a sorted targets file (chrom, 1-based pos) for the
    # streaming `-T` filter (no index required).
    targets_file = prepared_vcf.parent / 'family.targets.tsv'
    with targets_file.open('wt') as fh:
        for chrom, pos in sorted(family_positions):
            fh.write(f'{chrom}\t{pos}\n')

    lookup = dict()
    nonvariant = dict()
    for line in execute(f'bcftools view -H -T {targets_file} {prepared_vcf}', pipe=True):
        if not line.strip():
            continue
        items = line.rstrip('\n').split('\t')
        chrom, pos, _id, ref, alt = items[0], items[1], items[2], items[3], items[4]
        # Select on the record's own POS: a record merely spanning a family position
        # (e.g. an indel returned by span-overlap) is not a family site.
        if (chrom, int(pos)) not in family_positions:
            continue
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


def load_snv_family(snv_family_file):
    """Read the (gzip) snv-family TSV into coordinate records, with each record's `id`
    set to its family id (`fmid`). The TSV has a header row and carries
    `chrom, pos, ref, alt, fmid, msid, asid, psid`; the needed columns are selected by
    NAME (robust to added/reordered columns) and `msid/asid/psid` are read past — they
    are not needed for truth calling. Grouping by `id` then groups by `fmid`."""
    opener = gzip.open if is_gzip(snv_family_file) else open
    with opener(snv_family_file, 'rt') as fh:
        header = fh.readline().rstrip('\n').split('\t')
        i_chrom = header.index('chrom')
        i_pos = header.index('pos')
        i_ref = header.index('ref')
        i_alt = header.index('alt')
        i_fmid = header.index('fmid')
        rows = list()
        for line in fh:
            if not line.strip():
                continue
            items = line.rstrip('\n').split('\t')
            rows.append({
                'chrom': items[i_chrom],
                'pos': int(items[i_pos]),
                'id': items[i_fmid],
                'ref': items[i_ref].upper(),
                'alt': items[i_alt].upper(),
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
