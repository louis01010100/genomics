"""Backbone-driven exact-match SNV truth calling.

`genomics snv-truth` consumes an `axiom backbone` output VCF as `--backbone-file`
and builds a per-sample truth VCF (+ TSV) by exact `(chrom, pos, ref, ALT-set)`
matching against a normalized single-sample VCF. Backbone records that share an `ID`
form a SNP family: if one member matches the sample it is emitted with the sample's
genotype and its siblings are dropped; if none matches, every member is emitted with a
depth- and gender-aware homozygous-reference or no-call genotype.
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
    copy_vcf_header,
    execute,
)
from .vcf import Vcf, list_samples

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

    log_info('load samples, genders, depths')
    samples = load_list(samples_file)
    sample2gender = load_dict(genders_file)
    # Only the backbone positions are ever looked up, so random-search the
    # tabix-indexed depth tables for those alone instead of loading the whole
    # (genome-sized) tables into memory and pickling them into every worker.
    coord_keys = sorted({(row['chrom'], row['pos']) for row in coordinates})
    autosomes_depths = load_autosomes_depths(autosomes_depths_file, coord_keys, tmp_dir)
    sex_depths = load_sex_depths(sex_depths_file, coord_keys, tmp_dir)

    log_info('locate samples in input vcfs')
    sample2vcf = locate_samples(vcf_files, samples)

    ploidy_file = output_dir / f'{assembly}.ploidy'
    ploidy_file.write_text(load_ploidy(assembly))

    log_info('build per-sample truth')
    jobs = [
        {
            'sample': sample,
            'vcf_file': sample2vcf[sample],
            'gender': sample2gender[sample],
            'families': families,
            'genome_file': genome_file,
            'autosomes_depths': autosomes_depths,
            'sex_depths': sex_depths,
            'min_depth': min_depth,
            'ploidy_file': ploidy_file,
            'output_dir': output_dir,
            'tmp_dir': tmp_dir,
        }
        for sample in samples
    ]

    if n_threads > 1:
        with ProcessPool(n_threads) as pool:
            for _ in pool.uimap(process_sample, jobs):
                pass
    else:
        for job in jobs:
            process_sample(job)

    log_info('combine per-sample truth into one tsv')
    combine_tsv(samples, tmp_dir, output_dir / 'truth.tsv.gz')

    log_info('done')


def process_sample(job):
    sample = job['sample']
    vcf_file = job['vcf_file']
    gender = job['gender']
    families = job['families']
    genome_file = job['genome_file']
    autosomes_depths = job['autosomes_depths']
    sex_depths = job['sex_depths']
    min_depth = job['min_depth']
    ploidy_file = job['ploidy_file']
    output_dir = job['output_dir']

    tmp_base = job['tmp_dir'] / sample
    tmp_base.mkdir(parents=True, exist_ok=True)

    prepared = prepare_sample(vcf_file, sample, genome_file, tmp_base)
    prepared_name = list_samples(prepared)[0]

    lookup, nonvariant = build_lookup(prepared)

    matched_lines, fill_lines = call_families(
        families=families,
        lookup=lookup,
        nonvariant=nonvariant,
        gender=gender,
        autosomes_depths=autosomes_depths,
        sex_depths=sex_depths,
        min_depth=min_depth,
    )

    header_file = tmp_base / 'header.vcf'
    copy_vcf_header(prepared, header_file)

    fixed_fill_lines = fixploidy_lines(
        fill_lines=fill_lines,
        header_file=header_file,
        sample_name=prepared_name,
        gender=gender,
        ploidy_file=ploidy_file,
        tmp_base=tmp_base,
    )

    body_file = tmp_base / 'truth.vcf'
    with header_file.open('rt') as hfh, body_file.open('wt') as bfh:
        for line in hfh:
            bfh.write(line)
        for line in matched_lines:
            bfh.write(line + '\n')
        for line in fixed_fill_lines:
            bfh.write(line + '\n')

    vcf = Vcf(body_file, tmp_base).bgzip().sort().index()
    out_vcf = output_dir / f'{sample}.truth.vcf.bgz'
    vcf.move_to(out_vcf)

    write_body_tsv(out_vcf, tmp_base / 'body.tsv')

    return sample


def write_body_tsv(vcf_file, output_file):
    """bcftools query -> headerless TSV body: chrom pos id ref alt tgt. The header
    and sample_name column are added once when per-sample bodies are combined."""
    cmd = (''
           f"bcftools query"
           f"      -f '%CHROM\\t%POS\\t%ID\\t%REF\\t%ALT\\t[%TGT]\\n'"
           f"      {vcf_file}"
           f"      > {output_file}"
           '')
    execute(cmd)
    return output_file


def combine_tsv(samples, tmp_dir, output_file):
    """Concatenate the per-sample TSV bodies into one gzip TSV, writing the header
    once and prepending each row with its sample_name. Samples are emitted in sorted
    order so the combined output is deterministic regardless of processing order."""
    with gzip.open(output_file, 'wt') as out:
        out.write('sample_name\tchrom\tpos\tid\tref\talt\ttgt\n')
        for sample in sorted(samples):
            body_file = tmp_dir / sample / 'body.tsv'
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


def call_families(families, lookup, nonvariant, gender, autosomes_depths, sex_depths, min_depth):
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
                    depth = depth_for(member['chrom'], member['pos'], gender, autosomes_depths, sex_depths)
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


def locate_samples(vcf_files, samples):
    """Map each target sample to the single input VCF that contains it."""
    wanted = set(samples)
    mapping = dict()
    for vcf_file in vcf_files:
        present = set(list_samples(vcf_file))
        for sample in wanted & present:
            mapping[sample] = vcf_file
    missing = wanted - set(mapping)
    if missing:
        raise ValueError(f'target samples not found in any input vcf: {sorted(missing)}')
    return mapping


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


def _fetch_depths(depths_file, coord_keys, tmp_dir):
    """Stream only the backbone positions out of a tabix-indexed depth table.

    Writes the coordinates to a 1bp-region BED and runs `tabix -R`, so the amount
    read is bounded by the panel size rather than the genome. Yields the split
    columns of each matching row (tabix output carries no header)."""
    regions_file = tmp_dir / f'{depths_file.name}.regions.bed'
    with regions_file.open('wt') as fh:
        for chrom, pos in coord_keys:
            fh.write(f'{chrom}\t{pos - 1}\t{pos}\n')
    for line in execute(f'tabix -R {regions_file} {depths_file}', pipe=True):
        yield line.split('\t')


def load_autosomes_depths(depths_file, coord_keys, tmp_dir):
    """(chrom, pos) -> depth_mean for the backbone coordinates only, fetched via
    tabix. Producer schema (autosomes-depth): chrom, pos, depth_mean, n_samples."""
    depths = dict()
    for items in _fetch_depths(depths_file, coord_keys, tmp_dir):
        depths[(items[0], int(items[1]))] = _to_float(items[2])
    return depths


def load_sex_depths(depths_file, coord_keys, tmp_dir):
    """(chrom, pos) -> {'male', 'female'} for the backbone coordinates only, fetched
    via tabix. Producer schema (sex-depth): chrom, pos, mean_male, n_male,
    mean_female, n_female."""
    depths = dict()
    for items in _fetch_depths(depths_file, coord_keys, tmp_dir):
        depths[(items[0], int(items[1]))] = {
            'male': _to_float(items[2]),
            'female': _to_float(items[4]),
        }
    return depths


def depth_for(chrom, pos, sex, autosomes_depths, sex_depths):
    """Select the depth for a site: contig class picks the table; sex picks the
    sex-chromosome column. Missing key -> None (treated as below threshold)."""
    pos = int(pos)
    if chrom in AUTOSOMES or chrom in MITO:
        return autosomes_depths.get((chrom, pos))
    if chrom in SEX:
        record = sex_depths.get((chrom, pos))
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
