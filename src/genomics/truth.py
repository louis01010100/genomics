"""Backbone-driven exact-match SNV truth calling.

`genomics snv-truth` consumes an `axiom backbone` output VCF as `--coordinates-file`
and builds a per-sample truth VCF (+ TSV) by exact `(chrom, pos, ref, ALT-set)`
matching against a normalized single-sample VCF. Backbone records that share an `ID`
form a SNP family: if one member matches the sample it is emitted with the sample's
genotype and its siblings are dropped; if none matches, every member is emitted with a
depth- and gender-aware homozygous-reference or no-call genotype.
"""
import gzip
import math
import shutil
from collections import OrderedDict
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

# chr-prefixed GRCh38 ploidy for `bcftools +fixploidy` (CHROM FROM TO SEX PLOIDY).
# Unlisted regions default to ploidy 2, so autosomes and the PAR stay diploid. Male
# non-PAR chrX / chrY (MSY) and chrM are haploid; female chrY is absent (ploidy 0).
GRCH38_PLOIDY = (
    'chrX 2781480 155701382 M 1\n'
    'chrX 156030896 156040895 M 1\n'
    'chrY 1 10000 M 1\n'
    'chrY 2781480 56887902 M 1\n'
    'chrY 57217416 57227415 M 1\n'
    'chrY 1 57227415 F 0\n'
    'chrM 1 16569 M 1\n'
    'chrM 1 16569 F 1\n'
)

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
    prod: bool = True,
):
    output_dir.mkdir(parents=True, exist_ok=True)
    init_logging(output_dir / 'snv_truth.log')

    info = OrderedDict()
    info['coordinates-file'] = coordinates_file
    info['n_vcf_files'] = len(vcf_files)
    info['samples-file'] = samples_file
    info['genders-file'] = genders_file
    info['autosomes-depths-file'] = autosomes_depths_file
    info['sex-depths-file'] = sex_depths_file
    info['genome-file'] = genome_file
    info['output-dir'] = output_dir
    info['min-depth'] = min_depth
    info['n-threads'] = n_threads
    log_start(banner='SNV Truth Creation', info=info)

    log_info('load samples, genders, depths')
    samples = load_list(samples_file)
    sample2gender = load_dict(genders_file)
    autosomes_depths = load_autosomes_depths(autosomes_depths_file)
    sex_depths = load_sex_depths(sex_depths_file)

    log_info('load backbone coordinates')
    coordinates = load_coordinates(coordinates_file)
    families = group_families(coordinates)

    log_info('locate samples in input vcfs')
    sample2vcf = locate_samples(vcf_files, samples)

    tmp_dir = output_dir / 'tmp'
    ploidy_file = output_dir / 'grch38.ploidy'
    ploidy_file.write_text(GRCH38_PLOIDY)

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

    lookup = build_lookup(prepared)

    matched_lines, fill_lines = call_families(
        families=families,
        lookup=lookup,
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

    write_tsv(out_vcf, output_dir / f'{sample}.truth.tsv.gz', tmp_base)

    return sample


def write_tsv(vcf_file, output_file, tmp_dir):
    """bcftools query -> gzip TSV with columns: chrom pos id ref alt tgt."""
    tsv = tmp_dir / 'truth.tsv'
    with tsv.open('wt') as fh:
        fh.write('chrom\tpos\tid\tref\talt\ttgt\n')
    cmd = (''
           f"bcftools query"
           f"      -f '%CHROM\\t%POS\\t%ID\\t%REF\\t%ALT\\t[%TGT]\\n'"
           f"      {vcf_file}"
           f"      >> {tsv}"
           '')
    execute(cmd)
    with tsv.open('rb') as src, gzip.open(output_file, 'wb') as dst:
        shutil.copyfileobj(src, dst)
    return output_file


def prepare_sample(vcf_file, sample, genome_file, tmp_base):
    """Extract one sample, trim unused ALTs, left-align + parsimonious normalize
    (no multiallelic split) — mirrors the existing subset_samples chain."""
    vcf = (
        Vcf(vcf_file, tmp_base)
        .subset_samples({sample})
        .keep_format(fields=['GT'])
        .drop_info()
        .trim_alts()
        .normalize(genome_file)
        .uppercase()
        .exclude('ALT="."')
        .index()
    )
    return vcf.filepath


def build_lookup(prepared_vcf):
    """(chrom, pos, ref, ALT-set) -> genotype string, from the single-sample VCF."""
    lookup = dict()
    with gzip.open(prepared_vcf, 'rt') as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            items = line.rstrip('\n').split('\t')
            chrom, pos, _id, ref, alt = items[0], items[1], items[2], items[3], items[4]
            gt = items[9].split(':')[0]
            lookup[match_key(chrom, pos, ref, alt)] = gt
    return lookup


def call_families(families, lookup, gender, autosomes_depths, sex_depths, min_depth):
    """Per family: matched member(s) with the sample GT (siblings dropped), else one
    fill record per member (placeholder diploid GT, ploidy-corrected downstream)."""
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
            chrom = members[0]['chrom']
            pos = members[0]['pos']
            depth = depth_for(chrom, pos, gender, autosomes_depths, sex_depths)
            gt = '0/0' if is_homref(depth, min_depth) else './.'
            for member in members:
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


def load_autosomes_depths(depths_file):
    """(chrom, pos) -> depth_mean, from autosomes-depth.tsv (chr1..chr22 + chrM)."""
    depths = dict()
    opener = gzip.open if is_gzip(depths_file) else open
    with opener(depths_file, 'rt') as fh:
        idx = {c: i for i, c in enumerate(next(fh).rstrip('\n').split('\t'))}
        for line in fh:
            items = line.rstrip('\n').split('\t')
            key = (items[idx['chrom']], int(items[idx['pos']]))
            depths[key] = _to_float(items[idx['depth_mean']])
    return depths


def load_sex_depths(depths_file):
    """(chrom, pos) -> {'male': mean_male, 'female': mean_female}, from sex-depth.tsv."""
    depths = dict()
    opener = gzip.open if is_gzip(depths_file) else open
    with opener(depths_file, 'rt') as fh:
        idx = {c: i for i, c in enumerate(next(fh).rstrip('\n').split('\t'))}
        for line in fh:
            items = line.rstrip('\n').split('\t')
            key = (items[idx['chrom']], int(items[idx['pos']]))
            depths[key] = {
                'male': _to_float(items[idx['mean_male']]),
                'female': _to_float(items[idx['mean_female']]),
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
