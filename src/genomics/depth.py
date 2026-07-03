import shutil
import numpy as np
import pandas as pd
from pathlib import Path
from subprocess import PIPE, Popen
from collections import OrderedDict
from pathos.multiprocessing import ProcessPool
import multiprocess.context as ctx

from .utils import load_dict, init_logging, log_start, log_info


ctx._force_start_method('spawn')

AUTOSOMES = {f'chr{i}' for i in range(1, 23)}
MITO = 'chrM'
SEX = {'chrX', 'chrY'}
CONTIG_ORDER = [f'chr{i}' for i in range(1, 23)] + [MITO, 'chrX', 'chrY']


def export_cram_depths(
    crams_file: Path,
    genders_file: Path,
    genome_file: Path,
    output_dir: Path,
    n_threads: int,
    prod: bool = False,
):

    if prod:
        shutil.rmtree(output_dir, ignore_errors=True)

    output_dir.mkdir(parents=True, exist_ok=True)

    init_logging(output_dir / 'depth.log')

    info = OrderedDict()
    info['crams-file'] = crams_file
    info['genders-file'] = genders_file
    info['genome-file'] = genome_file
    info['output-dir'] = output_dir
    info['n-threads'] = n_threads

    log_start('CRAM-based Whole-Contig Read Depth', info)

    sample2cram = load_dict(crams_file)
    sample2gender = load_dict(genders_file)

    validate_genders(sample2cram, sample2gender)

    contigs = target_contigs(next(iter(sample2cram.values())))
    n_samples = len(sample2cram)

    autosomes_file = output_dir / 'autosomes-depth.tsv'
    sex_file = output_dir / 'sex-depth.tsv'

    with autosomes_file.open('wt') as a_fh, sex_file.open('wt') as s_fh:
        a_fh.write('chrom\tpos\tdepth_mean\tn_samples\n')
        s_fh.write('chrom\tpos\tmean_male\tn_male\tmean_female\tn_female\n')

        log_info('Compute depths')

        # One job per (sample, contig): all n_threads workers stay busy across
        # contigs instead of serializing one contig at a time. Each samtools call
        # is still scoped to a single contig (-r), so a contig's per-position
        # accumulation stays bounded to that contig. Results arrive unordered;
        # a contig is finalized once all its samples are in, and completed
        # contigs are flushed in CONTIG_ORDER for deterministic output.
        pending = {contig: [] for contig in contigs}
        done = {}
        write_idx = 0

        with ProcessPool(n_threads) as pool:
            for sample, gender, contig, depths in pool.uimap(
                    contig_depth,
                    cram_depth_jobs(sample2cram, sample2gender, genome_file, contigs)):

                bucket = pending[contig]
                bucket.append((gender, depths))

                if len(bucket) < n_samples:
                    continue

                if classify_contig(contig) == 'single':
                    mean, n = mean_single([vec for _, vec in bucket])
                    done[contig] = ('single', mean, n)
                else:
                    done[contig] = ('sex', *mean_per_gender(bucket))
                del pending[contig]

                while write_idx < len(contigs) and contigs[write_idx] in done:
                    c = contigs[write_idx]
                    payload = done.pop(c)
                    if payload[0] == 'single':
                        write_single(a_fh, c, payload[1], payload[2])
                    else:
                        write_sex(s_fh, c, payload[1], payload[2], payload[3], payload[4])
                    write_idx += 1

    log_info('Done')


def cram_depth_jobs(sample2cram, sample2gender, genome_file, contigs):
    for contig in contigs:
        for sample, cram_file in sample2cram.items():
            yield {
                'cram_file': cram_file,
                'genome_file': genome_file,
                'contig': contig,
                'sample': sample,
                'gender': sample2gender[sample],
            }


def classify_contig(chrom):
    if chrom in AUTOSOMES or chrom == MITO:
        return 'single'
    if chrom in SEX:
        return 'sex'
    return 'exclude'


def validate_genders(sample2cram, sample2gender):
    for sample in sample2cram:
        gender = sample2gender.get(sample)
        if gender not in ('male', 'female'):
            raise ValueError(
                f'sample {sample!r} has missing or invalid gender: {gender!r}')


def target_contigs(cram_file):
    cmd = f'samtools view -H {cram_file}'

    present = set()
    with Popen(cmd, shell=True, text=True, stdout=PIPE) as proc:
        for line in proc.stdout:
            if not line.startswith('@SQ'):
                continue
            for field in line.rstrip('\n').split('\t'):
                if field.startswith('SN:'):
                    present.add(field[3:])
        proc.wait()
        if proc.returncode:
            raise Exception(cmd)

    return [contig for contig in CONTIG_ORDER if contig in present]


def contig_depth(args):
    cram_file = args['cram_file']
    genome_file = args['genome_file']
    contig = args['contig']
    sample = args['sample']
    gender = args['gender']

    cmd = (''
           f'samtools depth -aa'
           f'    -r {contig}'
           f'    --reference {genome_file}'
           f'    {cram_file}'
           '')

    depths = []
    with Popen(cmd, shell=True, text=True, stdout=PIPE) as proc:
        for line in proc.stdout:
            depths.append(int(line.rstrip('\n').split('\t')[2]))
        proc.wait()
        if proc.returncode:
            raise Exception(cmd)

    return sample, gender, contig, np.array(depths, dtype=np.int64)


def mean_single(depth_vectors):
    total = None
    for vec in depth_vectors:
        total = vec.astype(np.float64) if total is None else total + vec
    return total / len(depth_vectors), len(depth_vectors)


def mean_per_gender(items):
    males = [vec for gender, vec in items if gender == 'male']
    females = [vec for gender, vec in items if gender == 'female']

    mean_male = mean_single(males)[0] if males else None
    mean_female = mean_single(females)[0] if females else None

    return mean_male, len(males), mean_female, len(females)


def write_single(fh, contig, mean, n_samples):
    for i in range(len(mean)):
        fh.write(f'{contig}\t{i + 1}\t{mean[i]}\t{n_samples}\n')


def write_sex(fh, contig, mean_male, n_male, mean_female, n_female):
    length = len(mean_male) if mean_male is not None else len(mean_female)
    for i in range(length):
        vmale = mean_male[i] if mean_male is not None else ''
        vfemale = mean_female[i] if mean_female is not None else ''
        fh.write(f'{contig}\t{i + 1}\t{vmale}\t{n_male}\t{vfemale}\t{n_female}\n')


# def export_gvcf_depths(
#     gvcf_files: list[Path],
#     coordinates_file: Path,
#     output_dir: Path,
#     n_threads: int,
#     prod: bool = True,
# ):
#
#     def jobs(gvcf_files, coordinates_file, depth_dir):
#         for gvcf_file in gvcf_files:
#             yield {
#                 'gvcf_file': gvcf_file,
#                 'output_dir': depth_dir / gvcf_file.name.split('.')[0],
#                 'coordinates_file': coordinates_file,
#             }
#
#     if prod:
#         shutil.rmtree(output_dir, ignore_errors=True)
#
#     output_dir.mkdir(exist_ok=True)
#
#     coordinates_file = export_coordinates(coordinates_file, output_dir)
#
#     depth_dir = output_dir / 'depth'
#     depth_dir.mkdir(exist_ok=True)
#
#     with ProcessPool(n_threads) as pool:
#         for future in pool.uimap(
#                 parse_gvcf,
#                 jobs(gvcf_files, coordinates_file, depth_dir),
#         ):
#             sample_name = future[0]
#             depth_df = pl.from_pandas(future[1])
#
#             output_file = depth_dir / f'{sample_name}.tsv'
#             depth_df.write_csv(output_file, has_header=True, separator='\t')
#
#     output_filename = coordinates_file.name.split('.')[0] + '-depths.tsv'
#     depth = summarize_depths(depth_dir)
#     depth.write_csv(
#         output_dir / f'{output_filename}',
#         has_header=True,
#         separator='\t',
#     )

def parse_gvcf(args):
    gvcf_file = args['gvcf_file']
    print(gvcf_file)
    unique_coordinates_file = args['unique_coordinates_file']
    output_dir = args['output_dir']
    output_dir.mkdir(exist_ok=True)

    cmd = (
        ''
        f'bcftools view '
        f'    -R {unique_coordinates_file}'
        f'    {gvcf_file}'
        f"    | bcftools query -f '[%CHROM\t%POS\t%END\t%REF\t%DP\t%MIN_DP\n]'"
        '')

    with Popen(cmd, shell=True, text=True, stdout=PIPE) as proc:
        bag = list()
        for line in proc.stdout:
            record = parse_gvcf_depth(line)
            bag.append(record)

        proc.wait()

        if proc.returncode:
            raise Exception(cmd)

    pd.DataFrame.from_records(bag).to_csv(
        output_dir / 'raw_depths.tsv',
        header=True,
        index=False,
        sep='\t',
    )

    chrom2coordinate = load_chrom2coordinate(unique_coordinates_file)
    chrom2depth = load_chrom2depth(output_dir / 'raw_depths.tsv')

    depths = merge(chrom2coordinate, chrom2depth)

    coordinates = pd.read_csv(
        unique_coordinates_file,
        header=None,
        names=['chrom', 'pos'],
        sep='\t',
        dtype='str',
    )

    print('coordinates', len(coordinates))
    depths = coordinates.merge(depths, on=['chrom', 'pos'], how='left')
    print('depths', len(depths))

    sample_name = gvcf_file.name.split('.')[0]
    return sample_name, depths


def load_chrom2coordinate(coordinates_file):
    bag = dict()
    with coordinates_file.open('rt') as fh:
        for line in fh:
            items = line.strip().split('\t')
            chrom = items[0]
            pos = items[1]

            if chrom not in bag:
                bag[chrom] = list()

            bag[chrom].append(pos)

    bag2 = dict()
    for k, v in bag.items():
        bag2[k] = sorted(v)
    return bag2


def load_chrom2depth(raw_depth_file):
    bag = dict()
    with raw_depth_file.open('rt') as fh:
        columns = {
            column: idx
            for idx, column in enumerate(next(fh).strip().split('\t'))
        }
        for line in fh:
            items = line.strip().split('\t')
            chrom = items[columns['chrom']]
            start = items[columns['start']]
            end = items[columns['end']]
            depth = items[columns['depth']]

            if chrom not in bag:
                bag[chrom] = list()

            bag[chrom].append({
                'chrom': chrom,
                'start': start,
                'end': end,
                'depth': depth,
            })
    bag2 = dict()

    for k, v in bag.items():
        bag2[k] = sorted(v, key=lambda x: x['start'])

    return bag2


def merge(chrom2coordinate, chrom2depth):

    bag = list()
    for chrom in chrom2coordinate.keys():

        if chrom not in chrom2depth:
            print(chrom)
            continue
        cs = chrom2coordinate[chrom]
        ds = chrom2depth[chrom]

        i, j = 0, 0

        while i < len(cs) and j < len(ds):
            pos = cs[i]
            d = ds[j]
            if pos >= d['start'] and pos <= d['end']:
                bag.append({
                    'chrom': chrom,
                    'pos': pos,
                    'depth': d['depth'],
                })
            if pos <= d['end']:
                i += 1
            else:
                j += 1

    return pd.DataFrame.from_records(bag)


def parse_gvcf_depth(line):
    items = line.strip().split('\t')
    chrom = items[0]
    start = items[1]
    end = items[2]
    ref = items[3]
    dp = items[4]
    min_dp = items[5]

    if 'm' in chrom.lower():
        depth = 30
    if dp.isnumeric():
        depth = dp
    elif min_dp.isnumeric():
        depth = min_dp
    else:
        depth = np.nan

    return {
        'chrom': chrom,
        'start': start,
        'end': end,
        'ref': ref,
        'depth': depth,
    }
