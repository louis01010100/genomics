import gzip
import shutil
import pandas as pd
from pathlib import Path
from subprocess import PIPE, STDOUT, Popen
from collections import OrderedDict
import logging
from icecream import ic
from .utils import df2tsv, is_gzip, load_dict, init_logging, log_start, log_stop, log_info
from .genome import Genome
from .variant import Variant

import numpy as np
import polars as pl
from pathos.multiprocessing import ProcessPool
from .vcf import Vcf
import multiprocess.context as ctx


ctx._force_start_method('spawn')

COORDINATES_FILENAME = 'coordinates.tsv'

def export_cram_depths(
    coordinates_file: Path,
    crams_file: Path,
    genders_file: Path,
    genome_file: Path,
    output_dir: Path,
    n_threads: int,
    prod: bool = False,
):


    if prod:
        shutil.rmtree(output_dir, ignore_errors=True)

    output_dir.mkdir(parents = True, exist_ok=True)

    log_file = output_dir / 'depth.log'

    init_logging(log_file)

    info = OrderedDict()

    info['coordinates-file'] = coordinates_file
    info['crams-file'] = crams_file
    info['genders-file'] = genders_file
    info['genome-file'] = genome_file
    info['ouput_dir'] = output_dir
    info['n-threads'] = n_threads

    log_start('CRAM-based Read Depth Export', info)

    sample2cram = load_dict(crams_file)
    sample2gender = load_dict(genders_file)

    log_info('Expand Coordinates')
    coordinates_expanded = expand_coordinates(coordinates_file, genome_file)
    coordinates_expanded = pl.read_csv(
            output_dir / 'coordinates_expanded_details.tsv', has_header = True, separator = '\t')
    coordinates_expanded_file = output_dir / 'coordinates_expanded.tsv'
    coordinates_expanded.select(['chrom', 'pos']).unique().sort(['chrom', 'pos']).write_csv(
            coordinates_expanded_file, include_header = False, separator = '\t')


    depths_dir = output_dir / 'depths'
    depths_dir.mkdir(exist_ok=True)

    log_info('Parse CRAM')
    with ProcessPool(n_threads) as pool:
        for future in pool.uimap(
                parse_cram,
                parse_cram_depths_jobs(sample2cram, sample2gender, genome_file, coordinates_expanded_file)
        ):
            sample_name = future[0]
            gender = future[1]
            depth_df = future[2]

            depth_df = coordinates_expanded.join(depth_df, on = ['chrom', 'pos'], how = 'left')\
                    .with_columns(
                            pl.col('depth').fill_null(0), 
                    )

            output_file = depths_dir / f'{sample_name}.tsv'
            depth_df.write_csv(output_file, include_header=True, separator='\t')

    log_info('Summarize Depths')
    depths = summarize_depths(depths_dir, sample2gender)

    output_filename = f'{coordinates_file.name.split(".")[0]}-depth.tsv'

    result = coordinates_expanded.join(depths, on=['chrom', 'pos'], how = 'left')

    df2tsv(result, output_dir / f'{output_filename}')


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

# male X and female Y are excluded from the summary
def summarize_depths(depth_dir, sample2gender):

    bag = dict()
    for f in depth_dir.glob('*.tsv'):

        with f.open('rt') as fh:
            sample = f.name[0:-4]
            gender = sample2gender[sample]
            col2idx = {col: idx for idx, col in enumerate(next(fh).strip().split('\t'))}

            n_fields = len(col2idx)

            for line in fh:
                record = line.strip().split('\t')

                record = pad(record, n_fields)
                chrom = record[col2idx['chrom']]
                pos = str2int(record[col2idx['pos']])
                ref = record[col2idx['ref']]
                start = str2int(record[col2idx['start']])
                end = str2int(record[col2idx['end']])
                depth = str2int(record[col2idx['depth']])

                key = (chrom, pos)

                if key not in bag:
                    bag[key] = list()

                if gender == 'male' and chrom == 'chrX':
                    continue
                if gender == 'female' and chrom == 'chrY':
                    continue

                bag[key].append(depth)

    result = list()

    for key, depths in bag.items():
        result.append({
            'chrom': key[0],
            'pos': key[1],
            'depth_median': int(np.median(depths)),
            'depth_max': np.max(depths),
            'depth_min': np.min(depths),
            'n_samples_zero_depth': len([x for x in depths if x == 0]),
            'n_samples': len(depths),
        })

    return pl.from_dicts(result)

def str2int(v):
    if v is None:
        return 0
    return int(v)

def pad(record: list, n_fields: int, padding = None):
    if len(record) == n_fields:
        return record

    record.extend([None] * (n_fields - len(record)))

    return record

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


def parse_cram(args):
    cram_file = args['cram_file']
    genome_file = args['genome_file']
    sample = args['sample']
    gender = args['gender']
    coordinates_file = args['coordinates_file']

    cmd = (''
           f'samtools mpileup'
           f'    -l {coordinates_file}'
           f'    -f {genome_file}'
           f'    {cram_file}'
           '')

    with Popen(cmd, shell=True, text=True, stdout=PIPE) as proc:
        bag = []
        for line in proc.stdout:
            items = line.strip().split('\t', 4)
            record = {
                'chrom': items[0],
                'pos': int(items[1]),
                'ref': items[2],
                'depth': int(items[3]),
            }

            bag.append(record)

        proc.wait()

        if proc.returncode:
            raise Exception(cmd)


    depths = pl.from_dicts(bag)

    return sample, gender, depths




def expand_coordinates(coordinates_vcf_file, genome_file):

    def create_coordinates(items, genome):
        chrom = items[0]
        pos = int(items[1])
        id_ = items[2]
        ref = items[3]
        alt = items[4]

        v = Variant(chrom = chrom, pos = pos, id_ = id_, ref = ref, alt = alt)
        region = v.expand(genome.chromosome(chrom)).region

        bag = list()

        for _pos in range(region.start, region.end + 1):

            bag.append({
                'chrom': chrom,
                'pos': _pos,
                'id': id_,
                'pos_original': pos,
                'ref_original': ref,
                'alt_original': alt,
                'start': region.start,
                'end': region.end,
            })

        return bag

    genome = Genome(genome_file)

    bag = list()
    if is_gzip(coordinates_vcf_file):
        with gzip.open(coordinates_vcf_file, 'rt') as fh:
            for line in fh:
                if line.startswith('##'):
                    continue
                break

            for line in fh:
                items = line.strip('\n').split('\t')

                record = create_coordinates(items, genome)
                bag.extend(record)
    else:
        with coordinates_vcf_file.open('rt') as fh:
            for line in fh:
                if line.startswith('##'):
                    continue
                break

            for line in fh:
                items = line.strip('\t').split('\t')

                record = create_coordinates(items, genome)
                bag.extend(record)

    return  pl.from_dicts(bag).sort(['chrom', 'pos'])

def parse_cram_depths_jobs(sample2cram, sample2gender, genome_file, coordinates_expanded_file):
    for sample, cram_file in sample2cram.items():
        yield {
            'cram_file': cram_file,
            'sample': sample,
            'gender': sample2gender[sample],
            'coordinates_file': coordinates_expanded_file,
            'genome_file': genome_file,
        }
