import gzip
import shutil
import pandas as pd
from pathlib import Path
from subprocess import PIPE, STDOUT, Popen
import logging
from icecream import ic
from .utils import df2tsv

import numpy as np
import polars as pl
from pathos.multiprocessing import ProcessPool
from .vcf import Vcf
import multiprocess.context as ctx

ctx._force_start_method('spawn')

COORDINATES_FILENAME = 'coordinates.tsv'

def export_cram_depths(
    crams_file: Path,
    genome_file: Path,
    coordinates_file: Path,
    output_dir: Path,
    n_threads: int,
    prod: bool = False,
):

    def jobs(crams, genome_file, atomized_coordinates_file):
        for cram in crams:
            yield {
                'cram_file': cram['cram_path'],
                'sample': cram['sample'],
                'gender': cram['gender'],
                'atomized_coordinates_file': atomized_coordinates_file,
                'genome_file': genome_file,
            }

    if prod:
        shutil.rmtree(output_dir, ignore_errors=True)

    output_dir.mkdir(exist_ok=True)

    crams = pl.read_csv(crams_file, has_header = True, separator = '\t').to_dicts()

    coordinates = atomize_coordinates(coordinates_file)
    coordinates.write_csv(output_dir / 'atomized_coordinates_full.tsv', include_header = True, separator = '\t')
    coordinates\
            .select(['chrom', 'pos'])\
            .write_csv(output_dir / 'atomized_coordinates.tsv', include_header = False, separator = '\t')

    depths_dir = output_dir / 'depths'
    depths_dir.mkdir(exist_ok=True)

    # with ProcessPool(n_threads) as pool:
    #     for future in pool.uimap(
    #             parse_cram,
    #             jobs(crams, genome_file, output_dir / 'atomized_coordinates.tsv'),
    #     ):
    #         sample_name = future[0]
    #         gender = future[1]
    #         depth_df = pl.from_pandas(future[2])
    #         # depth_df = future[1]
    #
    #
    #         depth_df = coordinates\
    #                 .join(depth_df, on = ['chrom', 'pos'], how = 'left')\
    #                 .with_columns(
    #                         pl.col('depth').fill_null(0), 
    #                         pl.lit(sample_name).alias('sample'), 
    #                         pl.lit(gender).alias('gender')
    #                 )
    #
    #
    #         output_file = depths_dir / f'{sample_name}.tsv'
    #         depth_df.write_csv(output_file, include_header=True, separator='\t')

    depths = summarize_depths(depths_dir)

    output_filename = f'{coordinates_file.name.split(".")[0]}-depth.tsv'

    result = pl.read_csv(coordinates_file, has_header = True, separator = '\t')
    print(result)
    print(depths)
    result = result.join(depths, on=['chrom', 'start', 'end'], how = 'left')

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

def summarize_depths(depth_dir):

    bag = dict()
    for f in depth_dir.glob('*.tsv'):

        with f.open('rt') as fh:
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
                gender = record[col2idx['gender']]

                key = (chrom, start, end)

                if key not in bag:
                    bag[key] = list()

                if gender == 'male' and chrom == 'chrX':
                    continue
                if gender == 'female' and chrom == 'chrY':
                    continue

                bag[key].append(depth)

    result = list()

    for key, depths in bag.items():
        print(depths)
        result.append({
            'chrom': key[0],
            'start': key[1],
            'end': key[2],
            'depth_median': int(np.median(depths)),
            'depth_max': np.max(depths),
            'depth_min': np.min(depths),
            'n_zero_depth': len([x for x in depths if x == 0]),
            'n_depths': len(depths),
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
    atomized_coordinates_file = args['atomized_coordinates_file']

    cmd = (''
           f'samtools mpileup'
           f'    -l {atomized_coordinates_file}'
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

    # coordinates = pd.read_csv(
    #     coordinates_file,
    #     header=None,
    #     names=['chrom', 'pos'],
    #     sep='\t',
    #     dtype='str',
    # )
    # depths = pd.DataFrame.from_records(bag)
    # depths = coordinates.merge(depths, on=['chrom', 'pos'], how='left')

    # coordinates = pl.read_csv(
    #     unique_coordinates_file,
    #     has_header=False,
    #     new_columns=['chrom', 'start', 'end'],
    #     separator='\t',
    #     infer_schema_length=0,
    # )

    # depths = pl.from_dicts(bag)

    depths = pd.DataFrame.from_records(bag)

    return sample, gender, depths

def atomize_coordinates(input_file):
    bag = list()
    with open(input_file, 'rt') as fh:
        header = next(fh)
        col2idx = {col: idx for idx, col in enumerate(header.strip().split('\t'))}

        for line in fh:
            record = line.strip().split('\t')

            chrom = record[col2idx['chrom']]
            start = str2int(record[col2idx['start']])
            end = str2int(record[col2idx['end']])

            for i in range(0, end - start + 1):
                bag.append({'chrom': chrom, 'pos': start + i, 'start': start, 'end': end})

    return pl.from_dicts(bag)



def export_unique_coordinates(input_file, output_dir):
    output_file = output_dir / COORDINATES_FILENAME
    data = pl.read_csv(input_file, has_header = True, separator = '\t')
    data = data.with_columns(pl.col('start') -1)
    data.select([ 'chrom', 'start' , 'end']) \
        .unique()\
        .sort(['chrom', 'start', 'end']) \
        .write_csv(
            output_file,
            include_header=False,
            separator='\t',
    )

    return output_file
