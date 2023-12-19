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


def export_gvcf_depths(
    gvcf_files: list[Path],
    coordinates_vcf_file: Path,
    output_dir: Path,
    n_threads: int,
    prod: bool = True,
):

    def jobs(gvcf_files, coordinates_file, depth_dir):
        for gvcf_file in gvcf_files:
            yield {
                'gvcf_file': gvcf_file,
                'output_dir': depth_dir / gvcf_file.name.split('.')[0],
                'coordinates_file': coordinates_file,
            }

    if prod:
        shutil.rmtree(output_dir, ignore_errors=True)

    output_dir.mkdir(exist_ok=True)

    coordinates_file = export_coordinates(coordinates_vcf_file, output_dir)

    depth_dir = output_dir / 'depth'
    depth_dir.mkdir(exist_ok=True)

    with ProcessPool(n_threads) as pool:
        for future in pool.uimap(
                parse_gvcf,
                jobs(gvcf_files, coordinates_file, depth_dir),
        ):
            sample_name = future[0]
            depth_df = pl.from_pandas(future[1])

            output_file = depth_dir / f'{sample_name}.tsv'
            depth_df.write_csv(output_file, has_header=True, separator='\t')

    output_filename = coordinates_vcf_file.name.split('.')[0] + '-depths.tsv'
    depth = summarize_depths(depth_dir)
    depth.write_csv(
        output_dir / f'{output_filename}',
        has_header=True,
        separator='\t',
    )


def export_cram_depths(
    cram_files: list[Path],
    genome_file: Path,
    coordinates_vcf_file: Path,
    output_dir: Path,
    n_threads: int,
    prod: bool = True,
):

    def jobs(cram_files, genome_file, coordinates_file):
        for cram_file in cram_files:
            yield {
                'cram_file': cram_file,
                'coordinates_file': coordinates_file,
                'genome_file': genome_file,
            }

    if prod:
        shutil.rmtree(output_dir, ignore_errors=True)

    output_dir.mkdir(exist_ok=True)

    coordinates_file = export_coordinates(coordinates_vcf_file, output_dir)

    depths_dir = output_dir / 'depths'
    depths_dir.mkdir(exist_ok=True)

    with ProcessPool(n_threads) as pool:
        for future in pool.uimap(
                parse_cram,
                jobs(cram_files, genome_file, coordinates_file),
        ):
            sample_name = future[0]
            depth_df = future[1]

            output_file = depths_dir / f'{sample_name}.tsv'
            depth_df.write_csv(output_file, has_header=True, separator='\t')
            # depth_df.to_csv(output_file, header=True, index=False, sep='\t')

    depths = summarize_depths(depths_dir)

    output_filename = f'{coordinates_vcf_file.name.split(".")[0]}-depth.tsv.gz'

    coordinates = Vcf(coordinates_vcf_file, output_dir).to_df(site_only=True)
    depths = coordinates.join(depths, on=['chrom', 'pos'])
    df2tsv(depths, output_dir / f'{output_filename}')


def summarize_depths(depth_dir):
    bag = []
    for f in depth_dir.glob('*.tsv'):
        d = pl.read_csv(
            f,
            has_header=True,
            separator='\t',
            infer_schema_length=0,
        ).with_columns(pl.col('depth').cast(pl.Float64))
        bag.append(d)

    data = pl.concat(bag)

    result = dict()

    bag = list()
    for key, df in data.groupby(['chrom', 'pos']):
        bag.append({
            'chrom': key[0],
            'pos': key[1],
            'depth': np.median(df['depth'])
        })

    return pl.from_dicts(bag)


def parse_gvcf(args):
    gvcf_file = args['gvcf_file']
    print(gvcf_file)
    coordinates_file = args['coordinates_file']
    output_dir = args['output_dir']
    output_dir.mkdir(exist_ok=True)

    cmd = (
        ''
        f'bcftools view '
        f'    -R {coordinates_file}'
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

    chrom2coordinate = load_chrom2coordinate(coordinates_file)
    chrom2depth = load_chrom2depth(output_dir / 'raw_depths.tsv')

    depths = merge(chrom2coordinate, chrom2depth)

    coordinates = pd.read_csv(
        coordinates_file,
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
                'pos': items[1],
                'ref': items[2],
                'depth': items[3],
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

    coordinates = pl.read_csv(
        coordinates_file,
        has_header=False,
        new_columns=['chrom', 'pos'],
        separator='\t',
        infer_schema_length=0,
    )
    depths = pl.from_dicts(bag)
    depths = coordinates.join(depths, on=['chrom', 'pos'], how='left')

    sample_name = cram_file.name.split('.')[0]
    return sample_name, depths


def export_coordinates(input_file, output_dir):
    output_file = output_dir / COORDINATES_FILENAME
    coordinates = Vcf(
        input_file,
        output_dir,
    ).to_df(site_only=True).select([
        'chrom',
        'pos',
    ]).unique()

    coordinates.write_csv(
        output_file,
        has_header=False,
        separator='\t',
    )

    return output_file
