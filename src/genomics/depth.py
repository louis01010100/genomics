import gzip
import shutil
from pathlib import Path
from subprocess import PIPE, STDOUT, Popen
import logging
from icecream import ic

import numpy as np
import polars as pl
from pathos.multiprocessing import ProcessPool
from .vcf import Vcf

COORDINATES_FILENAME = 'coordinates.tsv'


def export_gvcf_depth(
    gvcf_files: list[Path],
    coordinates_vcf_file: Path,
    output_dir: Path,
    n_threads: int,
    prod: bool = False,
):

    def jobs(gvcf_files, coordinates_file):
        for gvcf_file in gvcf_files:
            yield {
                'gvcf_file': gvcf_file,
                'coordinates_file': coordinates_file,
            }

    if not debug:
        shutil.rmtree(output_dir, ignore_errors=True)

    output_dir.mkdir(exist_ok=True)

    coordinates_file = export_coordinates(coordinates_vcf_file, output_dir)

    depth_dir = output_dir / 'depth'
    depth_dir.mkdir(exist_ok=True)

    with ProcessPool(n_threads) as pool:
        for future in pool.uimap(
                parse_gvcf,
                jobs(gvcf_files, coordinates_file),
        ):
            sample_name = future[0]
            depth_df = future[1]

            output_file = depth_dir / f'{sample_name}.tsv'
            depth_df.write_csv(output_file, has_header=True, separator='\t')

    depth = summarize_depth(depth_dir)
    depth.write_csv(
        output_dir / 'depth.tsv',
        has_header=True,
        separator='\t',
    )


def export_cram_depth(
    cram_files: list[Path],
    genome_file: Path,
    coordinates_vcf_file: Path,
    output_dir: Path,
    n_threads: int,
    prod: bool = False,
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

    depth_dir = output_dir / 'depth'
    depth_dir.mkdir(exist_ok=True)

    with ProcessPool(n_threads) as pool:
        for future in pool.uimap(
                parse_cram,
                jobs(cram_files, genome_file, coordinates_file),
        ):
            sample_name = future[0]
            depth_df = future[1]

            output_file = depth_dir / f'{sample_name}.tsv'
            depth_df.write_csv(output_file, has_header=True, separator='\t')

    depth = summarize_depth(depth_dir)
    depth.write_csv(
        output_dir / 'depth.tsv',
        has_header=True,
        separator='\t',
    )


def summarize_depth(depth_dir):
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
    coordinates_file = args['coordinates_file']

    cmd = (''
           f'bcftools view '
           f'    -R {coordinates_file}'
           f'    {gvcf_file}'
           f"    | bcftools query -f '[%CHROM\t%POS\t%REF\t%DP\t%MIN_DP\n]'"
           '')

    with Popen(cmd, shell=True, text=True, stdout=PIPE) as proc:
        bag = []
        for line in proc.stdout:
            items = line.strip().split('\t')
            chrom = items[0]
            pos = items[1]
            ref = items[2]
            dp = items[3]
            min_dp = items[4]

            if dp.isnumeric():
                depth = dp
            elif min_dp.isnumeric():
                depth = min_dp
            else:
                depth = np.nan

            bag.append({
                'chrom': chrom,
                'pos': pos,
                'ref': ref,
                'depth': depth
            })

        proc.wait()

        if proc.returncode:
            raise Exception(cmd)
    sample_name = gvcf_file.name.split('.')[0]
    return sample_name, pl.from_dicts(bag)


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
    sample_name = cram_file.name.split('.')[0]
    return sample_name, pl.from_dicts(bag)


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
