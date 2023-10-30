from pathlib import Path
import polars as pl
from .variant import Variant
from .vcf import Vcf
import gzip
from icecream import ic
import shutil


def export_coordinates(
    input_vcf_file: Path,
    genome_file: Path,
    genome_index_file: Path,
    output_dir: Path,
    prod: bool = True,
):

    if prod:
        shutil.rmtree(output_dir, ignore_errors=True)

    output_dir.mkdir(exist_ok=True)

    coordinates = preprocess(
        input_vcf_file,
        output_dir,
        genome_file,
        genome_index_file,
    )

    bag = list()
    for _, data in coordinates.groupby(['chrom', 'pos']):
        record = merge(data.to_dicts())
        bag.append(str(record))

    tmp_coordinates_file = output_dir / 'coordinates.vcf'

    copy_header(input_vcf_file, tmp_coordinates_file)

    with tmp_coordinates_file.open('at') as fh:
        for record in bag:
            fh.write(f'{record}\n')

    output_vcf_file = output_dir / (input_vcf_file.name \
            .replace('.vcf', '') \
            .replace('.bgz', '') \
            .replace('.gz', '',) + '-coordinates.vcf.bgz')

    postprocess(
        tmp_coordinates_file,
        genome_file,
        output_dir,
        output_vcf_file,
    )

    print(output_vcf_file)


def copy_header(input_file, output_file):
    with gzip.open(input_file, 'rt') as ifh, output_file.open('wt') as ofh:
        for line in ifh:
            if line.startswith('#'):
                ofh.write(line)
                continue
            break


def merge(records: list):

    tmp = None
    for record in records:

        variant = Variant(
            chrom=record['chrom'],
            pos=record['pos'],
            id_=record['id'],
            ref=record['ref'],
            alt=record['alt'],
        )

        if not tmp:
            tmp = variant
        else:
            tmp = tmp.sync_alleles(variant, site_only=True)
    return tmp


def preprocess(input_vcf_file, output_dir, genome_file, genome_index_file):

    preprocessed_vcf_file = Vcf(input_vcf_file, output_dir)\
            .bgzip()\
            .drop_qual() \
            .drop_filter()\
            .drop_info()\
            .drop_gt()\
            .fix_header(genome_index_file) \
            .normalize(genome_file) \
            .uppercase() \
            .filepath
    return pl.read_csv(
        preprocessed_vcf_file,
        comment_char='#',
        has_header=False,
        new_columns=[
            'chrom', 'pos', 'id', 'ref', 'alt', 'qual', 'filter', 'info'
        ],
        separator='\t',
    )


def postprocess(input_vcf_file, genome_file, output_dir, output_vcf_file):

    Vcf(input_vcf_file, output_dir) \
            .bgzip() \
            .sort() \
            .normalize(genome_file) \
            .uppercase() \
            .index() \
            .move_to(output_vcf_file)
