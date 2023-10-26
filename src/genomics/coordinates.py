from pathlib import Path
import polars as pl
from .variant import Variant
from .vcf import Vcf
import gzip
from icecream import ic


def export_coordinates(
    input_vcf_file: Path,
    output_vcf_file: Path,
):

    coordinates = pl.read_csv(
        input_vcf_file,
        comment_char='#',
        has_header=False,
        new_columns=[
            'chrom', 'pos', 'id', 'ref', 'alt', 'qual', 'filter', 'info'
        ],
        separator='\t',
    )

    bag = list()
    for _, data in coordinates.groupby(['chrom', 'pos']):

        result = merge(data.to_dicts())

        if result:
            bag.append(str(result))

    tmp_file = output_vcf_file.with_suffix('.tmp.vcf')

    with gzip.open(input_vcf_file, 'rt') as ifh, tmp_file.open('wt') as ofh:
        for line in ifh:
            if line.startswith('#'):
                ofh.write(line)
                continue
            break
        for record in bag:
            ofh.write(f'{record}\n')

    Vcf(
        tmp_file,
        tmp_file.parents[0],
    ).bgzip().sort().index().move_to(output_vcf_file)

    return output_vcf_file


def merge(records: list):

    result = None
    for record in records:

        variant = Variant(
            chrom=record['chrom'],
            pos=record['pos'],
            id_=record['id'],
            ref=record['ref'],
            alt=record['alt'],
        )

        if not result:
            result = variant
        else:
            result = result.sync_alleles(variant, site_only=True)
    return result
