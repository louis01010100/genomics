import shutil
from pathlib import Path

import polars as pl

from .vcf import Vcf, concat


def process(
    clinvar_vcf_file: Path,
    clinvar_papu_vcf_file: Path,
    genome_file: Path,
    genome_index_file: Path,
    output_dir: Path,
):

    shutil.rmtree(output_dir, ignore_errors=True)
    output_dir.mkdir(parents=True)

    tmp_dir = output_dir / 'tmp'
    tmp_dir.mkdir(parents=True)

    chrom_map = create_chrom_map()

    chrom_map_file = tmp_dir / 'chrom_map.tsv'

    chrom_map.write_csv(chrom_map_file, has_header=False, separator='\t')



    clinvar_vcf = Vcf(
        clinvar_vcf_file, tmp_dir,) \
                .bgzip()\
                .index()\
                .rename_chroms(chrom_map_file)\
                .include_chroms(target_chroms()) \
                .fix_header(genome_index_file)\
                .index()\
                .normalize(genome_file)\
                .index()

    clinvar_papu_vcf = Vcf(
        clinvar_papu_vcf_file, tmp_dir,) \
                .bgzip()\
                .index()\
                .rename_chroms(chrom_map_file)\
                .include_chroms(target_chroms()) \
                .fix_header(genome_index_file)\
                .index()\
                .normalize(genome_file)\
                .index()

    concat(
        vcf_files=[clinvar_vcf.filepath, clinvar_papu_vcf.filepath],
        output_file=output_dir / 'clinvar.vcf.bgz',
        tmp_dir=tmp_dir,
    )

    Vcf(output_dir / 'clinvar.vcf.bgz', output_dir).to_df(
        '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%INFO/CLNSIG\t%INFO/CLNREVSTAT\t%INFO/CLNSIGCONF\t%INFO/CLNDN\n',
        null_values=['.'],
    ).rename({
        'info_clnsig': 'clnsig',
        'info_clnsigconf': 'clnsigconf',
        'info_clnrevstat': 'clnrevstat',
        'info_clndn': 'clndn',
    }).write_csv(
        output_dir / 'clinvar.tsv',
        has_header=True,
        separator='\t',
    )


def target_chroms():
    bag = [f'chr{i}' for i in range(1, 23)]

    bag.append('chrX')
    bag.append('chrY')
    bag.append('chrM')

    return set(bag)


def create_chrom_map():
    bag = []
    for i in range(1, 23):
        bag.append({'name1': i, 'name2': f'chr{i}'})

    bag.append({'name1': 'X', 'name2': 'chrX'})
    bag.append({'name1': 'Y', 'name2': 'chrY'})
    bag.append({'name1': 'MT', 'name2': 'chrM'})

    return pl.from_dicts(bag)
