#!/usr/bin/env python
import shutil
from pathlib import Path

from .vcf import Vcf

CHROM_MAP_TEMPLATE = {
    'NC_000001': 'chr1',
    'NC_000002': 'chr2',
    'NC_000003': 'chr3',
    'NC_000004': 'chr4',
    'NC_000005': 'chr5',
    'NC_000006': 'chr6',
    'NC_000007': 'chr7',
    'NC_000008': 'chr8',
    'NC_000009': 'chr9',
    'NC_000010': 'chr10',
    'NC_000011': 'chr11',
    'NC_000012': 'chr12',
    'NC_000013': 'chr13',
    'NC_000014': 'chr14',
    'NC_000015': 'chr15',
    'NC_000016': 'chr16',
    'NC_000017': 'chr17',
    'NC_000018': 'chr18',
    'NC_000019': 'chr19',
    'NC_000020': 'chr20',
    'NC_000021': 'chr21',
    'NC_000022': 'chr22',
    'NC_000023': 'chrX',
    'NC_000024': 'chrY',
    'NC_012920': 'chrM',
}


def process(
    dbsnp_vcf_file: Path,
    output_dir: Path,
    genome_file: Path,
    genome_index_file: Path,
    n_threads: int,
):
    if output_dir.exists():
        shutil.rmtree(output_dir)

    tmp_dir = output_dir / 'tmp'
    tmp_dir.mkdir(parents=True)

    vcf = Vcf(dbsnp_vcf_file, tmp_dir, n_threads)

    vcf = vcf.bgzip().index()
    contigs = vcf.contigs

    chrom_map_file = tmp_dir / 'chrommap.tsv'

    export_chrom_map(contigs, CHROM_MAP_TEMPLATE, chrom_map_file)

    vcf.rename_chroms(chrom_map_file)\
        .include_chroms(CHROM_MAP_TEMPLATE.values())\
        .fix_header(genome_index_file)\
        .drop_info() \
        .normalize(genome_file) \
        .move_to(output_dir / 'dbsnp.vcf.bgz')


def export_chrom_map(contigs: set, template: dict, output_file):

    bag = dict()
    for contig in contigs:
        k = contig.split('.')[0]

        if k in template:
            bag[contig] = template[k]

    with output_file.open('wt') as fh:

        for k, v in bag.items():
            fh.write(f'{k}\t{v}\n')
