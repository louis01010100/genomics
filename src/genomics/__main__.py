#!/usr/bin/env python

import sys
from argparse import ArgumentParser
from pathlib import Path


import polars as pl

from . import acmg, clinvar, dbsnp, vcf, depth, gene, variants, genome

__VERSION__ = '0.6.0'


def main():

    parser = config_parsers()

    args = parser.parse_args()

    if args.subcommand == 'clinvar':
        clinvar.process(
            clinvar_vcf_file=Path(args.clinvar_vcf_file),
            clinvar_papu_vcf_file=Path(args.clinvar_papu_vcf_file),
            genome_file=Path(args.genome_file),
            genome_index_file=Path(args.genome_index_file),
            output_dir=Path(args.output_dir),
        )
    elif args.subcommand == 'varmatch':
        variants.match(
            data1_file = Path(args.data1_file),
            data2_file = Path(args.data2_file),
            genome_file = Path(args.genome_file),
            output_file = Path(args.output_file),
            n_threads=args.n_threads,
            batch_size=args.batch_size,
        )

    elif args.subcommand == 'export-gene':
        gene.export(
            genes_file = _new_path(args.genes_file),
            gene_names_file = _new_path(args.gene_names_file),
            output_file = _new_path(args.output_file),
        )
    elif args.subcommand == 'dbsnp-normalize':
        dbsnp.normalize(
            dbsnp_vcf_file=Path(args.dbsnp_vcf_file),
            genome_file=Path(args.genome_file),
            genome_index_file=Path(args.genome_index_file),
            output_dir=Path(args.output_dir),
            n_threads=args.n_threads,
        )
    elif args.subcommand == 'dbsnp-merged2map':
        dbsnp.merged2map(
            input_json_file= Path(args.input_json_file),
            output_tsv_file = Path(args.output_tsv_file),
            n_threads = args.n_threads,
        )
    elif args.subcommand == 'dbsnp-createdb':
        dbsnp.create_db(
            dbsnp_vcf_file=Path(args.dbsnp_vcf_file),
            genome_file=Path(args.genome_file),
            output_dir=Path(args.output_dir),
            n_threads=args.n_threads,
        )
    elif args.subcommand == 'acmg':
        acmg.process(
            input_file=Path(args.input_file),
            output_file=Path(args.output_file),
        )
    elif args.subcommand == 'cram-depth':
        depth.export_cram_depths(
            crams_file=Path(args.crams_file),
            genders_file=Path(args.genders_file),
            genome_file=Path(args.genome_file),
            output_dir=Path(args.output_dir),
            n_threads=args.n_threads,
        )

    # elif args.subcommand == 'gvcf-depth':
    #     gvcf_files = _load_files(args.gvcfs_file, args.gvcf_files)
    #     depth.export_gvcf_depths(
    #         gvcf_files=gvcf_files,
    #         output_dir=Path(args.output_dir),
    #         coordinates_file=Path(args.coordinates_file),
    #         n_threads=args.n_threads,
    #     )
    else:
        parser.print_help()
        sys.exit(1)

    pass


def config_parsers():

    parser = ArgumentParser()

    parser.add_argument(
        '--version',
        action='version',
        version=f'%(prog)s {__VERSION__}',
    )

    parsers = parser.add_subparsers(dest='subcommand')

    _config_varmatch_parser(parsers.add_parser('varmatch'))
    _config_clinvar_parser(parsers.add_parser('clinvar'))
    _config_dbsnp_merged2map_parser(parsers.add_parser('dbsnp-merged2map'))
    _config_dbsnp_normalize_parser(parsers.add_parser('dbsnp-normalize'))
    _config_dbsnp_createdb_parser(parsers.add_parser('dbsnp-createdb'))
    _config_acmg_parser(parsers.add_parser('acmg'))
    _config_cram_depth_parser(parsers.add_parser('cram-depth'))
    _config_gvcf_depth_parser(parsers.add_parser('gvcf-depth'))
    _config_gene_parser(parsers.add_parser('export-gene'))

    return parser

def _config_varmatch_parser(parser):
    parser.add_argument('--data1-file', required=True, help='must contain the VCF 4-tuple')
    parser.add_argument('--data2-file', required=True, help='must be an indexed VCF file with the ID column populated')
    parser.add_argument('--genome-file', required=True)
    parser.add_argument('--output-file', required=True)
    parser.add_argument('--n-threads', type=int, default=1)
    parser.add_argument('--batch-size', type=int, default=1)

def _config_dbsnp_merged2map_parser(parser):
    parser.add_argument('--input-json-file', required=True)
    parser.add_argument('--output-tsv-file', required=True)
    parser.add_argument('--n-threads', type=int, default=1)

def _config_dbsnp_normalize_parser(parser):
    parser.add_argument('--dbsnp-vcf-file', required=True)
    parser.add_argument('--genome-file', required=True)
    parser.add_argument('--genome-index-file', required=True)
    parser.add_argument('--output-dir', required=True)
    parser.add_argument('--n-threads', type=int, default=1)

def _config_dbsnp_createdb_parser(parser):
    parser.add_argument('--dbsnp-vcf-file', required=True)
    parser.add_argument('--genome-file', required=True)
    parser.add_argument('--output-dir', required=True)
    parser.add_argument('--n-threads', type=int, default=1)


def _config_clinvar_parser(parser):
    parser.add_argument('--clinvar-vcf-file', required=True)
    parser.add_argument('--clinvar-papu-vcf-file', required=True)
    parser.add_argument('--genome-file', required=True)
    parser.add_argument('--genome-index-file', required=True)
    parser.add_argument('--output-dir', required=True)


def _config_acmg_parser(parser):
    parser.add_argument('--input-file', required=True)
    parser.add_argument('--output-file', required=True)


def _config_cram_depth_parser(parser):
    parser.add_argument('--crams-file', required = True)
    parser.add_argument('--genders-file', required = True)
    parser.add_argument('--output-dir', required=True)
    parser.add_argument('--genome-file', required=True)
    parser.add_argument('--n-threads', type=int, default=1)


def _config_gvcf_depth_parser(parser):
    parser.add_argument('--gvcfs-file')
    parser.add_argument('--output-dir', required=True)
    parser.add_argument('--coordinates-file', required=True)
    parser.add_argument('--n-threads', type=int, default=1)
    parser.add_argument('gvcf_files', nargs='*')


def _config_gene_parser(parser):
    parser.add_argument('--genes-file', required=True)
    parser.add_argument('--gene-names-file', required=False)
    parser.add_argument('--output-file', required=False)

def _new_path(file_):
    if file_:
        return Path(file_)
    return None


if __name__ == '__main__':

    main()
