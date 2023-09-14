#!/usr/bin/env python

import sys
from argparse import ArgumentParser
from pathlib import Path

from . import acmg, clinvar, dbsnp, kgp, vcf

__VERSION__ = '0.2.0'


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
    elif args.subcommand == 'dbsnp':
        dbsnp.process(
            dbsnp_vcf_file=Path(args.dbsnp_vcf_file),
            genome_file=Path(args.genome_file),
            genome_index_file=Path(args.genome_index_file),
            output_dir=Path(args.output_dir),
            n_threads=args.n_threads,
        )
    elif args.subcommand == 'acmg':
        acmg.process(
            input_file=Path(args.input_file),
            output_file=Path(args.output_file),
        )
    elif args.subcommand == 'kgp':
        kgp.process(
            vcf_dir=Path(args.vcf_dir),
            cram_dir=Path(args.cram_dir),
            ref_file=Path(args.ref_file),
            output_dir=Path(args.output_dir),
            samples_file=Path(args.samples_file),
            coordinate_vcf_file=Path(args.coordinate_vcf_file),
            n_cram_samples=args.n_cram_samples,
            min_read_depth=args.min_read_depth,
            n_threads=args.n_threads,
        )
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

    _config_clinvar_parser(parsers.add_parser('clinvar'))
    _config_dbsnp_parser(parsers.add_parser('dbsnp'))
    _config_acmg_parser(parsers.add_parser('acmg'))
    _config_kgp_parser(parsers.add_parser('kgp'))

    return parser


def _config_dbsnp_parser(parser):
    parser.add_argument('--dbsnp-vcf-file', required=True)
    parser.add_argument('--genome-file', required=True)
    parser.add_argument('--genome-index-file', required=True)
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


def _config_kgp_parser(parser):
    parser.add_argument('--vcf-dir', required=True)
    parser.add_argument('--cram-dir', required=True)
    parser.add_argument('--output-dir', required=True)
    parser.add_argument('--ref-file', required=True)
    parser.add_argument('--samples-file', required=True)
    parser.add_argument('--coordinate-vcf-file', required=True)
    parser.add_argument('--min-read-depth', type=int, default=2)
    parser.add_argument('--n-cram-samples', type=int, default=10)
    parser.add_argument('--n-threads', type=int, default=1)


if __name__ == '__main__':
    main()
