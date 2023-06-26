#!/usr/bin/env python

import sys
from argparse import ArgumentParser
from pathlib import Path

from . import clinvar

__VERSION__ = '0.1.0'


def main():

    parser = config_parsers()

    args = parser.parse_args()

    if args.subcommand == 'clinvar':
        clinvar.process(
            clinvar_vcf_file=Path(args.clinvar_vcf_file),
            clinvar_papu_vcf_file=Path(args.clinvar_papu_vcf_file),
            genome_file=Path(args.genome_file),
            genome_index_file=Path(args.genome_index_file),
            tmp_dir=Path(args.tmp_dir),
            output_dir=Path(args.output_dir),
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

    return parser


def _config_clinvar_parser(parser):
    parser.add_argument('--clinvar-vcf-file', required=True)
    parser.add_argument('--clinvar-papu-vcf-file', required=True)
    parser.add_argument('--genome-file', required=True)
    parser.add_argument('--genome-index-file', required=True)
    parser.add_argument('--tmp-dir', required=True)
    parser.add_argument('--output-dir', required=True)


if __name__ == '__main__':
    main()
