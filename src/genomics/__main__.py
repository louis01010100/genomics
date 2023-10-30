#!/usr/bin/env python

import sys
from argparse import ArgumentParser
from pathlib import Path
import polars as pl

from . import acmg, clinvar, dbsnp, truth, vcf, coordinate, depth

__VERSION__ = '0.3.0'


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
    elif args.subcommand == 'coordinate':
        coordinate.export_coordinates(
            input_vcf_file=Path(args.input_vcf_file),
            genome_file=Path(args.genome_file),
            genome_index_file=Path(args.genome_index_file),
            output_dir=Path(args.output_dir),
        )
    elif args.subcommand == 'cram-depth':

        cram_files = _load_files(args.crams_file, args.cram_files)

        depth.export_cram_depths(
            cram_files=cram_files,
            genome_file=Path(args.genome_file),
            output_dir=Path(args.output_dir),
            coordinates_vcf_file=Path(args.coordinates_vcf_file),
            n_threads=args.n_threads,
        )

    elif args.subcommand == 'gvcf-depth':
        gvcf_files = _load_files(args.gvcfs_file, args.gvcf_files)
        depth.export_gvcf_depths(
            gvcf_files=gvcf_files,
            output_dir=Path(args.output_dir),
            coordinates_vcf_file=Path(args.coordinates_vcf_file),
            n_threads=args.n_threads,
        )
    elif args.subcommand == 'snv-truth':
        truth.export_snv_truth(
            vcfs_file=Path(args.vcfs_file),
            samples_file=Path(args.samples_file),
            coordinates_vcf_file=Path(args.coordinates_vcf_file),
            depths_file=Path(args.depths_file),
            genome_file=Path(args.genome_file),
            output_dir=Path(args.output_dir),
            min_depth=args.min_depth,
            explode=args.explode,
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
    _config_coordinate_parser(parsers.add_parser('coordinate'))
    _config_cram_depth_parser(parsers.add_parser('cram-depth'))
    _config_gvcf_depth_parser(parsers.add_parser('gvcf-depth'))
    _config_snv_truth_parser(parsers.add_parser('snv-truth'))

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


def _config_coordinate_parser(parser):
    parser.add_argument('--input-vcf-file', required=True)
    parser.add_argument('--genome-file', required=True)
    parser.add_argument('--genome-index-file', required=True)
    parser.add_argument('--output-dir', required=True)


def _config_cram_depth_parser(parser):
    parser.add_argument('--crams-file')
    parser.add_argument('--output-dir', required=True)
    parser.add_argument('--genome-file', required=True)
    parser.add_argument('--coordinates-vcf-file', required=True)
    parser.add_argument('--n-threads', type=int, default=1)
    parser.add_argument('cram_files', nargs='*')


def _config_gvcf_depth_parser(parser):
    parser.add_argument('--gvcfs-file')
    parser.add_argument('--output-dir', required=True)
    parser.add_argument('--coordinates-vcf-file', required=True)
    parser.add_argument('--n-threads', type=int, default=1)
    parser.add_argument('gvcf_files', nargs='*')


def _config_snv_truth_parser(parser):
    parser.add_argument('--vcfs-file', required=False)
    parser.add_argument('--samples-file', required=True)
    parser.add_argument('--coordinates-vcf-file', required=True)
    parser.add_argument('--depths-file', required=True)
    parser.add_argument('--genome-file', required=True)
    parser.add_argument('--output-dir', required=True)
    parser.add_argument('--min-depth', type=int, default=2)
    parser.add_argument('--explode', type=bool, default=True)
    parser.add_argument('--n-threads', type=int, default=1)


def _load_files(manifest_file, files):

    if manifest_file:
        files = list()
        with open(manifest_file, 'rt') as fh:
            next(fh)
            for line in fh:
                files.append(line.strip())

    files = list([Path(x) for x in files])

    return files


if __name__ == '__main__':

    main()
