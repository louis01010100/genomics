import gzip
import shutil
from pathlib import Path
from subprocess import PIPE, STDOUT, Popen
import logging
from icecream import ic

import numpy as np
import polars as pl
from pathos.multiprocessing import ProcessPool

from .utils import load, save
from .variant import Variant
from .vcf import Vcf, concat, fetch_variants, filter_variants, list_samples

COORDINATES_FILENAME = 'coordinates.tsv'

TRUTH_VCF_FILE = 'truth.vcf.bgz'
VCF_DIRNAME = 'vcf'


def export_snv_truth(
    vcf_files: list,
    samples_file: Path,
    coordinates_vcf_file: Path,
    depths_file: Path,
    genome_file: Path,
    output_dir: Path,
    n_threads: int = 1,
    min_depth: int = 4,
    prod: bool = False,
):
    if prod:
        shutil.rmtree(output_dir, ignore_errors=True)

    output_dir.mkdir(exist_ok=True)

    vcf_tmp_dir = output_dir / VCF_DIRNAME

    vcf_tmp_dir.mkdir(exist_ok=True)

    print('load coordinates')
    coordinates_file = export_coordinates(coordinates_vcf_file, output_dir)

    print('load depths')
    depths = load_depths(depths_file)

    print('load samples')
    samples = load_samples(samples_file)

    print('subset samples')
    vcf_file = subset_samples(
        vcf_files=vcf_files,
        genome_file=genome_file,
        samples=samples,
        output_dir=vcf_tmp_dir,
        trim_alts=len(depths) > 0,
        n_threads=n_threads,
    )

    vcf_file = vcf_tmp_dir / 'samples.vcf.bgz'

    coordinates = pl.read_csv(
        coordinates_file,
        has_header=False,
        new_columns=[
            'chrom', 'pos', 'id', 'ref', 'alt', 'qual', 'filter', 'info'
        ],
        separator='\t',
    )

    print('fetch snvs')
    result_vcf = subset_snvs(
        vcf_file=vcf_file,
        coordinates=coordinates,
        depths=depths,
        min_depth=min_depth,
        output_dir=vcf_tmp_dir,
        n_threads=n_threads,
    )

    kgp_truth_vcf = result_vcf.move_to(output_dir / TRUTH_VCF_FILE)
    print('done')


def subset_snvs(
    vcf_file: Path,
    coordinates: pl.DataFrame,
    depths: dict,
    min_depth: int,
    output_dir: Path,
    n_threads: int,
):

    def jobs(coordinates, vcf_file):
        for record in coordinates.to_dicts():
            yield {
                'chrom': record['chrom'],
                'pos': int(record['pos']),
                'id': record['id'],
                'ref': record['ref'],
                'alt': record['alt'],
                'vcf_file': vcf_file,
            }

    def process(job):
        chrom = job['chrom']
        pos = int(job['pos'])
        id_ = job['id']
        ref = job['ref']
        alt = job['alt']
        vcf_file = job['vcf_file']

        refsnv = Variant(chrom=chrom, pos=pos, ref=ref, alt=alt, id_=id_)

        snvs = fetch_variants(
            refsnv.chrom,
            refsnv.pos,
            vcf_file,
            regions_overlap=1,
        )

        if len(snvs) == 0:
            result = refsnv
            note = 'missing'
        elif len(snvs) == 1:
            if '*' in snvs[0].alts:
                result = refsnv
                note = 'complex'
            else:
                note = 'done'
                result = refsnv.sync_alleles(snvs[0])
        else:
            target = None
            for snv in snvs:
                if '*' in snv.alts:
                    continue
                if snv.pos == refsnv.pos:
                    target = snv
            if target:
                result = refsnv.sync_alleles(target)
                note = 'done'
            else:
                result = refsnv
                note = 'complex'

        return {'snv': result, 'note': note}

    n_samples = len(list_samples(vcf_file))

    output_file = output_dir / vcf_file.name.replace('.vcf.bgz', '-snv.vcf')

    copy_header(vcf_file, output_file)

    n_total = len(coordinates)
    n_done = 0
    with output_file.open('at') as fh:
        with ProcessPool(n_threads) as pool:
            for result in pool.uimap(process, jobs(coordinates, vcf_file)):

                for variant in new_variants(
                        snv=result['snv'],
                        note=result['note'],
                        depths=depths,
                        n_samples=n_samples,
                        min_depth=min_depth,
                ):
                    fh.write(f'{variant}\n')
                n_done += 1
                print(f'{n_done / n_total}', end='\r')

    return Vcf(output_file, output_dir, n_threads).bgzip().drop_qual(
    ).drop_filter().drop_info().fill_tags().sort().index()


def new_variants(snv: Variant, note: str, depths: dict, n_samples: int,
                 min_depth: int):

    if note == 'done':
        pass
    elif note == 'complex':
        snv.format = 'GT'
        snv.calls = '\t'.join(['./.' for i in range(0, n_samples)])
    elif note == 'missing':
        snv.format = 'GT'
        key = (snv.chrom, snv.pos)

        if not depths:
            snv.calls = '\t'.join(['./.' for i in range(0, n_samples)])
        elif key not in depths:
            snv.calls = '\t'.join(['./.' for i in range(0, n_samples)])
        else:
            depth = depths[key]

            if depth >= min_depth:
                snv.calls = '\t'.join(['0/0' for i in range(0, n_samples)])
            else:
                print(key, f'{depth} < {min_depth}')
                snv.calls = '\t'.join(['./.' for i in range(0, n_samples)])
    else:
        raise Exception(note)

    variants = list()

    for id_ in snv.id.split(','):
        variants.append(
            Variant(
                chrom=snv.chrom,
                pos=snv.pos,
                ref=snv.ref,
                alt=snv.alt,
                id_=id_,
                qual=snv.qual,
                filter_=snv.filter,
                info=snv.info,
                format_=snv.format,
                calls=snv.calls,
            ))

    return variants


def copy_header(input_file, output_file):
    with gzip.open(input_file, 'rt') as ifh, output_file.open('wt') as ofh:
        for line in ifh:
            if line.startswith('##'):
                ofh.write(line)
                continue
            ofh.write(line)
            break


def subset_samples(vcf_files, genome_file, samples, output_dir, trim_alts,
                   n_threads):

    def jobs(vcf_files, outupt_dir, samples, trim_alts):
        for vcf_file in vcf_files:
            yield {
                'vcf_file': vcf_file,
                'output_dir': output_dir,
                'samples': samples,
                'trim_alts': trim_alts,
            }

    def process(job):
        vcf_file = job['vcf_file']
        output_dir = job['output_dir']
        samples = job['samples']
        trim_alts = job['trim_alts']
        if samples:
            v = Vcf(vcf_file, output_dir) \
                    .subset_samples(samples) \
                    .keep_format(fields = ['GT'])
        else:
            v = Vcf(vcf_file, output_dir) \
                    .keep_format(fields = ['GT'])

        if trim_alts:
            v = v.drop_info() \
                    .trim_alts() \
                    .normalize(genome_file) \
                    .exclude('ALT="."') \
                    .index()
        else:
            v = v.drop_info() \
                    .normalize(genome_file) \
                    .index()

        return v.filepath

    bag = []

    with ProcessPool(n_threads) as pool:
        for vcf_file in pool.uimap(
                process,
                jobs(vcf_files, output_dir, samples, trim_alts),
        ):
            bag.append(vcf_file)

    output_file = output_dir / 'samples.vcf.bgz'

    result_vcf = concat(
        vcf_files=bag,
        output_file=output_file,
        tmp_dir=output_dir,
        n_threads=n_threads,
    )

    return output_file


def load_samples(file_):
    if not file_:
        return set()
    samples = set(
        pl.read_csv(
            file_,
            has_header=True,
            separator='\t',
        )['sample_id'])


def load_depths(depths_file):
    bag = dict()

    if not depths_file:
        return bag

    with depths_file.open('rt') as fh:
        col2idx = {
            column: idx
            for idx, column in enumerate((next(fh)).strip().split('\t'))
        }

        for line in fh:
            items = line.strip().split('\t')
            chrom = items[col2idx['chrom']]
            pos = int(items[col2idx['pos']])
            depth = items[col2idx['depth']]

            try:
                depth = float(depth)
            except ValueError:
                depth = np.nan

            key = (chrom, pos)

            bag[key] = depth

    return bag


def export_coordinates(input_file, output_dir):
    output_file = output_dir / COORDINATES_FILENAME
    coordinates = Vcf(input_file, output_dir).to_df(site_only=True)
    coordinates.write_csv(output_file, has_header=False, separator='\t')

    return output_file
