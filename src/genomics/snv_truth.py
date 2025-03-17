import gzip
import shutil
from pathlib import Path
from subprocess import PIPE, STDOUT, Popen
from collections import OrderedDict
from icecream import ic

import numpy as np
import polars as pl
from pathos.multiprocessing import ProcessPool

from .utils import load, save, create_col2idx, is_gzip, load_dict, load_list, log_start, log_stop, log_info
from .variant import Variant, sync
from .vcf import Vcf, concat, fetch_variants, filter_variants, list_samples, list_contigs, standardize
from .gregion import GenomicRegion

COORDINATES_FILENAME = 'coordinates.tsv'

VCF_DIRNAME = 'vcf'

X_PAR_1 = GenomicRegion('chrX',10000, 2781479)
X_PAR_2 = GenomicRegion('chrX', 155701382, 156030895)

X_PARS = [X_PAR_1, X_PAR_2]


## Use male samples as the reference for read depth, as females are expected to have 0 depth on chrY.
# 10 hrs
def export_snv_truth(
    coordinates_file: Path,
    vcf_files: list,
    samples_file: Path,
    genders_file: Path,
    depths_file: Path,
    genome_file: Path,
    output_dir: Path,
    n_threads: int = 1,
    min_depth: int = 4,
    chrm_missing_as_homref: bool = False,
    merge_vcf: bool = False,
    prod: bool = True,
):

    if prod:
        shutil.rmtree(output_dir, ignore_errors=True)
    output_dir.mkdir(exist_ok=True)

    log_file = output_dir / 'snv_truth.log'

    info = OrderedDict()
    info['coordinates-file'] = coordinates_file
    info['samples-file'] = samples_file
    info['genders-file'] = genders_file
    info['depths-file'] = depths_file
    info['genome-file'] = genome_file
    info['output-dir'] = output_dir
    info['min-depth'] = min_depth
    info['chrm-missing-as-homref'] = chrm_missing_as_homref
    info['merge-vcf'] = merge_vcf
    info['n-threads'] = n_threads

    log_start(banner = 'SNV Truth Creation', info = info)

    vcf_tmp_dir = output_dir / VCF_DIRNAME

    vcf_tmp_dir.mkdir(exist_ok=True)

    log_info('load coordinates')
    coordinates = load_coordinates(coordinates_file)

    log_info('load depths')
    depths = load_depths(depths_file)

    log_info('load samples')
    samples = load_list(samples_file)
    sample2gender = load_dict(genders_file)

    log_info('subset samples')
    vcf_files = subset_samples(
        vcf_files=vcf_files,
        genome_file=genome_file,
        samples=samples,
        output_dir=vcf_tmp_dir,
        n_threads=n_threads,
    )

    log_info('fetch snvs')

    vcf_files = subset_snvs(
        coordinates = coordinates,
        vcf_files=vcf_files,
        sample2gender = sample2gender,
        genome_file = genome_file,
        depths=depths,
        min_depth=min_depth,
        output_dir=vcf_tmp_dir,
        chrm_missing_as_homref=chrm_missing_as_homref,
        n_threads=n_threads,
    )

    if merge_vcf:
        one_vcf = concat(
            vcf_files=vcf_files,
            output_file=output_dir / 'truth.vcf.bgz',
            tmp_dir=vcf_tmp_dir,
            n_threads=n_threads,
            preprocess=False,
        )
        print(one_vcf.filepath)
    else:
        for vcf_file in vcf_files:
            Vcf(vcf_file, vcf_tmp_dir).index().move_to(output_dir)

    print('done')


def subset_snvs(
    coordinates: list[dict],
    vcf_files: list[Path],
    sample2gender: dict,
    genome_file: Path,
    depths: dict,
    min_depth: int,
    output_dir: Path,
    chrm_missing_as_homref: bool,
    n_threads: int,
):

    genome = Genome(genome_file)

    bag = list()
    for vcf_file in vcf_files:

        output_file = output_dir / vcf_file.name.replace(
            '.vcf.bgz', '-snv.vcf')

        copy_header(vcf_file, output_file)

        samples = list_samples(vcf_file)
        n_samples = len(samples)

        genders = [sample2gender(sample) for sample in samples]


        with output_file.open('at') as fh:
            with ProcessPool(n_threads) as pool:
                for result in pool.uimap(
                        _subset_snvs, _subset_snv_jobs(
                            coordinates,
                            genome,
                            vcf_file,
                        )):

                    snv = fill_calls(
                        snv=result['snv'],
                        note=result['note'],
                        depths=depths,
                        genders = genders,
                        min_depth=min_depth,
                        chrm_missing_as_homref=
                        chrm_missing_as_homref,
                        genome = genome,
                    )

                    fh.write(f'{snv}\n')
                    n_done += 1
                    print(f'{n_done / n_total}', end='\r')

        bag.append(output_file)

    bag2 = list()
    with ProcessPool(n_threads) as pool:
        for result in pool.uimap(
                standardize_vcf,
                standardize_vcf_jobs(bag, output_dir),
        ):
            bag2.append(result)

    return bag2


def _subset_snv_jobs(coordinates, genome, vcf_file):
    contigs = list_contigs(vcf_file)

    for coordinate in coordinates:

        chrom = coordinate['chrom']
        pos = int(coordinate['pos']),
        id_ = coordinate['id']
        ref = coordinate['ref'],
        alt = coordinate['alt'],

        if id_ == '.':
            id_ = f'{chrom}_{pos}_{ref}_{alt}'

        snv = Variant(
                chrom = chrom,
                pos = int(coordinate['pos']),
                id_ =  id_,
                ref = ref,
                alt = alt,
        )

        yield {
                'snv': snv,
                'vcf_file': vcf_file,
                'chrom': genome.chromosome(chrom),
        }

def _subset_snvs(job):

    def _new_record(coordindate, coordinate_synced, variant, note):

        if variant is None:
            variant = Variant(
                chrom = '.', 
                pos = '.', 
                id_ = '.', 
                ref = '.', 
                alt = '.'
                qual = '.'
                filter_ = '.'
                info = '.'
            )
        return {
            'chrom': coordinate.chrom,
            'pos': coordinate.pos,
            'id': coordinate.id,
            'ref': coordinate.ref,
            'alt': coordinate.alt,
            'pos_synced': coordinate_synced.pos,
            'ref_synced': coordinate_synced.ref,
            'alt_synced': coordiante_synced.alt,
            'variant': variant,
            'note': note,
        }


    snv = job['snv']
    vcf_file = job['vcf_file']
    chromosome = job['chrom']

    max_region = snv.max_region(chromosome)
    candidates = fetch_variants(
        chrom = snv.chrom,
        pos = max_region.start,
        end = max_region.end,
        vcf_file = vcf_file,
        regions_overlap = 1,
    )

    if len(candidates) == 0:
        result = _new_record(coordinate, coordinate, None, 'missing')
    elif len(candidates) == 1:
        if '*' in candidates[0].alts:
            result = _new_record(coordinate, coordinate, None, 'complex')
        else:
            vx, vy = sync(coordinate, candidate, chromosome)
            if not result :
                result = refsnv
                note = 'complex'
            note = 'done'
    else:
        target = None
        candidates = [snv.align(refsnv, genome) for snv in candidates]

        for snv in candidates:
            if snv.same_coordinate(refsnv):
                target = snv
        if target:
            note = 'done'
        else:
            result = refsnv
            note = 'complex'
    result.id = id_

    return {'snv': result, 'note': note}

def _subset_snvs(job):

    refsnv = job['refsnv']
    span_start = job['span_start']
    span_end = job['span_end']
    vcf_file = job['vcf_file']
    genome = job['genome']


    snvs = fetch_variants(
        chrom = refsnv.chrom,
        pos = span_start,
        end = span_end,
        vcf_file = vcf_file,
        regions_overlap = 1,
    )

    if len(snvs) == 0:
        result = refsnv
        note = 'missing'
    elif len(snvs) == 1:
        if '*' in snvs[0].alts:
            result = refsnv
            note = 'complex'
        else:
            result = snvs[0].align(refsnv, genome)
            if not result :
                result = refsnv
                note = 'complex'
            note = 'done'
    else:
        target = None
        snvs = [snv.align(refsnv, genome) for snv in snvs]

        for snv in snvs:
            if snv.same_coordinate(refsnv):
                target = snv
        if target:
            note = 'done'
        else:
            result = refsnv
            note = 'complex'
    result.id = id_

    return {'snv': result, 'note': note}



def standardize_vcf(job):
    vcf_file = job['vcf_file']
    output_dir = job['output_dir']

    return standardize(vcf_file, output_dir)


def standardize_vcf_jobs(vcf_files, output_dir):
    for vcf_file in vcf_files:
        yield {
            'vcf_file': vcf_file,
            'output_dir': output_dir,
        }


def fill_calls(
        snv: Variant,
        note: str,
        depths: dict,
        genders: list,
        min_depth: int,
        chrm_missing_as_homref: bool,    #chrM always has high coverage
):


    if note == 'done':
        pass
    elif note == 'complex':
        snv.format = 'GT'
        snv.calls = create_dummy_calls(snv.chrom, snv.pos, '.', x_pars = X_PARS)
    elif note == 'missing':
        snv.format = 'GT'
        key = (snv.chrom, snv.pos)

        if not depths:
            snv.calls = create_dummy_calls(snv.chrom, snv.pos, '.', x_pars = X_PARS)
        elif chrm_missing_as_homref and 'm' in snv.chrom.lower():
            snv.calls = create_dummy_calls(snv.chrom, snv.pos, '0', x_pars = X_PARS)
        elif key not in depths:
            snv.calls = create_dummy_calls(snv.chrom, snv.pos, '.', x_pars = X_PARS)
        else:
            depth = depths[key]

            if depth >= min_depth:
                snv.calls = create_dummy_calls(snv.chrom, snv.pos, '0', x_pars = X_PARS)
            else:
                print(key, f'{depth} < {min_depth}')
                snv.calls = create_dummy_calls(snv.chrom, snv.pos, '.', x_pars = X_PARS)
    else:
        raise Exception(note)

    return snv


def copy_header(input_file, output_file):
    with gzip.open(input_file, 'rt') as ifh, output_file.open('wt') as ofh:
        for line in ifh:
            if line.startswith('##'):
                ofh.write(line)
                continue
            ofh.write(line)
            break


def subset_samples(vcf_files, genome_file, samples, output_dir, n_threads):

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
                    .uppercase() \
                    .exclude('ALT="."') \
                    .index()
        else:
            v = v.drop_info() \
                    .normalize(genome_file) \
                    .uppercase() \
                    .index()

        return v.filepath

    bag = []

    with ProcessPool(n_threads) as pool:
        for vcf_file in pool.uimap(
                process,
                jobs(vcf_files, output_dir, samples, trim_alts),
        ):
            bag.append(vcf_file)

    return bag



def load_depths(depths_file):

    def process(fh):
        bag = dict()
        col2idx = {
            column: idx
            for idx, column in enumerate((next(fh)).strip().split('\t'))
        }

        for line in fh:
            items = line.strip().split('\t')
            id_ = items[col2idx['chrom']]
            chrom = items[col2idx['chrom']]
            pos = int(items[col2idx['pos']])
            depth = items[col2idx['depth_min']]

            try:
                depth = float(depth)
            except ValueError:
                depth = np.nan

            key = (chrom, pos)

            bag[key] = depth
        return bag


    if not depths_file:
        return bag

    if is_gzip(depths_file):
        with gzip.open(depths_file, 'rt') as fh:
            return process(fh)
    else:
        with depths_file.open('rt') as fh:
            return process(fh)


def load_coordinates(coordiantes_file):
    coordinates = pl.read_csv(coordinates_file, comment_prefix = '##', has_header = True, separator = '\t')

    coordinates = coordinates.rename({'#CHROM': 'chrom'})
    coordinates.columns = [x.lowercase() for x in coordinates.columns]

    return coordinates.to_dicts()





def export_coordinates(input_file, output_dir):
    output_file = output_dir / COORDINATES_FILENAME
    coordinates.write_csv(output_file, has_header=False, separator='\t')

    return output_file

def create_dummy_calls(chrom: str, pos: int, genders: list, ref = True, x_pars = X_PARS):

    if ref:
        allele = '0'
    else:
        allele = '.'

    if 'chrX' == chrom:
        in_par = False
        for par in x_pars:
            if par.contains(chrom, pos):
                in_par = True

        if in_par:
            return '\t'.join([f'{allele}/{allele}' for x in genders])
        else:
            return '\t'.join([f'{allele}' if x == 'male' else f'{allele}/{allele}' for x in genders])

    elif 'chrY' == chrom:
        return '\t'.join([f'{allele}' if x == 'male' else f'.' for x in genders])

    elif 'chrM' ==  chrom:
        return '\t'.join([f'{allele}' for x in genders])
    else:
        return '\t'.join([f'{allele}/{allele}' for x in genders])

