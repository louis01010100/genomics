import gzip
import shutil
from pathlib import Path
from subprocess import PIPE, STDOUT, Popen
from collections import OrderedDict
from icecream import ic
from .genome import Genome

import numpy as np
import polars as pl
from pathos.multiprocessing import ProcessPool

from .utils import is_gzip, load_dict, load_list, init_logging, log_start, log_stop, log_info,copy_vcf_header
from .variant import Variant, sync
from .vcf import Vcf, concat, fetch_variants, filter_variants, list_samples, list_contigs, standardize
from .gregion import GenomicRegion
from .genome import Genome

COORDINATES_FILENAME = 'coordinates.tsv'

SUBSET_SAMPLES_DIR = 'samples'
SUBSET_SNVS_DIR = 'snv'

X_PAR_1 = GenomicRegion('chrX',10000, 2781479)
X_PAR_2 = GenomicRegion('chrX', 155701382, 156030895)

X_PARS = [X_PAR_1, X_PAR_2]

JOB_SIZE=10000


## Use male samples as the reference for read depth, as females are expected to have 0 depth on chrY.
# 10 hrs
def export_snv_truth(
    coordinates_file: Path,
    vcf_files: list,
    depths_file: Path,
    samples_file: Path,
    genders_file: Path,
    genome_file: Path,
    output_dir: Path,
    n_threads: int = 1,
    min_depth: int = 4,
    chrm_missing_as_homref: bool = False,
    merge_vcf: bool = False,
    prod: bool = True,
):

    # if prod:
    #     shutil.rmtree(output_dir, ignore_errors=True)
    output_dir.mkdir(exist_ok=True)

    log_file = output_dir / 'snv_truth.log'
    init_logging(log_file)

    info = OrderedDict()
    info['coordiantes-file'] = coordinates_file
    info['n_vcf_files'] = len(vcf_files)
    info['samples-file'] = samples_file
    info['genders-file'] = genders_file
    info['genome-file'] = genome_file
    info['output-dir'] = output_dir
    info['min-depth'] = min_depth
    info['chrm-missing-as-homref'] = chrm_missing_as_homref
    info['merge-vcf'] = merge_vcf
    info['n-threads'] = n_threads

    log_start(banner = 'SNV Truth Creation', info = info)

    sample_vcf_dir = output_dir / 'samples'
    sample_vcf_dir.mkdir(exist_ok=True)

    log_info('load samples')
    samples = load_list(samples_file)
    sample2gender = load_dict(genders_file)

    log_info('subset samples')
    vcf_files = subset_samples(
        vcf_files=vcf_files,
        genome_file=genome_file,
        samples=samples,
        trim_alts = True,
        output_dir=sample_vcf_dir,
        n_threads=n_threads,
    )

    vcf_files = [x for x in sample_vcf_dir.glob('*/*info-norm-uppercase-samples-format-info-trim_alt-norm-uppercase-ex.vcf.bgz')]

    # log_info('groupd spanning deletions')
    #
    # spandel_dir = output_dir / 'spandels'
    # spandel_dir.mkdir(exist_ok=True)
    #
    # vcf_files = group_spanning_deletions(
    #         vcf_files = vcf_files,
    #         output_dir = spandel_dir,
    #         n_threads = n_threads,
    # )

    log_info('load coordinates')
    coordinates = load_coordinates(coordinates_file)

    log_info('load depths')
    depths = load_depths(depths_file)

    log_info('fetch snvs')
    snv_vcf_dir = output_dir / 'snv'

    result = subset_snvs(
        coordinates = coordinates,
        vcf_files=vcf_files,
        sample2gender = sample2gender,
        genome_file = genome_file,
        depths=depths,
        min_depth=min_depth,
        output_dir=snv_vcf_dir,
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

    for vcf_file in vcf_files:

        task_dir = output_dir / vcf_file.name.replace('.vcf.bgz', '')
        task_dir.mkdir(parents = True, exist_ok = True)

        truth_vcf_file = task_dir / 'truth.vcf'

        truth_synced_vcf_file = task_dir / 'truth_synced.vcf'

        snv_profile_file = task_dir / 'snv_profile.tsv'

        copy_vcf_header(vcf_file, truth_vcf_file)
        copy_vcf_header(vcf_file, truth_synced_vcf_file)

        with snv_profile_file.open('wt') as fh:
            snv_profile_header = '\t'.join([
                    'chrom',
                    'pos',
                    'id',
                    'ref',
                    'alt',
                    'pos_synced',
                    'ref_synced',
                    'alt_synced',
                    'note',
            ])
            fh.write(f'{snv_profile_header}\n')

        samples = list_samples(vcf_file)
        n_samples = len(samples)

        genders = [sample2gender[sample] for sample in samples]

        with truth_vcf_file.open('at') as fh_t, \
                truth_synced_vcf_file.open('at') as fh_s, \
                snv_profile_file.open('at') as fh_p:

            # for job in fetch_calls_jobs(
            #             coordinates = coordinates,
            #             genome = genome,
            #             depths = depths,
            #             genders = genders,
            #             min_depth = min_depth,
            #             chrm_missing_as_homref = chrm_missing_as_homref,
            #             vcf_file = vcf_file,
            #         ):
            #
            #     results = fetch_calls(job)
            #     for result in results:
            #
            #         variant= result['variant']
            #         fh_t.write(f'{variant}\n')
            #
            #         variant_synced = result['variant_synced']
            #         fh_s.write(f'{variant_synced}\n')
            #
            #         record = '\t'.join([
            #                 result['chrom'],
            #                 str(result['pos']),
            #                 result['id'],
            #                 result['ref'],
            #                 result['alt'],
            #                 str(result['pos_synced']),
            #                 result['ref_synced'],
            #                 result['alt_synced'],
            #                 result['note'],
            #         ])
            #         fh_p.write(f'{record}\n')
            #

            with ProcessPool(n_threads) as pool:
                for results in pool.uimap(
                        fetch_calls, fetch_calls_jobs(
                            coordinates = coordinates,
                            genome = genome,
                            depths = depths,
                            genders = genders,
                            min_depth = min_depth,
                            chrm_missing_as_homref = chrm_missing_as_homref,
                            vcf_file = vcf_file,
                        )):
                    for result in results:

                        variant= result['variant']
                        fh_t.write(f'{variant}\n')

                        variant_synced = result['variant_synced']
                        fh_s.write(f'{variant_synced}\n')

                        record = '\t'.join([
                                result['chrom'],
                                str(result['pos']),
                                result['id'],
                                result['ref'],
                                result['alt'],
                                str(result['pos_synced']),
                                result['ref_synced'],
                                result['alt_synced'],
                                result['note'],
                        ])
                        fh_p.write(f'{record}\n')

    return {
        'truth_vcf_files' : list(task_dir.glob('*/truth.vcf')),
        'truth_synced_vcf_files' : list(task_dir.glob('*/truth_synced.vcf')),
        'snv_profile_files' : list(task_dir.glob('*/snv_profile.tsv')),
    }


def format_report(data):
        return {
            'chrom': coordinate.chrom,
            'pos': coordinate.pos,
            'id': coordinate.id,
            'ref': coordinate.ref,
            'alt': coordinate.alt,
            'pos_synced': coordinate_synced.pos,
            'ref_synced': coordinate_synced.ref,
            'alt_synced': coordinate_synced.alt,
            'variant': variant,
            'variant_synced': variant_synced,
            'note': note,
        }


def fetch_calls_jobs(coordinates, genome, depths, genders, min_depth, chrm_missing_as_homref, vcf_file):
    contigs = list_contigs(vcf_file)

    bag = list()

    for coordinate in coordinates:

        chrom = coordinate['chrom']

        if chrom not in contigs:
            continue

        pos = int(coordinate['pos']),
        id_ = coordinate['id']
        ref = coordinate['ref']
        alt = coordinate['alt']

        if id_ == '.':
            id_ = f'{chrom}_{pos}_{ref}_{alt}'

        coordinate = Variant(
                chrom = chrom,
                pos = int(coordinate['pos']),
                id_ =  id_,
                ref = ref,
                alt = alt,
        )

        job =  {
                'coordinate': coordinate,
                'vcf_file': vcf_file,
                'chrom': genome.chromosome(chrom),
                'depths': depths,
                'genders': genders,
                'min_depth': min_depth,
                'chrm_missing_as_homref' : chrm_missing_as_homref,
        }

        bag.append(job)

        if len(bag) > JOB_SIZE:
            yield bag
            bag = list()

    if len(bag) > 0:
        yield bag


def fetch_calls(jobs):

    def _new_record(coordinate, coordinate_synced, variant, variant_synced, note):

        if variant is None:
            variant = coordinate.clone()
            variant_synced = coordinate_synced.clone()
        return {
            'chrom': coordinate.chrom,
            'pos': coordinate.pos,
            'id': coordinate.id,
            'ref': coordinate.ref,
            'alt': coordinate.alt,
            'pos_synced': coordinate_synced.pos,
            'ref_synced': coordinate_synced.ref,
            'alt_synced': coordinate_synced.alt,
            'variant': variant,
            'variant_synced': variant_synced,
            'note': note,
        }


    bag = list()

    for job in jobs:
        coordinate = job['coordinate']
        vcf_file = job['vcf_file']
        chromosome = job['chrom']
        genders=job['genders']
        depths=job['depths']
        min_depth=job['min_depth']
        chrm_missing_as_homref=job['chrm_missing_as_homref']

        coordinate_expanded = coordinate.expand(chromosome)

        candidates = fetch_variants(
            chrom = coordinate_expanded.chrom,
            pos = coordinate_expanded.region.start,
            end = coordinate_expanded.region.end,
            vcf_file = vcf_file,
            regions_overlap = 1,
        )

        if len(candidates) == 0:
            record = _new_record(coordinate, coordinate_expanded, None, None, 'MISSING')
        else:
            assert  len([ x for x in candidates if '*' in x.alts]) == 0

            if len(candidates) == 0:
                record = _new_record(coordinate, coordinate_expanded, None, None, 'COMPLEX')
            elif len(candidates) == 1:
                candidate = candidates[0]
                coordinate_synced, candidate_synced = sync(coordinate, candidate, chromosome)
                record  = _new_record(coordinate, coordinate_synced, candidate, candidate_synced, 'DONE')
            else:
                assert len(candidates) > 1

                best_coordinate_synced = None
                best_candidate = None
                best_candidate_synced = None
                n_match_alts = None

                for candidate in candidates:
                    coordinate_synced, candidate_synced = sync(coordinate, candidate, chromosome)

                    n_match_alts_current = len(set(coordinate_synced.alts) & set(candidate_synced.alts))

                    if n_match_alts is None:
                        n_match_alts = n_match_alts_current
                        best_coordinate_synced = coordinate_synced
                        best_candidate = candidate
                        best_candidate_synced = candidate_synced
                    elif n_match_alts_current > n_match_alts:
                        n_match_alts = n_match_alts_current
                        best_coordinate_synced = coordinate_synced
                        best_candidate = candidate
                        best_candidate_synced = candidate_synced
                    else:
                        continue
                record = _new_record(coordinate, best_coordinate_synced, best_candidate, best_candidate_synced, 'DONE')

        record = fill_missing_calls(
            record,
            depths=depths,
            genders = genders,
            min_depth=min_depth,
            chrm_missing_as_homref=
            chrm_missing_as_homref,
        )

        bag.append(record)

    return bag


def fill_missing_calls(
        data: dict,
        depths: dict,
        genders: list,
        min_depth: int,
        chrm_missing_as_homref: bool,    #chrM always has high coverage
):

    # 'chrom': coordinate.chrom,
    # 'pos': coordinate.pos,
    # 'id': coordinate.id,
    # 'ref': coordinate.ref,
    # 'alt': coordinate.alt,
    # 'pos_synced': coordinate_synced.pos,
    # 'ref_synced': coordinate_synced.ref,
    # 'alt_synced': coordinate_synced.alt,
    # 'variant': variant,
    # 'variant_synced': variant_synced,
    # 'note': note,

    note = data['note']

    if note == 'DONE':
        return data

    chrom = data['chrom']
    pos = data['pos']
    ref = data['ref']

    pos_synced = data['pos_synced']
    ref_synced = data['ref_synced']

    region = GenomicRegion(chrom, pos, pos + len(ref) - 1)
    region_synced = GenomicRegion(chrom, pos_synced, pos_synced + len(ref_synced) - 1)

    # if data['id'] == 'AX-656564988':
    #     print(note, region_synced)


    if note == 'COMPLEX':
        ref = False

    elif note == 'MISSING':

        keys = [(chrom, pos) for pos in range(region_synced.start, region_synced.end + 1)]


        if not depths:
            ref = False
        elif len(keys & depths.keys()) == 0:
            ref = False
        elif chrm_missing_as_homref and 'm' in chrom.lower():
            ref = True
        else:
            bag = list()
            for key in keys:
                if key not in depths:
                    assert False, f'{key};{data}'

            depth = min([depths[key] for key in keys])

            if depth >= min_depth:
                ref = True
                note = 'DONE'
            else:
                ref = False
                print(key, f'{depth} < {min_depth}')

        # if data['id'] == 'AX-656564988':
        #     print(keys, region_synced, depth)
    else:
        raise Exception(note)

    data['note'] = note
    data['variant'].format = 'GT'
    data['variant_synced'].format = 'GT'
    data['variant'].calls = create_dummy_calls(chrom,  region, ref = ref, genders = genders, x_pars = X_PARS)
    data['variant_synced'].calls = create_dummy_calls(chrom, region_synced, ref = ref, genders = genders, x_pars = X_PARS)

    return data


def subset_samples(vcf_files, genome_file, samples, output_dir, trim_alts, n_threads):

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
                    .subset_samples(samples, ) \
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


def load_coordinates(coordinates_file):
    coordinates = pl.read_csv(coordinates_file, comment_prefix = '##', has_header = True, separator = '\t')

    coordinates = coordinates.rename({'#CHROM': 'chrom'})
    coordinates.columns = [x.lower() for x in coordinates.columns]

    return coordinates.to_dicts()





def export_coordinates(input_file, output_dir):
    output_file = output_dir / COORDINATES_FILENAME
    coordinates.write_csv(output_file, has_header=False, separator='\t')

    return output_file

def create_dummy_calls(chrom: str, region: GenomicRegion, genders: list, ref = True, x_pars = X_PARS):

    if ref:
        allele = '0'
    else:
        allele = '.'

    if 'chrX' == chrom:
        in_par = False
        for par in x_pars:
            if par.overlaps(region):
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



def group_spanning_deletions(
            vcf_files,
            output_dir,
            n_threads,
    ):

    def process(job):
        vcf_file = job['vcf_file']
        output_dir = job['output_dir']
        return Vcf(vcf_file, tmp_dir = output_dir).group_spanning_deletions().filepath

    def jobs(vcf_files, output_dir):
        for vcf_file in vcf_files:
            print(vcf_file)
            yield {
                'vcf_file': vcf_file,
                'output_dir': output_dir,
            }


    with ProcessPool(n_threads) as pool:
        for future in pool.uimap(process, jobs(vcf_files, output_dir)):
            pass

    return [x for x in output_dir.glob('**/*.vcf.bgz')]


