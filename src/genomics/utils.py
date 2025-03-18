import gzip
import logging
import humanfriendly
import pickle
import shutil
from pathlib import Path
import polars as pl
from subprocess import PIPE, STDOUT, Popen
from typing import TextIO, Tuple, Union
import sys

import pandas as pd
from icecream import ic

COMPLEMENT_BASES = str.maketrans("ACGT", "TGCA")
BANNER_WIDTH = 50


class _AllelePairs():

    def __init__(self):
        self.allele_pairs = dict()

    def update(self, other) -> None:
        self.allele_pairs.update(other.allele_pairs)

    def add_allele_pair(self, ref, alt) -> None:
        key = (ref, alt)
        self.allele_pairs[key] = {'ref': ref, 'alt': alt}

        longest_ref = max([x['ref'] for x in self.allele_pairs.values()],
                          key=len)

        delta = dict()

        for key, value in self.allele_pairs.items():
            ref = value['ref']
            alt = value['alt']

            if ref == longest_ref:
                continue

            suffix = longest_ref[len(ref):]

            new_ref = ref + suffix
            new_alt = ','.join([x + suffix for x in alt.split(',')])

            delta[key] = {'ref': new_ref, 'alt': new_alt}

        self.allele_pairs.update(delta)

    def get_updated_allele_pair(self, ref, alt) -> dict:
        key = (ref, alt)
        return self.allele_pairs[key]

    def __str__(self) -> str:
        return __repr__(self)

    def __repr__(self) -> str:

        bag = []

        for k, v in self.allele_pairs.items():
            value = '{' + f'{k}:{v}' + '}'
            bag.append(value)
            output = ','.join(bag)

        return f'[{output}]'

def init_logging(log_file: Path):

    h0 = logging.StreamHandler(sys.stderr)
    h0.setLevel(logging.INFO)

    h1 = logging.FileHandler(log_file, 'w')
    h1.setLevel(logging.DEBUG)

    logging.basicConfig(
        level=logging.DEBUG,
        handlers=[h0, h1],
        format='[%(levelname)s] %(asctime)s\t%(message)s',
        datefmt='%Y-%m-%d %H:%M:%S',
    )

def log_info(message: str):
    logging.info('')
    logging.info(message)

def log_start(
    banner: str,
    info: dict,
):
    logging.info('#' * BANNER_WIDTH)
    logging.info(banner.center(BANNER_WIDTH))
    logging.info('#' * BANNER_WIDTH)

    for k, v in info.items():
        logging.info(f"{k}: {v}")


def log_stop(banner: str, start_time, stop_time):
    logging.info('')
    logging.info('#' * BANNER_WIDTH)
    logging.info(banner.center(BANNER_WIDTH))
    logging.info('#' * BANNER_WIDTH)
    logging.info(f"start time: {start_time}")
    logging.info(f"stop time: {stop_time}")
    logging.info(
        f"duration: {humanfriendly.format_timespan(stop_time - start_time)}")


def load_list(file_, dedup = True):

    samples = pl.read_csv(
            file_,
            has_header=True,
            separator='\t',
        )['sample']

    if dedup:
        samples = list(set(samples))
    else:
        samples = list(samples)

    return samples

def load_dict(file_):
    bag = dict()
    with file_.open('rt') as fh:
        next(fh)

        for line in fh:
            items = line.strip().split('\t')
            bag[items[0]] = items[1]

    return bag

def is_gzip(filepath):
    with open(filepath, 'rb') as fd:
        magic_number = fd.read(2)
        if magic_number == b'\x1f\x8b':
            return True
        return False


def vcf2dict(*vcf_files) -> dict:
    result = None

    for vcf_file in vcf_files:
        if is_gzip(vcf_file):
            with gzip.open(vcf_file, 'rt') as fh:
                output = _vcf2dict(fh)
        else:
            with vcf_file.open('rt') as fh:
                output = _vcf2dict(fh)

        if not result:
            result = output

        for coordinate, allele_pairs in output.items():
            if coordinate not in result:
                continue
            result[coordinate].update(allele_pairs)

    return result


def _vcf2dict(fh: TextIO) -> dict:
    bag = dict()
    for line in fh:
        if line.startswith('#'):
            continue

        chrom, pos, id_, ref, alt, rest = line.strip().split('\t', 5)

        coordinate = (chrom, pos)

        if coordinate not in bag:
            bag[coordinate] = _AllelePairs()
        bag[coordinate].add_allele_pair(ref, alt)

    return bag


def index(filepath: Path):
    cmd = f'bcftools index {filepath}'
    execute(cmd)


def subset(
    input_vcf_file: Path,
    output_vcf_file: Path,
    regions_file: Path,
):
    cmd = f'bcftools view --regions-file {regions_file} {input_vcf_file} -O z -o {output_vcf_file}'
    execute(cmd)


def bgzip(filepath: Path, index=True, n_threads=1):
    cmd = f'bgzip  --threads {n_threads}   {filepath}'
    execute(cmd)

    if index:
        index(filepath.with_suffix('.vcf.bgz'))


def df2tsv(df, file_, header=True, separator='\t'):
    if isinstance(df, pl.DataFrame):
        df.write_csv(file_, include_header=header, separator=separator)
    elif isinstance(df, pd.DataFrame):
        df.to_csv(
            file_,
            header=header,
            index=index,
            sep=sep,
            na_rep=na_rep,
        )




def execute(cmd, debug=False, pipe=False):

    if debug:
        print(cmd)

    bag = []

    with Popen(cmd, shell=True, text=True, stdout=PIPE) as proc:
        for line in proc.stdout:
            if pipe:
                bag.append(line.strip())

        proc.wait()

        if proc.returncode:
            raise Exception(cmd)

    return bag


def tsv2df(input_file: Path,
           comment='#',
           header=0,
           sep='\t',
           dtype='str',
           keep_default_na: bool = True) -> pd.DataFrame:
    return pd.read_csv(
        input_file,
        comment=comment,
        header=header,
        sep=sep,
        dtype=dtype,
        keep_default_na=True,
    )



def chroms():
    x = ['chr' + str(x) for x in range(1,23)]
    x.extend(['chrX', 'chrY', 'chrM'])

    return sorted(x)


def revcom(seq):
    return seq.translate(COMPLEMENT_BASES)[::-1]


def save(obj, file_):
    with gzip.open(file_, 'wb') as fh:
        pickle.dump(obj, fh)


def load(file_):
    with gzip.open(file_, 'rb') as fh:
        return pickle.load(fh)


