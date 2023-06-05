import gzip
from pathlib import Path
from subprocess import PIPE, STDOUT, Popen
from typing import TextIO, Tuple, Union


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


def is_gzipped(filepath):
    with open(filepath, 'rb') as fd:
        magic_number = fd.read(2)
        if magic_number == b'\x1f\x8b':
            return True
        return False


def vcf2dict(*vcf_files) -> dict:
    result = None

    for vcf_file in vcf_files:
        if is_gzipped(vcf_file):
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


def execute(cmd):
    with Popen(cmd, shell=True, text=True, stdout=PIPE,
               stderr=STDOUT) as proc:
        for line in proc.stdout:
            print(line)

        proc.wait()

        if proc.returncode:
            raise Exception(cmd)
