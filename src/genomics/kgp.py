import shutil
from pathlib import Path

import polars as pl
from pathos.multiprocessing import ProcessPool

from .vcf import Vcf, concat


def process(
    vcf_dir,
    samples_file: Path,
    coordinates_file: Path,
    cram_dir: Path,
    ref_file: Path,
    output_dir: Path,
    n_threads: int = 1,
    min_read_depth: int = 2,
    n_cram_samples: int = 10,
):

    if output_dir.exists():
        shutil.rmtree(output_dir)
    output_dir.mkdir(exist_ok=True)

    kgp_vcf_tmp_dir = output_dir / 'vcf'
    kgp_cram_tmp_dir = output_dir / 'cram'

    kgp_vcf_tmp_dir.mkdir(exist_ok=True)
    kgp_cram_tmp_dir.mkdir(exist_ok=True)

    coordinates = pl.read_csv(
        coordinates_file,
        has_header=True,
        separator='\t',
    )

    samples = set(
        pl.read_csv(
            samples_file,
            has_header=True,
            separator='\t',
        )['sample_id'])

    result_vcf = extract_variants(
        vcf_dir=vcf_dir,
        output_dir=kgp_vcf_tmp_dir,
        samples=samples,
        coordinates=coordinates,
        n_threads=n_threads,
    )

    print(result_vcf.filepath)

    result = extract_depth(
        cram_dir,
        ref_file,
        coordinates,
        kgp_cram_tmp_dir,
        n_threads,
        n_samples=n_cram_samples,
    )

    print(result)


def extract_variants(
    vcf_dir: Path,
    output_dir: Path,
    samples: set,
    coordinates: pl.DataFrame,
    n_threads: int,
):

    def jobs(vcf_dir, outupt_dir, samples, coordinates):
        coordinates = coordinates.to_pandas()
        for vcf_file in vcf_dir.glob('*.vcf.bgz'):
            yield {
                'vcf_file': vcf_file,
                'output_dir': output_dir,
                'samples': samples,
                'coordinates': coordinates,
            }

    def process(job):
        vcf_file = job['vcf_file']
        output_dir = job['output_dir']
        samples = job['samples']
        coordinates = job['coordinates']
        return Vcf( vcf_file, output_dir) \
                .subset_variants(coordinates) \
                .trim_alts() \
                .subset_samples(samples) \
                .filepath

    bag = []
    with ProcessPool(n_threads) as pool:
        for vcf_file in pool.uimap(
                process, jobs(
                    vcf_dir,
                    output_dir,
                    samples,
                    coordinates,
                )):
            bag.append(vcf_file)

    result_vcf = concat(bag, output_dir, 'kgp.vcf.gz', n_threads)

    return result_vcf


def extract_depth(
    cram_dir: Path,
    ref_file: Path,
    coordinates: pl.DataFrame,
    output_dir: Path,
    n_threads: int,
    n_samples: int,
):

    def jobs(cram_dir, ref_file, coordinates_file):

        n = 0

        for cram_file in cram_dir.glob('*'):
            yield {
                'cram_file': cram_file,
                'coordinates_file': coordinates_file,
                'ref_file': ref_file,
            }
            n += 1
            if n >= n_samples:
                break

    def parse_cram(args):
        cram_file = args['cram_file']
        ref_file = args['ref_file']
        coordinates_file = args['coordinates_file']

        cmd = (''
               f'samtools mpileup'
               f'    -l {coordinates_file}'
               f'    -f {ref_file}'
               f'    {cram_file}'
               '')

        with Popen(cmd, shell=True, text=True, stdout=PIPE) as proc:
            bag = []
            for line in proc.stdout:
                items = line.strip().split('\t')
                chrom = items[0]
                pos = items[1]
                ref = items[2]
                depth = items[3]

                bag.append({
                    'chrom': chrom,
                    'pos': pos,
                    'ref': ref,
                    'depth': depth
                })

            proc.wait()

            if proc.returncode:
                raise Exception(cmd)
        return cram_file.name, pl.from_dicts(bag)

    coordinates_file = output_dir / 'coordinates.tsv'

    coordinates.write_csv(coordinates_file, has_header=False, separator='\t')

    depth_bag = {}
    with ProcessPool(n_threads) as pool:
        for future in pool.uimap(
                parse_cram,
                jobs(cram_dir, ref_file, coordinates_file),
        ):
            sample_name = future[0]
            depth_df = future[1]

            output_file = output_dir / f'{sample_name}.tsv'
            depth_df.write_csv(output_file, has_header=True, separator='\t')

            depth_bag[sample_name] = output_file

    return depth_bag
