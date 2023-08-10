import shutil
from pathlib import Path
from subprocess import PIPE, STDOUT, Popen

import numpy as np
import polars as pl
from pathos.multiprocessing import ProcessPool

from .vcf import Vcf, concat


def process(
    vcf_dir,
    samples_file: Path,
    coordinate_vcf_file: Path,
    cram_dir: Path,
    ref_file: Path,
    output_dir: Path,
    n_threads: int = 1,
    min_read_depth: int = 2,
    n_cram_samples: int = 10,
):

    # if output_dir.exists():
    #     shutil.rmtree(output_dir)
    output_dir.mkdir(exist_ok=True)

    kgp_vcf_tmp_dir = output_dir / 'vcf'
    kgp_cram_tmp_dir = output_dir / 'cram'

    kgp_vcf_tmp_dir.mkdir(exist_ok=True)
    kgp_cram_tmp_dir.mkdir(exist_ok=True)

    coordinates = pl.from_pandas(
        Vcf(coordinate_vcf_file, kgp_vcf_tmp_dir,
            n_threads).to_df(site_only=True))

    samples = set(
        pl.read_csv(
            samples_file,
            has_header=True,
            separator='\t',
        )['sample_id'])

    variants = extract_variants(
        vcf_dir=vcf_dir,
        output_dir=kgp_vcf_tmp_dir,
        samples=samples,
        coordinates=coordinates,
        n_threads=n_threads,
    )

    tmp = result.join(
        coordinates.with_columns(pl.lit(True).alias('coord')),
        on=['chrom', 'pos'],
        how='left',
    )

    coordinate_missing = tmp.filter(pl.col('coord') == False).select(
        pl.col(['chrom', 'pos']))

    result = extract_depth(
        cram_dir,
        ref_file,
        coordinate_missing,
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

    merged_vcf_file = subset_samples(vcf_dir, output_dir, samples, n_threads)

    merged_vcf = Vcf(merged_vcf_file, output_dir)
    merged_vcf.index()

    variants = subset_variants(
        merged_vcf,
        output_dir,
        coordinates,
        n_threads,
    )

    variants = variants.with_columns(
        pl.col('pos').cast(pl.Int64).alias('pos'))

    variants_file = output_dir / 'variants.tsv'

    variants.write_csv(variants_file, has_header=True, separator='\t')

    variants = pl.read_csv(variants_file, has_header=True, separator='\t')

    return variants


def subset_variants(merged_vcf, output_dir, coordinates, n_threads):

    def jobs(coordinates, merged_vcf):
        for coord in coordinates.to_dicts():
            yield {
                'chrom': coord['chrom'],
                'pos': coord['pos'],
                'id': coord['id'],
                'ref': coord['ref'],
                'alt': coord['alt'],
                'vcf': merged_vcf,
            }

    def process(job):
        chrom = job['chrom']
        pos = job['pos']
        id_ = job['id']
        ref = job['ref']
        alt = job['alt']

        print(chrom, pos, id_, ref, alt)

        merged_vcf = job['vcf']

        records = merged_vcf.fetch_variant(chrom, pos)

        return records

    bag = []
    with ProcessPool(n_threads) as pool:
        for records in pool.uimap(process, jobs(coordinates, merged_vcf)):
            bag.extend(records)

    data = pl.from_dicts(bag)

    return data


def subset_samples(vcf_dir, output_dir, samples, n_threads):

    def jobs(vcf_dir, outupt_dir, samples):
        for vcf_file in vcf_dir.glob('*.vcf.bgz'):
            yield {
                'vcf_file': vcf_file,
                'output_dir': output_dir,
                'samples': samples,
            }

    def process(job):
        vcf_file = job['vcf_file']
        output_dir = job['output_dir']
        samples = job['samples']
        return Vcf( vcf_file, output_dir) \
                .subset_samples(samples) \
                .trim_alts() \
                .filepath

    bag = []

    # with ProcessPool(n_threads) as pool:
    #     for vcf_file in pool.uimap(
    #             process,
    #             jobs(vcf_dir, output_dir, samples),
    #     ):
    #         bag.append(vcf_file)

    output_file = output_dir / 'kgp.vcf.bgz'
    result_vcf = concat(
        vcf_files=bag,
        output_file=output_file,
        tmp_dir=output_dir,
        n_threads=n_threads,
    )

    return output_file


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

        for cram_file in cram_dir.glob('*.cram'):
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
        sample_name = cram_file.name.split('.')[0]
        return sample_name, pl.from_dicts(bag)

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

    bag = []
    for f in output_dir.glob('*.tsv'):
        if f.name == 'coordinates.tsv':
            continue
        d = pl.read_csv(f, has_header=True, separator='\t')
        bag.append(d)

    data = pl.concat(bag)

    bag = []
    for keys, df in data.groupby(['chrom', 'pos', 'ref']):
        median = np.median(df['depth'])
        bag.append({
            'chrom': keys[0],
            'pos': keys[1],
            'ref': keys[2],
            'depth': median
        })

    return pl.from_dicts(bag)
