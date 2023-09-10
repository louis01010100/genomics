import gzip
import shutil
from pathlib import Path
from subprocess import PIPE, STDOUT, Popen

import numpy as np
import polars as pl
from pathos.multiprocessing import ProcessPool

from .utils import load, save
from .variant import Variant
from .vcf import Vcf, concat, fetch_variants, filter_variants, list_samples

MIN_READ_DEPTH = 10


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

    coordinates = Vcf(coordinate_vcf_file, kgp_vcf_tmp_dir,
                      n_threads).to_df(site_only=True)

    # depths = extract_depth(
    #     cram_dir,
    #     ref_file,
    #     coordinates,
    #     kgp_cram_tmp_dir,
    #     n_threads,
    #     n_samples=n_cram_samples,
    # )
    #
    # save(depths, output_dir / 'depths.gz')

    depths = load(output_dir / 'depths.gz')

    samples = set(
        pl.read_csv(
            samples_file,
            has_header=True,
            separator='\t',
        )['sample_id'])

    vcf_file = subset_samples(
        vcf_dir,
        kgp_vcf_tmp_dir,
        samples,
        n_threads,
    )

    subset_snvs(
        vcf_file=vcf_file,
        coordinates=coordinates,
        depths=depths,
        output_dir=kgp_vcf_tmp_dir,
        n_threads=n_threads,
    ).write_csv(
        output_dir / 'kgp_truth.tsv',
        has_header=True,
        separator='\t',
    )


def subset_snvs(
    vcf_file: Path,
    coordinates: pl.DataFrame,
    depths: dict,
    output_dir: Path,
    n_threads: int,
):

    n_samples = len(list_samples(vcf_file))

    output_file = output_dir / vcf_file.name.replace('.vcf.bgz', '-snv.vcf')

    with gzip.open(vcf_file, 'rt') as ifh, output_file.open('wt') as ofh:
        for line in ifh:
            if line.startswith('##'):
                ofh.write(line)
                continue
            ofh.write(line)
            break

    with output_file.open('at') as fh:
        for coord in coordinates.to_dicts():
            chrom = coord['chrom']
            pos = int(coord['pos'])
            ref = coord['ref']
            alt = coord['alt']

            refsnv = Variant(chrom, pos, ref, alt)

            snvs = fetch_variants(refsnv.chrom, refsnv.pos, vcf_file)

            snv = filter_variants(refsnv, snvs)

            if snv:
                fh.write(f'{snv.data}\n')
                continue
            key = (chrom, pos)

            if key not in depths:
                gt = '\t'.join(['./.' for i in range(0, n_samples)])
                continue

            depth = depths[key]

            if depth > MIN_READ_DEPTH:
                gt = '\t'.join(['0/0' for i in range(0, n_samples)])
            else:
                gt = '\t'.join(['./.' for i in range(0, n_samples)])

            fh.write(f'{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\t.\t.\tGT\t{gt}\n')

    return Vcf(output_file, output_dir, n_threads).bgzip().index().to_df()


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
        return Vcf(vcf_file, output_dir) \
                .subset_samples(samples) \
                .drop_info() \
                .drop_format(fields = ['AB','AD', 'DP', 'GQ', 'PGT', 'PID','PL'], ) \
                .trim_alts() \
                .index() \
                .filepath

    bag = []

    # with ProcessPool(n_threads) as pool:
    #     for vcf_file in pool.uimap(
    #             process,
    #             jobs(vcf_dir, output_dir, samples),
    #     ):
    #         bag.append(vcf_file)

    bag = [x for x in output_dir.glob('*/*trim_alt.vcf.bgz')]

    output_file = output_dir / 'kgp.vcf.bgz'

    # result_vcf = concat(
    #     vcf_files=bag,
    #     output_file=output_file,
    #     tmp_dir=output_dir,
    #     n_threads=n_threads,
    # )

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

    with ProcessPool(n_threads) as pool:
        for future in pool.uimap(
                parse_cram,
                jobs(cram_dir, ref_file, coordinates_file),
        ):
            sample_name = future[0]
            depth_df = future[1]

            output_file = output_dir / f'{sample_name}.tsv'
            depth_df.write_csv(output_file, has_header=True, separator='\t')

    bag = []
    for f in output_dir.glob('*.tsv'):
        if f.name == 'coordinates.tsv':
            continue
        d = pl.read_csv(
            f,
            has_header=True,
            separator='\t',
            infer_schema_length=0,
        ).with_columns(pl.col('depth').cast(pl.Float64))
        bag.append(d)

    data = pl.concat(bag)

    result = dict()
    for key, df in data.groupby(['chrom', 'pos', 'ref']):
        result[key] = np.median(df['depth'])

    return result
