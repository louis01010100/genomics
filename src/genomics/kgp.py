import gzip
import shutil
from pathlib import Path
from subprocess import PIPE, STDOUT, Popen
import logging

import numpy as np
import polars as pl
from pathos.multiprocessing import ProcessPool

from .utils import load, save
from .variant import Variant
from .vcf import Vcf, concat, fetch_variants, filter_variants, list_samples, merge

MIN_READ_DEPTH = 10

TMP_DEPTH_FILENAME = 'depth.gz'

KGP_TRUTH_VCF_FILENAME = 'kgp_truth.vcf.bgz'
VCF_DIRNAME = 'vcf'
CRAM_DIRNAME = 'cram'


def export_snv_truth(
    vcf_dir,
    samples_file: Path,
    coordinates_vcf_file: Path,
    cram_dir: Path,
    genome_file: Path,
    output_dir: Path,
    n_threads: int = 1,
    min_read_depth: int = 2,
    n_cram_samples: int = 10,
):

    # shutil.rmtree(output_dir, ignore_errors=True)

    output_dir.mkdir(exist_ok=True)

    kgp_vcf_tmp_dir = output_dir / VCF_DIRNAME
    kgp_cram_tmp_dir = output_dir / CRAM_DIRNAME

    kgp_vcf_tmp_dir.mkdir(exist_ok=True)
    kgp_cram_tmp_dir.mkdir(exist_ok=True)

    coordinates = mung_coordinates(coordinates_vcf_file, output_dir)

    depths = extract_depth(
        cram_dir,
        genome_file,
        coordinates,
        kgp_cram_tmp_dir,
        n_threads,
        n_samples=n_cram_samples,
    )

    save(depths, output_dir / TMP_DEPTH_FILENAME)

    depths = load(output_dir / TMP_DEPTH_FILENAME)

    print('depths loaded')

    samples = set(
        pl.read_csv(
            samples_file,
            has_header=True,
            separator='\t',
        )['sample_id'])

    # vcf_file = subset_samples(
    #     vcf_dir,
    #     kgp_vcf_tmp_dir,
    #     samples,
    #     n_threads,
    # )

    vcf_file = kgp_vcf_tmp_dir / 'kgp-samples.vcf.bgz'

    result_vcf = subset_snvs(
        vcf_file=vcf_file,
        coordinates=coordinates,
        depths=depths,
        output_dir=kgp_vcf_tmp_dir,
        n_threads=n_threads,
    )

    kgp_truth_vcf = result_vcf.move_to(output_dir / 'kgp_truth.vcf.bgz')


def mung_coordinates(coordinates_vcf_file, output_dir) -> pl.DataFrame:

    vcf_file1 = output_dir / 'file1.vcf'
    vcf_file2 = output_dir / 'file2.vcf'

    with gzip.open(coordinates_vcf_file, 'rt') as fh, vcf_file1.open(
            'wt') as fh1, vcf_file2.open('wt') as fh2:
        for line in fh:
            if line.startswith('#'):
                fh1.write(line)
                fh2.write(line)
                continue
            fh1.write(line)
            fh2.write(line)
            break

        coordinates = pl.read_csv(
            coordinates_vcf_file,
            comment_char='#',
            has_header=False,
            new_columns=[
                'chrom', 'pos', 'id', 'ref', 'alt', 'qual', 'filter', 'info'
            ],
            separator='\t',
        )

        for coordinate, data in coordinates.groupby(['chrom', 'pos']):
            data2 = data.drop('id').unique().to_dicts()

            record = data2[0]

            fh1.write(
                f"{record['chrom']}\t{record['pos']}\t\t{record['ref']}\t{record['alt']}\t\t\t\n"
            )

            for record in data2[1:]:
                fh2.write(
                    f"{record['chrom']}\t{record['pos']}\t\t{record['ref']}\t{record['alt']}\t\t\t\n"
                )

    output_vcf_file = output_dir / 'coordinates.vcf.bgz'

    merge([vcf_file2, vcf_file2], output_dir, output_vcf_file, flag='all')

    return Vcf(output_vcf_file, output_dir).to_df(site_only=True)


def subset_snvs(
    vcf_file: Path,
    coordinates: pl.DataFrame,
    depths: dict,
    output_dir: Path,
    n_threads: int,
):

    def jobs(coordinates, vcf_file):
        for record in coordinates.to_dicts():
            yield {
                'chrom': record['chrom'],
                'pos': int(record['pos']),
                'ref': record['ref'],
                'alt': record['alt'],
                'vcf_file': vcf_file,
            }

    def process(job):
        chrom = job['chrom']
        pos = int(job['pos'])
        ref = job['ref']
        alt = job['alt']
        vcf_file = job['vcf_file']

        refsnv = Variant(chrom, pos, ref, alt)

        snvs = fetch_variants(refsnv.chrom, refsnv.pos, vcf_file)

        snvs = filter_variants(refsnv, snvs)

        return {
            'snvs': snvs,
            'chrom': chrom,
            'pos': pos,
            'ref': ref,
            'alt': alt,
        }

    n_samples = len(list_samples(vcf_file))

    output_file = output_dir / vcf_file.name.replace('.vcf.bgz', '-snv.vcf')

    with gzip.open(vcf_file, 'rt') as ifh, output_file.open('wt') as ofh:
        for line in ifh:
            if line.startswith('##'):
                ofh.write(line)
                continue
            ofh.write(line)
            break

    n_total = len(coordinates)
    n_done = 0
    with output_file.open('at') as fh:

        with ProcessPool(n_threads) as pool:
            for result in pool.uimap(process, jobs(coordinates, vcf_file)):
                n_done += 1
                print(f'{n_done / n_total}', end='\r')
                snvs = result['snvs']

                if len(snvs) == 1:
                    snv = snvs[0]
                    fh.write(f'{snv.data}\n')
                    continue

                if len(snvs) > 1:
                    logging.warning(f'{snvs}')
                    snv = snvs[0]
                    fh.write(f'{snv.data}\n')
                    continue

                chrom = result['chrom']
                pos = str(result['pos'])
                ref = result['ref']
                alt = result['alt']

                key = (chrom, pos)

                if key not in depths:
                    gt = '\t'.join(['./.' for i in range(0, n_samples)])
                else:
                    depth = depths[key]

                    if depth > MIN_READ_DEPTH:
                        gt = '\t'.join(['0/0' for i in range(0, n_samples)])
                    else:
                        print(key)
                        gt = '\t'.join(['./.' for i in range(0, n_samples)])

                fh.write(
                    f'{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\t.\t.\tGT\t{gt}\n')

    # output_file = 'workspace/vcf/kgp-samples-snv.vcf'
    return Vcf(output_file, output_dir, n_threads).bgzip().sort().index()


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

    with ProcessPool(n_threads) as pool:
        for vcf_file in pool.uimap(
                process,
                jobs(vcf_dir, output_dir, samples),
        ):
            bag.append(vcf_file)

    # bag = [x for x in output_dir.glob('*/*trim_alt.vcf.bgz')]

    output_file = output_dir / 'kgp-samples.vcf.bgz'

    result_vcf = concat(
        vcf_files=bag,
        output_file=output_file,
        tmp_dir=output_dir,
        n_threads=n_threads,
    )

    return output_file


def extract_depth(
    cram_dir: Path,
    genome_file: Path,
    coordinates: pl.DataFrame,
    output_dir: Path,
    n_threads: int,
    n_samples: int,
):

    def jobs(cram_dir, genome_file, coordinates_file):

        n = 0

        for cram_file in cram_dir.glob('*.cram'):
            yield {
                'cram_file': cram_file,
                'coordinates_file': coordinates_file,
                'genome_file': genome_file,
            }
            n += 1
            if n >= n_samples:
                break

    def parse_cram(args):
        cram_file = args['cram_file']
        genome_file = args['genome_file']
        coordinates_file = args['coordinates_file']

        cmd = (''
               f'samtools mpileup'
               f'    -l {coordinates_file}'
               f'    -f {genome_file}'
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
                jobs(cram_dir, genome_file, coordinates_file),
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
    for key, df in data.groupby(['chrom', 'pos']):
        result[key] = np.median(df['depth'])

    return result
