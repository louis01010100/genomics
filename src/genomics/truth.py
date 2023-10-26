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
from .coordinates import export_coordinates

MIN_READ_DEPTH = 10

TMP_DEPTH_FILENAME = 'depth.gz'

TRUTH_VCF_FILE = 'truth.vcf.bgz'
VCF_DIRNAME = 'vcf'
CRAM_DIRNAME = 'cram'
COORDINATES_FILENAME = 'coordinates.vcf.bgz'


def export_snv_truth(
    vcfs_file: Path,
    crams_file: Path,
    samples_file: Path,
    coordinates_vcf_file: Path,
    genome_file: Path,
    output_dir: Path,
    n_threads: int = 1,
    min_read_depth: int = 2,
    debug=True,
):

    if not debug:
        shutil.rmtree(output_dir, ignore_errors=True)

    output_dir.mkdir(exist_ok=True)

    vcf_tmp_dir = output_dir / VCF_DIRNAME
    cram_tmp_dir = output_dir / CRAM_DIRNAME

    vcf_tmp_dir.mkdir(exist_ok=True)
    cram_tmp_dir.mkdir(exist_ok=True)

    export_coordinates(
        coordinates_vcf_file,
        output_dir / COORDINATES_FILENAME,
    )

    coordinates = Vcf(output_dir / COORDINATES_FILENAME,
                      output_dir).to_df(site_only=True)

    depths = extract_depth(
        crams_file,
        genome_file,
        coordinates,
        cram_tmp_dir,
        n_threads,
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

    vcf_file = subset_samples(
        vcfs_file=vcfs_file,
        genome_file=genome_file,
        samples=samples,
        output_dir=vcf_tmp_dir,
        n_threads=n_threads,
    )

    # vcf_file = vcf_tmp_dir / 'samples.vcf.bgz'

    result_vcf = subset_snvs(
        vcf_file=vcf_file,
        coordinates=coordinates,
        depths=depths,
        output_dir=vcf_tmp_dir,
        n_threads=n_threads,
    )

    kgp_truth_vcf = result_vcf.move_to(output_dir / TRUTH_VCF_FILE)


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
                snv = result['snv']
                note = result['note']

                if note == 'done':
                    pass
                elif note == 'complex':
                    snv.format = 'GT'
                    snv.calls = '\t'.join(
                        ['./.' for i in range(0, n_samples)])
                elif note == 'missing':
                    snv.format = 'GT'
                    key = (snv.chrom, str(snv.pos))

                    if key not in depths:
                        snv.calls = '\t'.join(
                            ['./.' for i in range(0, n_samples)])
                    else:
                        depth = depths[key]

                        if depth > MIN_READ_DEPTH:
                            snv.calls = '\t'.join(
                                ['0/0' for i in range(0, n_samples)])
                        else:
                            print(key)
                            snv.calls = '\t'.join(
                                ['./.' for i in range(0, n_samples)])
                else:
                    raise Exception(note)

                for id_ in snv.id.split(','):
                    v = Variant(
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
                    )
                    fh.write(f'{v}\n')

    return Vcf(output_file, output_dir, n_threads).bgzip().drop_qual(
    ).drop_filter().drop_info().fill_tags().sort().index()


def subset_samples(vcfs_file, genome_file, samples, output_dir, n_threads):

    def jobs(vcf_files, outupt_dir, samples):
        for vcf_file in vcf_files:
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
                .normalize(genome_file) \
                .drop_info() \
                .drop_format(fields = ['AB','AD', 'DP', 'GQ', 'PGT', 'PID','PL'], ) \
                .trim_alts() \
                .exclude('ALT="."') \
                .index() \
                .filepath

    bag = []

    vcf_files = set([
        Path(x) for x in pl.read_csv(
            vcfs_file, has_header=True, separator='\t')['vcf_path']
    ])

    with ProcessPool(n_threads) as pool:
        for vcf_file in pool.uimap(
                process,
                jobs(vcf_files, output_dir, samples),
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


def extract_depth(
    crams_file: Path,
    genome_file: Path,
    coordinates: pl.DataFrame,
    output_dir: Path,
    n_threads: int,
):

    def jobs(cram_files, genome_file, coordinates_file):
        for cram_file in cram_files:
            yield {
                'cram_file': cram_file,
                'coordinates_file': coordinates_file,
                'genome_file': genome_file,
            }

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

    cram_files = list(
        set([
            Path(x) for x in pl.read_csv(
                crams_file, has_header=True, separator='\t')['cram_path']
        ]))

    with ProcessPool(n_threads) as pool:
        for future in pool.uimap(
                parse_cram,
                jobs(cram_files, genome_file, coordinates_file),
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
