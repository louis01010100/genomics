from pathos.multiprocessing import ProcessPool
from .vcf import Vcf, concat
import polars as pl
from pathlib import Path

def process(
        sample_file: Path,
        target_population:str,
        workspace_dir:Path,
        vcf_dir:Path,
        genome_file:Path,
        genome_index_file:Path,
        n_threads:int,
    ):

    workspace_dir.mkdir(parents=True, exist_ok=True)
    tmp_dir = workspace_dir / 'tmp'

    target_samples = generate_sample_list(
        sample_file,
        target_population,
        workspace_dir,
    )

    vcf_files = [Vcf(x, tmp_dir) for x in vcf_dir.glob('*.vcf.gz')]

    bag = list()
    with ProcessPool(n_threads) as pool:
        for job in pool.uimap(
                standardize,
                jobs(vcf_files, target_samples, genome_file,
                         genome_index_file)):
            job.move_to(workspace_dir)
            print(f'DONE:\t{job.filepath}')

def jobs(vcf_files, target_samples, genome_file, genome_index_file):
    for vcf_file in vcf_files:
        yield {
            'vcf_file': vcf_file,
            'target_samples': target_samples,
            'genome_file': genome_file,
            'genome_index_file': genome_index_file,
        }


def standardize(args):

    vcf = args['vcf_file']
    target_samples = args['target_samples']
    genome_file = args['genome_file']
    genome_index_file = args['genome_index_file']

    return vcf.bgzip() \
            .subset_samples(target_samples) \
            .fix_header(genome_index_file) \
            .standardize(normalize = True)

def create_sample2gender_file(sample_info_file, sample2gender_file):

    samples = pl.read_excel(sample_info_file)

    samples.columns = [x.lower() for x in samples.columns]

    samples.select(['sample', 'gender']) \
            .write_csv(sample2gender_file, include_header = True, separator = '\t')

def generate_sample_list(
    sample_file,
    target_population,
    workspace_dir,
):
    target_sample_file = workspace_dir / 'samples_file.tsv'

    samples = pl.read_csv(
        sample_file,
        has_header = True,
        separator = '\t',
        infer_schema_length = None,
    )

    if target_population:
        samples = samples.filter(pl.col('Population') == target_population)

    samples.write_csv(target_sample_file, include_header = True, separator = '\t')


    return set(samples['Sample'])


