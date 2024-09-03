#!/usr/bin/env python
import pgzip
import shutil
from pathlib import Path
from .genome import Genome
from .variant import Variant, sync
import gzip
from ncls import NCLS
from pathos.multiprocessing import ProcessPool
from .utils import chroms

from .vcf import Vcf

DB_DIRNAME='dbsnp_db'
BATCH_SIZE=10000000

PGZIP_BLOCK_SIZE= 16 * 1024 * 1024  # 16MB
BUFFER_SIZE = 1024 * 1024

DETAILS_HEADER = f'rsid\tchrom\tpos\tref\talt\tstart\tend'
SNVS_FILENAME = 'snvs.tsv.gz'
INTERVALS_FILENAME = 'intervals.bed.gz'
DETAILS_FILENAME = 'details.tsv.gz'
INDEX_FILENAME = 'index.tsv.gz'

CHROM_MAP_TEMPLATE = {
    'NC_000001': 'chr1',
    'NC_000002': 'chr2',
    'NC_000003': 'chr3',
    'NC_000004': 'chr4',
    'NC_000005': 'chr5',
    'NC_000006': 'chr6',
    'NC_000007': 'chr7',
    'NC_000008': 'chr8',
    'NC_000009': 'chr9',
    'NC_000010': 'chr10',
    'NC_000011': 'chr11',
    'NC_000012': 'chr12',
    'NC_000013': 'chr13',
    'NC_000014': 'chr14',
    'NC_000015': 'chr15',
    'NC_000016': 'chr16',
    'NC_000017': 'chr17',
    'NC_000018': 'chr18',
    'NC_000019': 'chr19',
    'NC_000020': 'chr20',
    'NC_000021': 'chr21',
    'NC_000022': 'chr22',
    'NC_000023': 'chrX',
    'NC_000024': 'chrY',
    'NC_012920': 'chrM',
}


class DbSnp():

    def __init__(self, chrom, db_dir, chromosome):

        self._chrom = chrom
        self.vcf_fh = gzip.open(vcf_file, 'rt')

        # self.vcf_file = vcf_file
        self.chromosome = chromosome
        self.index = _load_index(db_dir / chrom / INDEX_FILENAME)
        self.intervals = load_intervals(db_dir / chrom/ INTERVALS_FILENAME)

    @property
    def chrom(self):
        return self._chrom


    def query(self, variant, start, end):

        assert self._chrom == variant.chrom

        chrom = self._chrom

        rsids = list()
        results = self.intervals.find_overlap(start -1, end)

        for result in results:
            rsids.append(result[-1])


        results = list()
        details = list()

        for rsid in rsids:
            offset = self.index[rsid]
            self.vcf_fh.seek(offset)
            record = self.vcf_fh.readline().strip().split('\t')
            candidate = Variant(
                chrom = chrom,
                pos = int(record[1]), 
                id_ = record[2],
                ref = record[3],
                alt = record[4]
            )

            vx, vy = sync(variant, candidate, self.chromosome)
            if set(vx.alts) & set(vy.alts):
                results.append({
                    'chrom': chrom,
                    'pos': candidate.pos,
                    'ref': candidate.ref,
                    'alt': candidate.alt,
                    'rsid': vy.id,
                })

                details.append({
                    'chrom': chrom,
                    'pos': variant.pos,
                    'ref': variant.ref,
                    'alt': variant.alt,
                    'pos_sync': vx.pos,
                    'ref_sync': vx.ref,
                    'alt_sync': vx.alt,
                    'alt_dbsnp': candidate.alt,
                    'alt_dbsnp_sync': vy.alt,
                    'rsid': vy.id,
                })

        return {'result': results, 'details': details}



def _load_index(file_):
    index = dict()
    with gzip.open(file_, 'rt') as fh:
        for line in fh:
            record = line.strip().split('\t')
            rsid = int(record[0])
            offset = int(record[1])
            index[rsid] = offset
    return index


def _create_index(output_dir, n_threads):

    def jobs(output_dir):
        for chrom in chroms():
            chrom_dir = output_dir / chrom
            yield chrom_dir

    def process(job):
        chrom_dir = job
        snv_file = chrom_dir / SNVS_FILENAME
        index_file = chrom_dir / INDEX_FILENAME

        buffer_ = list()
        with gzip.open(snv_file,'rt') as fh, gzip.open(index_file, 'wt') as ifh:

            offset = fh.tell()
            line = fh.readline()
            n_records = 0


            while line:
                if not line.startswith('#'):
                    k = line.split('\t',3)[2][2:]
                    buffer_.append(f'{k}\t{offset}')
                    n_records += 1


                if n_records > BUFFER_SIZE:
                    data = '\n'.join(buffer_)
                    ifh.write(f'{data}\n')
                    n_records = 0
                    buffer_ = list()

                offset = fh.tell()
                line = fh.readline()
                

            if buffer_:
                data = '\n'.join(buffer_)
                ifh.write(f'{data}\n')
        return chrom_dir

    with ProcessPool(n_threads) as pool:
        for result in pool.uimap(process, jobs(output_dir)):
            print(result)


def _create_intervals(
        genome_file: Path,
        output_dir: Path,
        n_threads: int,
    ):

    def jobs(output_dir, genome_file):
        genome = Genome(genome_file)
        for chrom in chroms():
            chrom_dir = output_dir / chrom

            yield {
                'chrom_dir': chrom_dir,
                'chromosome': genome.chromosome(chrom),
            }

    def process(job):
        chrom_dir = job['chrom_dir']
        chromosome = job['chromosome']

        snvs_file = chrom_dir / SNVS_FILENAME
        intervals_file = chrom_dir / INTERVALS_FILENAME
        details_file = chrom_dir / DETAILS_FILENAME

        bed_buffer = list()
        csv_buffer = list()
        n_records = 0

        with gzip.open(snvs_file, 'rt') as fh, \
                gzip.open(intervals_file, 'wt') as ifh, \
                gzip.open(details_file, 'wt') as dfh:

            for line in fh:
                record = line.strip().split('\t')
                chrom = record[0]
                pos = int(record[1])
                rsid = record[2][2:]
                ref = record[3]
                alt = record[4]

                if alt == '.':
                    continue
                result = process_snv(chrom, pos, rsid, ref, alt, chromosome)
                start = result['start']
                end = result['end']
                bed_record = '\t'.join([str(x) for x in [chrom, start -1, end, rsid]])
                csv_record = '\t'.join([str(x) for x in [rsid, chrom, pos, ref, alt, start, end]])

                bed_buffer.append(bed_record)
                csv_buffer.append(csv_record)
                n_records += 1

                if n_records > BUFFER_SIZE:
                    bed_data = '\n'.join(bed_buffer)
                    csv_data = '\n'.join(csv_buffer)
                    ifh.write(f'{bed_data}\n')
                    dfh.write(f'{csv_data}\n')
                    n_records = 0
                    bed_buffer = list()
                    csv_buffer = list()

            if n_records:
                bed_data = '\n'.join(bed_buffer)
                csv_data = '\n'.join(csv_buffer)
                ifh.write(f'{bed_data}\n')
                dfh.write(f'{csv_data}\n')

        return chrom_dir.name

    def process_snv(chrom, pos, rsid, ref, alt, chromosome):
        start = None
        end = None

        for a in alt.split(','):
            v = Variant(chrom = chrom, pos = pos, ref = ref, alt = a)
            normalized = v.normalize(chromosome)
            denormalized = v.denormalize(chromosome)

            current_start = normalized.region.start
            current_end = denormalized.region.end

            if not start:
                start = current_start
                end = current_end
                continue

            if current_start < start:
                start = current_start
            if current_end > end:
                end = current_end

        return {'rsid': rsid, 'start': start, 'end': end}

    with ProcessPool(n_threads) as pool:
        for result in pool.uimap(process, jobs(output_dir, genome_file)):
            print(result)


def create_db(
    dbsnp_vcf_file: Path,
    genome_file: Path,
    output_dir: Path,
    n_threads: int,
):

    output_dir.mkdir(exist_ok = True)

    # # 64 threads, 25 min
    # _chop(dbsnp_vcf_file, output_dir, n_threads)
    #
    # # 64 threads, 60 min
    # _create_intervals(genome_file, output_dir, n_threads)

    _create_index(output_dir, n_threads)


    # intervals_file = output_dir / (dbsnp_vcf_file.name.split('.')[0] + '.bed.gz')
    # details_file = output_dir / (dbsnp_vcf_file.name.split('.')[0] + '.csv.gz')
    #
    # # 64 threads 20 min
    # with pgzip.open(intervals_file, 'wt', thread = n_threads, blocksize = PGZIP_BLOCK_SIZE) as ofh:
    #     for chrom in chroms():
    #         chrom_dir = output_dir / chrom
    #
    #         print(chrom_dir)
    #
    #         with pgzip.open(
    #                 (chrom_dir / INTERVALS_FILENAME), 'rt', 
    #                 thread = n_threads, blocksize = PGZIP_BLOCK_SIZE) as fh:
    #             for line in fh:
    #                 ofh.write(line)
    #
    # with pgzip.open(details_file, 'wt', thread = n_threads, blocksize = PGZIP_BLOCK_SIZE) as ofh:
    #
    #     ofh.write(f'{DETAILS_HEADER}\n')
    #
    #     for chrom in chroms():
    #         chrom_dir = output_dir / chrom
    #
    #         print(chrom_dir)
    #
    #         with pgzip.open(
    #                 (chrom_dir / DETAILS_FILENAME), 'rt', 
    #                 thread = n_threads, blocksize = PGZIP_BLOCK_SIZE) as fh:
    #             next(fh)
    #             for line in fh:
    #                 ofh.write(line)


def load_intervals(file_):

    bag = dict()

    with gzip.open(file_, 'rt') as fh:
        for line in fh:
            record = line.strip().split('\t')

            chrom = record[0]
            start = int(record[1])
            end = int(record[2])
            rsid = int(record[3])

            bag = {'start': list(), 'end': list(), 'rsid': list()}

            bag['start'].append(start)
            bag['end'].append(end)
            bag['rsid'].append(rsid)

    intervals = NCLS( bag['start'], bag['end'], bag['rsid'])

    return intervals



def normalize(
    dbsnp_vcf_file: Path,
    output_dir: Path,
    genome_file: Path,
    genome_index_file: Path,
    n_threads: int,
):
    if output_dir.exists():
        shutil.rmtree(output_dir)

    tmp_dir = output_dir / 'tmp'
    tmp_dir.mkdir(parents=True)

    vcf = Vcf(dbsnp_vcf_file, tmp_dir, n_threads)

    vcf = vcf.bgzip().index()
    contigs = vcf.contigs

    chrom_map_file = tmp_dir / 'chrommap.tsv'

    export_chrom_map(contigs, CHROM_MAP_TEMPLATE, chrom_map_file)

    vcf.rename_chroms(chrom_map_file)\
        .include_chroms(CHROM_MAP_TEMPLATE.values())\
        .fix_header(genome_index_file)\
        .drop_info() \
        .normalize(genome_file) \
        .move_to(output_dir / 'dbsnp.vcf.bgz')


def export_chrom_map(contigs: set, template: dict, output_file):

    bag = dict()
    for contig in contigs:
        k = contig.split('.')[0]

        if k in template:
            bag[contig] = template[k]

    with output_file.open('wt') as fh:

        for k, v in bag.items():
            fh.write(f'{k}\t{v}\n')


def _chop(dbsnp_vcf_file, tmp_dir, n_threads, blocksize = PGZIP_BLOCK_SIZE):

    class Writer():
        def __init__(self, tmp_dir, chrom, n_threads, blocksize):
            self.tmp_dir = tmp_dir / chrom
            self.tmp_dir.mkdir(exist_ok = True)
            self.output_file = self.tmp_dir / SNVS_FILENAME
            self.fh = pgzip.open(self.output_file,  'wt', thread = n_threads, blocksize = blocksize)
            self._chrom = chrom
            self._buffer = list()

        @property
        def chrom(self):
            return self._chrom

        def write(self, record):
            self._buffer.append(record)
            if len(self._buffer) > BUFFER_SIZE:
                chunk = '\n'.join(self._buffer)
                self.fh.write(f'{chunk}\n')
                self._buffer = list()

        def close(self):
            if self._buffer:
                chunk = '\n'.join(self._buffer)
                self.fh.write(f'{chunk}\n')
            try:
                self.fh.close()
            except Exception as e:
                raise e

    writer = None

    with gzip.open(dbsnp_vcf_file, 'rt') as fh:
        for line in fh:
            if line.startswith('##'): 
                continue
            break

        bag = list()
        for line in fh:
            record = line.strip().split('\t')[0:5]
            chrom = record[0]
            if not writer:
                writer = Writer(tmp_dir, chrom, n_threads, blocksize)

            if writer.chrom != chrom:
                writer.close()
                writer = Writer(tmp_dir, chrom, n_threads, blocksize)

            writer.write('\t'.join(record))
    if writer:
        writer.close()

