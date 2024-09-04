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

BATCH_SIZE=10000000

PGZIP_BLOCK_SIZE= 16 * 1024 * 1024  # 16MB
BUFFER_SIZE = 1024 * 1024

DETAILS_HEADER = f'rsid\tchrom\tpos\tref\talt\tstart\tend'
SNVS_FILENAME = 'snvs.tsv.gz'
INTERVALS_FILENAME = 'intervals.tsv.gz'
DETAILS_FILENAME = 'details.tsv.gz'
RSID2OFFSET_FILENAME = 'rsid2offset.tsv.gz'
INDEX_FILENAME = 'idx.tsv.gz'
DB_FILENAME = 'db.tsv.gz'

DB_COLUMNS = ['chrom', 'pos', 'rsid', 'ref', 'alt', 'start_bed', 'end_bed', 'idx']

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
        self.chromosome = chromosome
        self._db = _load_db(db_dir / chrom / DB_FILENAME)
        self._idxs = _load_idxs(db_dir / chrom/ INDEX_FILENAME)
        self._db_fh = gzip.open(db_dir / chrom / DB_FILENAME, 'rt')
        self._col2idx = {c : i for i, c in enumerate(DB_COLUMNS)}

    @property
    def chrom(self):
        return self._chrom


    def query(self, variant, start, end):

        assert self._chrom == variant.chrom

        chrom = self._chrom

        intervals = self._db.find_overlap(start -1, end)
        idxs = [interval[-1] for interval in intervals]

        for idx in idxs:
            self._db_fh.seek(self._idxs[idx])
            snv = self._db_fh.readline().strip().split('\t')

            candidate = Variant(
                chrom = chrom,
                pos = int(snv[self._col2idx['pos']]), 
                id_ = snv[self._col2idx['rsid']], 
                ref = snv[self._col2idx['ref']], 
                alt = snv[self._col2idx['alt']], 
            )

            vx, vy = sync(variant, candidate, self.chromosome)
            if set(vx.alts) & set(vy.alts):

                results.append({
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

        return results

# 64 threads, 1hr32min
def create_db(
    dbsnp_vcf_file: Path,
    genome_file: Path,
    output_dir: Path,
    n_threads: int,
):

    output_dir.mkdir(exist_ok = True)

    # # 64 threads, 25 min
    _chop(dbsnp_vcf_file, output_dir, n_threads)

    # 64 threads, 60 min
    _create_db(genome_file, output_dir, n_threads)

    _create_idxs(output_dir, n_threads)


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

    _export_chrom_map(contigs, CHROM_MAP_TEMPLATE, chrom_map_file)

    vcf.rename_chroms(chrom_map_file)\
        .include_chroms(CHROM_MAP_TEMPLATE.values())\
        .fix_header(genome_index_file)\
        .drop_info() \
        .normalize(genome_file) \
        .move_to(output_dir / 'dbsnp.vcf.bgz')

def _create_idxs(output_dir, n_threads):

    def jobs(output_dir):
        for chrom in chroms():
            chrom_dir = output_dir / chrom
            yield chrom_dir

    def process(job):
        chrom_dir = job
        db_file = chrom_dir / DB_FILENAME
        idx_file = chrom_dir / INDEX_FILENAME

        buffer_ = list()
        with gzip.open(db_file,'rt') as fh, gzip.open(idx_file, 'wt') as ifh:

            header = fh.readline()

            offset = fh.tell()
            line = fh.readline()
            n_records = 0

            ifh.write('idx\toffset\n')


            while line:
                record = line.strip().split('\t')
                idx = record[-1]
                buffer_.append(f'{idx}\t{offset}')
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


def _create_db(
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
        db_file = chrom_dir / DB_FILENAME

        buffer_ = list()
        n_records = 0

        idx = 0

        with gzip.open(snvs_file, 'rt') as fh, \
                gzip.open(db_file, 'wt') as ifh:

            header = '\t'.join(DB_COLUMNS)

            ifh.write(f'{header}\n')

            for line in fh:
                record = line.strip().split('\t')
                chrom = record[0]
                pos = int(record[1])
                rsid = record[2]
                ref = record[3]
                alt = record[4]

                if alt == '.':
                    continue
                
                result = process_snv(chrom, pos, rsid, ref, alt, chromosome)
                start = result['start']
                end = result['end']
                record = '\t'.join([str(x) for x in [chrom, pos, rsid, ref, alt, start - 1, end, idx]])

                buffer_.append(record)
                n_records += 1
                idx += 1

                if n_records > BUFFER_SIZE:
                    bed_data = '\n'.join(buffer_)
                    ifh.write(f'{bed_data}\n')
                    n_records = 0
                    buffer_ = list()

            if n_records:
                bed_data = '\n'.join(buffer_)
                ifh.write(f'{bed_data}\n')

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



def _load_idxs(file_):
    bag = dict()
    with gzip.open(file_, 'rt') as fh:

        col2idx = {col: idx for idx, col in enumerate(next(fh).strip().split('\t'))}
        for line in fh:
            record = line.strip().split('\t')
            idx = int(record[col2idx['idx']])
            offset = int(record[col2idx['offset']])
            bag[idx] = offset
    return bag



def _load_db(file_):

    starts = list()
    ends = list()
    indices = list()

    with gzip.open(file_, 'rt') as fh:

        col2idx = {col: idx for idx, col in enumerate(next(fh).strip().split('\t'))}
        for line in fh:
            record = line.strip().split('\t')
            idx = int(record[col2idx['idx']])
            start_bed = int(record[col2idx['start_bed']])
            end_bed = int(record[col2idx['end_bed']])

            starts.append(start_bed)
            ends.append(end_bed)
            indices.append(idx)


    intervals = NCLS( starts, ends, indices)

    return intervals





def _export_chrom_map(contigs: set, template: dict, output_file):

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

