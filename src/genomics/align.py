import gzip
import pickle
import zlib
from datetime import datetime
from pathlib import Path

import polars as pl
from icecream import ic
from pathos.multiprocessing import ProcessPool

COMPLEMENT_DNA = str.maketrans('ACGT', 'TGCA')

FRAGMENT_LENGTH = 100    # 18 min
FRAGMENT_LENGTH = 1000    #  5 min
FRAGMENT_LENGTH = 10000    #  3 min
FRAGMENT_LENGTH = 100000    #  3 min

TMP_DIR = Path('tmp')


def align(query_file: Path, genome_file: Path, output_file: Path,
          n_threads: int):

    TMP_DIR.mkdir(exist_ok=True)

    time_begin = datetime.now()

    queries = load_queries(query_file)

    print(f'queries_loaded\t{datetime.now()}')
    read_length = len(queries[0]['read1'])
    genome = load_genome(genome_file)
    print(f'genome loaded\t{datetime.now()}')
    genome_size = get_genome_size(genome)

    genome_index = index(
        genome=genome,
        genome_length=genome_size,
        fragment_length=FRAGMENT_LENGTH,
        read_length=read_length,
        n_threads=n_threads,
    )

    # print(f'genome_index done\t{datetime.now()}')
    # save(genome_index, TMP_DIR / 'genome_index.obj')
    #
    # genome_index = load(TMP_DIR / 'genome_index.obj')

    print(f'genome_index created\t{datetime.now()}')

    result = _align(queries, genome_index, n_threads)

    result.write_csv(output_file, has_header=True, separator='\t')
    print(f'query ended\t{datetime.now()}')


# def index_old(
#     genome: dict,
#     read_length: int,
#     genome_size: int,
#     batch_size: int,
#     n_threads: int,
# ):
#
#     def new_jobs(genome, read_length, batch_size):
#
#         bag = list()
#
#         for chrom, seq in genome.items():
#             for i in range(0, len(seq) - read_length + 1, 1):
#                 job = {
#                     'chrom': chrom,
#                     'seed': seq[i:i + read_length],
#                     'start': i,
#                     'read_length': read_length,
#                 }
#
#                 bag.append(job)
#
#                 if i % batch_size == 0:
#                     yield bag
#                     bag = list()
#         if len(bag):
#             yield bag
#
#     def process(jobs):
#         bag = dict()
#         n_jobs = len(jobs)
#         for job in jobs:
#             chrom = job['chrom']
#             seed = job['seed']
#             start = job['start']
#             read_length = job['read_length']
#
#             end = start + read_length
#
#             zeed = zlib.compress(seed.encode())
#
#             if zeed not in bag:
#                 bag[zeed] = list()
#             bag[zeed].append((chrom, '+', start, end))
#
#             seed_revcom = seed[::-1].translate(COMPLEMENT_DNA)
#
#             zeed_revcom = zlib.compress(seed_revcom.encode())
#
#             if zeed_revcom not in bag:
#                 bag[zeed_revcom] = list()
#             bag[zeed_revcom].append((chrom, '-', start, end))
#         return bag, n_jobs
#
#     bag = dict()
#
#     n_total = genome_size * read_length
#     n_done = 0
#     with ProcessPool(n_threads) as pool:
#         for result, n in pool.uimap(
#                 process,
#                 new_jobs(
#                     genome,
#                     read_length,
#                     batch_size,
#                 ),
#         ):
#             bag = merge(bag, result)
#             n_done += n
#
#             # print(f'{n_done / n_total}', end='\r', flush=True)
#     return bag


def index(
    genome: dict,
    read_length: int,
    fragment_length: int,
    genome_length: int,
    n_threads: int,
):

    def iter_seq(seq, fragment_length, read_length):
        i = 0

        seq_length = len(seq)

        while True:
            j = i + fragment_length

            if j >= seq_length:
                j = seq_length

            yield {
                'fragment': seq[i:j],
                'fragment_start': i,
                'fragment_end': j,
            }

            if j == seq_length:
                break

            i = j - read_length + 1

    def new_jobs(genome, fragment_length, read_length):

        bag = list()

        for chrom, seq in genome.items():

            for record in iter_seq(seq, fragment_length, read_length):
                yield {
                    'chrom': chrom,
                    'fragment': record['fragment'],
                    'fragment_start': record['fragment_start'],
                    'fragment_end': record['fragment_end'],
                    'read_length': read_length,
                }

    def process(job):
        chrom = job['chrom']
        fragment = job['fragment']
        fragment_start = job['fragment_start']
        fragment_end = job['fragment_end']
        read_length = job['read_length']

        bag = dict()

        for i in range(0, len(fragment) - read_length + 1, 1):
            seed_start = i + fragment_start
            seed_end = seed_start + read_length
            seed = fragment[i:i + read_length]

            zeed = zlib.compress(seed.encode())

            if zeed not in bag:
                bag[zeed] = list()
            bag[zeed].append((chrom, '+', seed_start, seed_end))

            seed_revcom = seed[::-1].translate(COMPLEMENT_DNA)

            zeed_revcom = zlib.compress(seed_revcom.encode())

            if zeed_revcom not in bag:
                bag[zeed_revcom] = list()
            bag[zeed_revcom].append((chrom, '-', seed_start, seed_end))

        return bag

    bag = dict()

    n_total = genome_length * read_length
    n_done = 0
    with ProcessPool(n_threads) as pool:
        for result in pool.uimap(
                process,
                new_jobs(
                    genome,
                    fragment_length,
                    read_length,
                ),
        ):
            bag.update(result)

            print(f'{n_done}', end='\r', flush=True)
            # print(f'{n_done / n_total}', end='\r', flush=True)
            # print(f'{(n_done / n_total) * 100:0.9f}%',
            #       end='\r',
            #       flush=True)
    return bag


def _align(queries, genome_index, n_threads):

    def jobs(queries, genome_index):
        for record in queries:
            ic(record['id'])
            yield {
                'id': record['id'],
                'read1': record['read1'],
                'read2': record['read2'],
                'insert_size': record['insert_size'],
                'genome_index': genome_index,
            }

    def process(job):
        id_ = job['id']
        ic(id_)
        read1 = job['read1']
        read2 = job['read2']
        insert_size = job['insert_size']
        genome_index = job['genome_index']

        zread1 = zlib.compress(read1.encode())
        zread2 = zlib.compress(read2.encode())

        if zread1 not in genome_index:
            return [hit(id_, read1, read2, insert_size)]
        if zread2 not in genome_index:
            return [hit(id_, read1, read2, insert_size)]

        hits1 = genome_index[zread1]
        hits2 = genome_index[zread2]

        if len(hits1) == len(hits2) == 1:
            ic()

            hit1 = hits1[0]
            hit2 = hits2[0]
            hit1_chrom = hit1[0]
            hit2_chrom = hit2[0]
            ic(hit1)
            ic(hit2)

            if hit1_chrom != hit2_chrom:
                return [hit(id_, read1, read2, insert_size)]

            chrom = hit1_chrom

            hit1_strand = hit1[1]
            hit2_strand = hit2[1]
            if hit1_strand != hit2_strand:
                return [hit(id_, read1, read2, insert_size)]

            strand = hit1_strand

            read1_start = hit1[2] + 1
            read1_end = hit1[3]
            read2_start = hit2[2] + 1
            read2_end = hit2[3]

            if strand == '+':
                actual_insert_size = read2_end - read1_start + 1
            else:
                actual_insert_size = read1_end - read2_start + 1

            if actual_insert_size != insert_size:
                return [hit(id_, read1, read2, insert_size)]

            return [
                hit(
                    id_=id_,
                    read1=read1,
                    read2=read2,
                    insert_size=insert_size,
                    chrom=chrom,
                    strand=strand,
                    read1_start=read1_start,
                    read1_end=read1_end,
                    read2_start=read2_start,
                    read2_end=read2_end,
                )
            ]

        hits = list()

        ic()
        for hit1 in hits1:
            chrom = hit1[0]
            strand = hit1[1]
            read1_start = hit1[2]
            read1_end = hit1[3]

            bag = list()
            for hit2 in hits2:
                if hit2[0] != chrom:
                    continue
                if hit2[1] != strand:
                    continue
                read2_end = hit2[3]

                if read2_end - read1_start == insert_size:
                    bag.append(hit2)

            for item in bag:
                read2_start = item[2]
                read2_end = item[3]

                hits.append(
                    hit(
                        id_=id_,
                        read1=read1,
                        read2=read2,
                        insert_size=insert_size,
                        chrom=chrom,
                        strand=strand,
                        read1_start=read1_start,
                        read1_end=read1_end,
                        read2_start=read2_start,
                        read2_end=read2_end,
                    ))

        return hits

    bag = list()

    for job in jobs(queries, genome_index):
        result = process(job)
        bag.extend(result)

    # else:
    #     with ProcessPool(n_threads) as pool:
    #         for result in pool.uimap(process, jobs(queries, genome_index)):
    #             bag.extend(result)

    result = pl.from_dicts(bag)

    return result


def hit(
    id_,
    read1,
    read2,
    insert_size,
    chrom=None,
    strand=None,
    read1_start=None,
    read1_end=None,
    read2_start=None,
    read2_end=None,
):
    return {
        'id': id_,
        'read1': read1,
        'read2': read2,
        'insert_size': insert_size,
        'chrom': chrom,
        'read1_start': read1_start,
        'read1_end': read1_end,
        'read2_start': read2_start,
        'read2_end': read2_end,
        'strand': strand,
    }


def load_queries(query_file):
    bag = []
    with query_file.open('rt') as fh:
        next(fh)    # skip header
        for line in fh:
            items = line.strip().split('\t')
            bag.append({
                'id': items[0],
                'read1': items[1],
                'read2': items[2],
                'insert_size': int(items[3]),
            })

    return bag


def chop_fasta(data):

    bag = dict()

    chrom = None

    while True:
        idx_comment = data.find('>')

        idx_seq_start = data.find('\n')

        chrom = data[idx_comment + 1:idx_seq_start]
        data = data[idx_seq_start + 1:]

        idx_seq_end = data.find('>')

        if idx_seq_end == -1:
            bag[chrom] = data.upper().replace('\n', '')
            break
        else:
            seq = data[0:idx_seq_end]
            bag[chrom] = seq.upper().replace('\n', '')

        data = data[idx_seq_end:]

    return bag


def load_genome(genome_file):

    with gzip.open(genome_file, 'rt') as fh:
        contents = fh.read()

        return chop_fasta(contents)


def get_genome_size(genome):
    n = 0

    for chrom, seq in genome.items():
        n += len(seq)

    return n


def merge(bag1, bag2):

    for k, v in bag2.items():
        if k in bag1:
            bag1[k].extend(v)
        else:
            bag1[k] = bag2[k]
    return bag1


def load(file_):
    with gzip.open(file_, 'rb') as fh:
        return pickle.load(fh)


def save(obj, file_):
    with gzip.open(file_, 'wb') as fh:
        pickle.dump(obj, fh)
