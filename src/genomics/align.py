COMPLEMENT_DNA = str.maketrans('ACGT', 'TGCA')


def align(query_file, genome_file, output_dir, read_length: int, n_threads):

    genome = load_genome(genome_file)
    genome_index = index(genome, read_length)

    queries = load_queries(query_file)

    with ProcessPool(n_threads) as pool:
        for result in pool.uimap(process, jobs(queries, genome_index)):
            pass
def hit(probeset_id, read1, read2, insert_size, chrom = None, strand = None, read1_start = None, read1_end = None, read2_start = None, read2_end = None,):
    return {
            'probeset_id': probeset_id,
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

def process(job):
    probeset_id = job['probeset_id']
    read1 =job['read1']
    read2 = job['read2']
    insert_size = job['insert_size']
    genome_index = job['genome_index']

    if read1 not in genome_index:
        return hit(probeset_id, read1, read2, insert_size)
    if read2 not in genome_index:
        return hit(probeset_id, read1, read2, insert_size)

    hits1 = genome_index[read1]
    hits2 = genome_index[read2]

    if len(hits1) == len(hits2) == 1:
        hit1 = hits1[0]
        hit2 = hits2[0]
        hit1_chrom = hit1[0]
        hit2_chrom = hit2[0]

        if hit1_chrom != hit2_chrom
            return hit(probeset_id, read1, read2, insert_size)

        chrom = hit1_chrom

        hit1_strand = hit1[1]
        hit2_strand = hit2[1]
        if hit1_strand != hit2_strand:
            return hit(probeset_id, read1, read2, insert_size)

        strand = hit1_strand

        read1_start = hit1[2] + 1
        read1_end = hit1[3] 
        read2_start = hit2[2] + 1
        read2_end = hit2[3] 

        if read2_end - read1_end + 1 != insert_size:
            return hit(probeset_id, read1, read2, insert_size)


        return hit(probeset_id = probeset_id, read1 = read1, read2 = read2, insert_size = insert_size, chrom = chrom, strand = strand, read1_start = read1_start, read1_end = read1_end, read2_start = read2_start, read2_end = read2_end,)

    for hit1 in hits1:
        chrom = hit1['chrom']
        strand = hit1['chrom']



def jobs(queries, genome_index):
    for record in queries.to_dicts():
        yield {
                'probeset_id': record['probeset_id'],
                'read1': record['read1'],
                'read2': record['read2'],
                'insert_size': record['insert_size'],
                'genome_index': genome_index,
                }


def load_queries(query_file):
    erturn pl.read_csv(query_file, has_header = True, separator = '\t')


def load_genome(genome_file):

    genome = dict()

    chrom = None

    tmp = []

    with genome_file.open('rt') as fh:
        for line in fh:
            line = line.strip()

            if line.startswith('>'):
                if not chrom:
                    chrom = line.split(' ')[0][1:]
                else:
                    genome[chrom] = ''.join(tmp)
                    chrom = line.split(' ')[0][1:]
            else:
                tmp.append(line)


def index(genome: dict, read_length):

    bag = dict()

    for chrom, seq in genome.items():
        for i in range(0, len(seq) - read_length + 1, 1):
            seed = seq[i, i + read_length]
            start = i
            end = i + read_length

            if seed not in bag:
                bag[seed] = list()
            bag[seed].append((chrom, '+', start, end))

            seed_revcom = seed[::-1].translate(COMPLEMENT_DNA)

            if seed_revcom not in bag:
                bag[seed_revcom] = list()
            bag[seed_revcom].append((chrom, '-', start, end))

    return bag
