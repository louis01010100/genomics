class Genome():

    def __init__(self, genome_fh):
        chroms = {}
        id_ = None
        seq = None

        for line in genome_fh:
            line = line.strip()
            if line.startswith('>'):
                if id_:
                    chroms[id_] = seq
                match = re.match(r'>([^ ]+)(:? .+)?$', line)
                id_ = match.group(1)
                seq = ''
            else:
                seq += line.upper()
        if id_:
            chroms[id_] = seq

        self._chroms = chroms

    def slice(self, chrom, start, stop):
        return self._chroms[chrom][start:stop]

    def length(self, chrom):
        return len(self._chroms[chrom])
