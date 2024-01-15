from .gregion import GenomicRegion

HG38_X_PAR_1 = GenomicRegion('chrX', 10000, 2781479)
HG38_X_PAR_2 = GenomicRegion('chrX', 2781479, 156030895)
HG38_Y_PAR_1 = GenomicRegion('chrY', 10000, 2781479)
HG38_Y_PAR_2 = GenomicRegion('chrY', 56887902, 57217415)

HG19_X_PAR_1 = GenomicRegion('chrX', 60001, 2699520)
HG19_X_PAR_2 = GenomicRegion('chrX', 154931044, 155260560)
HG19_Y_PAR_1 = GenomicRegion('chrY', 10001, 2649520)
HG19_Y_PAR_2 = GenomicRegion('chrY', 59034050, 59363566)


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

    @property
    def version(self):
        return self._version

    def slice(self, chrom, start, stop):
        return self._chroms[chrom][start:stop]

    def length(self, chrom):
        return len(self._chroms[chrom])

    # def ploidy(self, chrom, pos):
    #     v = GenomicRegion(chrom.upper, pos, pos)
