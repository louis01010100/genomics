class GenomicRegion():

    def __init__(self, chrom, start, end):
        self._chrom = chrom
        self._start = start
        self._stop = end

    def contains(self, chrom, pos):
        if self.chrom != chrom:
            return False
        return self.start <= pos and self.end >= pos

    def overlaps(self, other):
        if self.chrom != other.chrom:
            return False
        return self.start <= other.end and self.end >= other.start


    def merge(self, other):
        if self.chrom != other.chrom:
            raise Exception('{self} and {other} cannot be merged')
        return GenomicRegion(self.chrom, min(self.start, other.start),
                             max(self.end, other.end))

    @property
    def chrom(self):
        return self._chrom

    @property
    def start(self):
        return self._start

    @property
    def end(self):
        return self._stop

    def __eq__(self, other):
        return self.chrom == other.chrom and self.start == other.start and self.end == other.end

    def __str__(self):
        return f'{self.chrom}:{self.start}-{self.end}'

    def __repr__(self):
        return self.__str__()
