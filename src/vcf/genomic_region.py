class GenomicRegion():
    def __init__(self, chrom, start, stop):
        self._chrom = chrom
        self._start = start
        self._stop = stop

    def overlap(self, other):
        if self.chrom != other.chrom:
            return False
        return self.start <= other.stop and self.stop >= other.start

    def merge(self, other):
        if self.chrom != other.chrom:
            raise Exception('{self} and {other} cannot be merged')
        return GenomicRegion(self.chrom, min(self.start, other.start), max(self.stop, other.stop))

    @property
    def chrom(self):
        return self._chrom

    @property
    def start(self):
        return self._start

    @property
    def stop(self):
        return self._stop

    def __eq__(self, other):
        return self.chrom == other.chrom and self.start == other.start and self.stop == other.stop

    def __str__(self):
        return f'{self.chrom}:{self.start}-{self.stop}'

    def __repr__(self):
        return self.__str__()
