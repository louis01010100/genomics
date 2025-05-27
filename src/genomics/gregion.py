from pathlib import Path
from ncls import NCLS

class GenomicRegions():
    def __init__(self, regions, registory):
        self._regions = regions
        self._registory = registory

    @property
    def chroms(self):
        return set([str(x) for x in self._regions.keys()])

    def find_overlap(self, chrom:str, start:int, end:int):

        bag = list()
        if chrom not in self._regions.keys():
            return bag

        records = self._regions[chrom].find_overlap(start, end)

        for record in records:
            idx = record[2]
            bag.append(self._registory[idx])

        return bag

def create_genomic_regions(records):

    registry = dict()

    idx = 0

    bag = dict()

    for record in records:
        chrom = record['chrom']
        start = record['start']
        end = record['end']

        registry[idx] = record

        if chrom not in bag:
            bag[chrom] = {'idx': list(), 'starts': list(), 'ends': list()}

        bag[chrom]['idx'].append(idx)
        bag[chrom]['starts'].append(start)
        bag[chrom]['ends'].append(end)

        idx += 1

    bag2 = dict()

    for k, v in bag.items():
        bag2[k] = NCLS(v['starts'], v['ends'], v['idx'])

    return GenomicRegions(bag2, registry)

class GenomicRegion():

    def __init__(self, chrom, start, end, name = None):
        self._chrom = chrom
        self._start = start
        self._stop = end
        self._name = name

    def intersects(self, other):
        if not self.overlaps(other):
            return None

        start = min(self.start, other.start)
        end = min(self.end, other.end)

        return GenomicRegion(self.chrom, start, end)

    def contains(self, chrom, pos):
        if self.chrom != chrom:
            return False
        return self.start <= pos and self.end >= pos

    def overlaps(self, other):
        if self.chrom != other.chrom:
            return False
        return self.start < other.end and self.end > other.start


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
    @property
    def name(self):
        return self._name

    def __eq__(self, other):
        return self.chrom == other.chrom and self.start == other.start and self.end == other.end

    def __str__(self):
        return f'{self.chrom}:{self.start}-{self.end}'

    def __repr__(self):
        return self.__str__()
