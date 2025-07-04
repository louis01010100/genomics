from pathlib import Path
from ncls import NCLS

class GenomicRegionDatabase():
    def __init__(self, regions, registory):
        self._regions = regions
        self._registory = registory

    @property
    def chroms(self):
        return set([str(x) for x in self._regions.keys()])

    def find_overlap(self, region:dict):

        bag = list()
        if region['chrom'] not in self._regions.keys():
            return bag

        records = self._regions[region['chrom']].find_overlap(region['start'], region['end'])

        for record in records:
            idx = record[2]
            bag.append(self._registory[idx])

        return bag

def create_database(records: list[dict]):

    registry = dict()

    idx = 0

    bag = dict()

    for record in records:
        chrom = record['chrom']
        start = record['start']
        end = record['end']

        # assert type(start) == int, record
        # assert type(end) == int, record

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

    return GenomicRegionDatabase(bag2, registry)

class GenomicRegion():

    def __init__(self, chrom, start, end, name = None):
        self._chrom = chrom
        self._start = start
        self._end = end
        self._name = name

    def intersects(self, other):
        if not self.overlaps(other):
            return None

        start = max(self.start, other.start)
        end = min(self.end, other.end)

        return GenomicRegion(self.chrom, start, end)

    def contains(self, chrom, pos):
        if self.chrom != chrom:
            return False
        return self.start <= pos and self.end >= pos

    def get_reciprocal_overlap(self, other):
        return _get_reciprocal_overlap(self, other)

    def get_max_boundary_difference(self, other):
        if self.chrom != other.chrom:
            return None

        return max(abs(self.start - other.start), abs(self.end - other.end))

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

    def __len__(self):
        return self._end - self._start

    @property
    def start(self):
        return self._start

    @property
    def end(self):
        return self._end
    @property
    def name(self):
        return self._name

    def __eq__(self, other):
        return self.chrom == other.chrom and self.start == other.start and self.end == other.end

    def __str__(self):
        return f'{self.chrom}:{self.start}-{self.end}'

    def __repr__(self):
        return self.__str__()




def _get_reciprocal_overlap(region_1, region_2):

    if region_1.chrom != region_2.chrom:
        return None

    intersect = region_1.intersects(region_2)

    if intersect is None:
        return 0

    return max (len(intersect) / len(region_1), len(intersect) / len(region_2))
