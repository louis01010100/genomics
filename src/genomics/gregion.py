from pathlib import Path
from ncls import NCLS

class GenomicRegions():
    def __init__(self, regions, idx2name):
        self.regions = regions
        self.idx2name = idx2name

    @property
    def chroms(self):
        return set([str(x) for x in self.regions.keys()])

    def find_overlap(self, chrom:str, start:int, end:int):

        bag = list()
        if chrom not in self.regions.keys():
            return bag

        records = self.regions[chrom].find_overlap(start, end)

        for record in records:
            start = record[0]
            end = record[1]
            idx = record[2]
            name = self.idx2name[idx]
            bag.append(GenomicRegion(chrom, start, end, name))

        return bag


def create_genomic_regions(df):


    bag = dict()

    idx = 0
    idx2name = dict()

    for record in df.to_dicts():
        name = record['name']
        chrom = record['chrom']
        start = int(record['start'])
        end = int(record['end'])

        if chrom not in bag:
            bag[chrom] = {'idx': list(), 'starts': list(), 'ends': list()}

        bag[chrom]['idx'].append(idx)
        bag[chrom]['starts'].append(start)
        bag[chrom]['ends'].append(end)

        idx2name[idx]= name
        idx += 1

    

    bag2 = dict()

    for k, v in bag.items():
        bag2[k] = NCLS(v['starts'], v['ends'], v['idx'])

    return GenomicRegions( bag2, idx2name)


class GenomicRegion():

    def __init__(self, chrom, start, end, name = None):
        self._chrom = chrom
        self._start = start
        self._stop = end
        self._name = name

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
    @property
    def name(self):
        return self._name

    def __eq__(self, other):
        return self.chrom == other.chrom and self.start == other.start and self.end == other.end

    def __str__(self):
        return f'{self.chrom}:{self.start}-{self.end}'

    def __repr__(self):
        return self.__str__()
