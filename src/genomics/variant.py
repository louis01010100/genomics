class Variant():

    def __init__(self, chrom, pos, ref, alt, id_='.', data=''):
        self.chrom = chrom
        self.pos = pos
        self.id = id_
        self.ref = ref
        self.alt = alt
        self._alts = set()
        self._data = data

    @property
    def data(self):
        return self._data

    @property
    def alts(self):
        if not self._alts:
            self._alts = set(self.alt.split(','))
        return self._alts

    @property
    def is_mnv(self):
        mnv = True
        for alt in self._alts:
            if len(alt) != len(self.ref):
                mnv = False
        return mnv

    def __str__(self):

        return f'{self.chrom}\t{self.id}\t{self.pos}\t{self.ref}\t{self.alt}'

    def __repr__(self):
        return self.__str__()

    def __hash__(self):
        return hash((self.chrom, self.pos, self.ref, self.alt))

    def __eq__(self, other):
        if self is other:
            return True

        if self.chrom != other.chrom:
            return False

        if self.pos != other.pos:
            return False

        if self.ref != other.ref:
            return False

        if self.alt != other.alt:
            return False

        return True
