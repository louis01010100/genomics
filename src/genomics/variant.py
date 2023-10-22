import re
from icecream import ic

REGULAR_BASE = re.compile(r'^[ACGT]+$')

from .genome import Genome
from .gregion import GenomicRegion


class Variant():

    def __init__(
        self,
        chrom: str,
        pos: int,
        ref: str,
        alt: str,
        id_: str = None,
        data: str = None,
    ):

        self.id = id_
        self.chrom = chrom
        self.pos = pos
        self.ref = ref
        self.alt = alt
        self.data = data

        # if not REGULAR_BASE.match(self.ref):
        #     print(
        #         f'[WARN]\tWEIRD REF: {self.chrom}:{self.pos} {self.ref}>{self.alt}'
        #     )
        #
        # for alt in self.alt.split(','):
        #     if not REGULAR_BASE.match(alt):
        #         print(
        #             f'[WARN]\tWEIRD ALT: {self.chrom}:{self.pos} {self.ref}>{self.alt}'
        #         )

        self._region = get_region(
            self.chrom,
            self.pos,
            self.ref,
            self.alt,
        )

    @property
    def region(self):
        return self._region

    @property
    def alts(self):
        return self.alt.split(',')

    def fix_ref(self, genome: Genome):
        g_ref = genome.slice(self.chrom, self.pos - 1,
                             self.pos - 1 + len(self.ref))
        if g_ref != self.ref:
            print(f'[WARN]\tFix Inconsistent REF: {self.ref}/{g_ref}')

        self.ref = g_ref

        return Variant(
            chrom=self.chrom,
            pos=self.pos,
            id_=self.id_,
            ref=g_ref,
            alt=self.alt,
            data=self.data,
        )

    def norm(self, genome, left_justified=True):
        g_ref = genome.slice(
            self.chrom,
            self.pos - 1,
            self.pos - 1 + len(self.ref),
        )
        if g_ref != self.ref:
            print(f'[WARN]\tInconsistent REF: {self.ref}/{g_ref}')

        self.ref = g_ref

        if left_justified:
            return self._norm_left(genome)
        return self._norm_right(genome)

    def _norm_right(self, genome):
        alts = [alt for alt in self.alt.split(',')]
        pos, ref, alts = trim_right_bases(
            self.chrom,
            self.pos,
            self.ref,
            alts,
            genome,
        )

        pos, ref, alts = trim_left_bases(
            pos,
            ref,
            alts,
        )

        return Variant(
            chrom=self.chrom,
            pos=pos,
            id=self.id,
            ref=ref,
            alt=','.join(alts),
            data=self.data,
        )

    def merge(self, other):

        if self == other:
            return Variant(
                chrom=self.chrom,
                pos=self.pos,
                id_=self.id,
                ref=self.ref,
                alt=self.alt,
                data=self.data,
            )

        sregion = self.region
        oregion = other.region

        if not sregion.is_overlapping(oregion):
            raise Exception(sregion, oregion)

        if self.pos == other.pos and self.ref == other.ref:
            bag = set()
            for alt in [*self.alts, *other.alts]:
                bag.add(alt)
            alt = ','.join(sorted(list(bag)))
            return Variant(self.chrom, self.pos, self.ref, alt)

        pos, s_ref, s_alts, o_ref, o_alts = sync_prefix(
            start1=sregion.start,
            ref1=self.ref,
            alts1=self.alts,
            start2=oregion.start,
            ref2=other.ref,
            alts2=other.alts,
        )

        s_ref, s_alts, o_ref, o_alts = sync_suffix(
            end1=sregion.end,
            ref1=s_ref,
            alts1=s_alts,
            end2=oregion.end,
            ref2=o_ref,
            alts2=o_alts,
        )

        assert s_ref == o_ref, f'{s_ref} != {o_ref}'

        ref = s_ref
        alts = ','.join(sorted(list(set([*s_alts, *o_alts]))))

        return Variant(self.chrom, pos, ref, alts)

    @property
    def is_snv(self):
        return is_snv(self.ref, self.alt)

    @property
    def is_ins(self):
        return is_ins(self.ref, self.alt)

    @property
    def is_del(self):
        return is_del(self.ref, self.alt)

    @property
    def is_ma(self):
        return is_ma(self.ref, self.alt)

    @property
    def is_mnv(self):
        return is_mnv(self.ref, self.alt)

    def _norm_left(self, genome):
        alts = [alt for alt in self.alt.split(',')]
        pos, ref, alts = trim_right_bases(
            self.chrom,
            self.pos,
            self.ref,
            alts,
            genome,
        )

        pos, ref, alts = trim_left_bases(
            pos,
            ref,
            alts,
        )

        return Variant(
            chrom=self.chrom,
            pos=pos,
            id=self.id,
            ref=ref,
            alt=','.join(alts),
            data=self.data,
        )

    def __eq__(self, other):
        if self is other:
            return True

        return (self.chrom == other.chrom and self.pos == other.pos
                and self.id == other.id and self.ref == other.ref
                and self.alt == other.alt)

    def __hash__(self):
        return hash((self.chrom, self.pos, self.id, self.ref, self.alt))

    def __str__(self):
        return f'{self.chrom}\t{self.pos}\t{self.id}\t{self.ref}\t{self.alt}\t.\t.\t{self.data}'

    def __repr__(self):
        return f'{self.chrom}\t{self.pos}\t{self.id}\t{self.ref}\t{self.alt}\t.\t.\t{self.data}'


def get_region(chrom, pos, ref, alt):
    if is_snv(ref, alt):
        return GenomicRegion(chrom, pos, pos)
    if is_ins(ref, alt):
        return GenomicRegion(chrom, pos, pos)
    if is_del(ref, alt):
        return GenomicRegion(chrom, pos, pos + len(ref) - 1)
    if is_ma(ref, alt):
        return GenomicRegion(chrom, pos, pos + len(ref) - 1)
    if is_mnv(ref, alt):
        return GenomicRegion(chrom, pos, pos + len(ref) - 1)
    raise Exception(f'{pos} {ref} {alt}')


def is_snv(ref, alt):

    if (len(ref) == len(alt) == 1):
        return True

    return False


def is_ins(ref, alt):

    if alt.startswith(ref) and len(ref) < len(alt):
        return True

    return False


def is_del(ref, alt):

    if ref.startswith(alt) and len(ref) > len(alt):
        return True

    return False


def is_ma(ref, alt):
    return ',' in alt


def is_mnv(ref, alt):
    if is_snv(ref, alt):
        return False
    if is_ins(ref, alt):
        return False

    if is_del(ref, alt):
        return False

    if is_ma(ref, alt):
        return False

    if len(ref) == len(alt):
        return True

    return True


def trim_right_bases(chrom, pos, ref, alts, genome):
    if len(ref) == 0:
        pos -= 1
        base = genome.slice(chrom, pos - 1, pos)
        ref = base
        alts = [base + x for x in alts]

        return trim_right_bases(
            chrom,
            pos,
            ref,
            alts,
            genome,
        )

    if any([len(x) == 0 for x in alts]):
        pos -= 1
        base = genome.slice(chrom, pos - 1, pos)
        ref = base + ref
        alts = [base + x for x in alts]

        return trim_right_bases(
            chrom,
            pos,
            ref,
            alts,
            genome,
        )

    ref_end_base = ref[-1]

    trim_base = True

    for alt in alts:
        if alt[-1] != ref_end_base:
            trim_base = False

    if trim_base:
        ref = ref[0:-1]
        alts = [x[0:-1] for x in alts]
        return trim_right_bases(
            chrom,
            pos,
            ref,
            alts,
            genome,
        )
    return pos, ref, alts


def trim_left_bases(pos, ref, alts):

    if len(ref) == 1:
        return pos, ref, alts

    if any([len(x) == 1 for x in alts]):
        return pos, ref, alts

    base = ref[0]

    trim = True

    for alt in alts:
        if alt[0] != base:
            trim = False

    if trim:
        ref = ref[1:]
        alts = [x[1:] for x in alts]
        pos += 1

        return trim_left_bases(pos, ref, alts)
    else:
        return pos, ref, alts


def is_vcf(pos_start, pos_stop, ref_allele, alt_allele):
    if ref_allele == '-':
        return False

    if alt_allele == '-':
        return False

    return True

    # # Multiallelic SNP is always in VCF format for Axiom arrays
    # if ',' in alt_allele:
    #     return True
    #
    #
    # #SNV/MNV
    # if pos_start == pos_stop and len(ref_allele) == len(alt_allele):
    #     return True
    #
    # # INDEL
    # if ref_allele.startswith(alt_allele):
    #     return True
    #
    # if alt_allele.startswith(ref_allele):
    #     return True
    #
    # return False


def sync_prefix(
    start1: int,
    ref1: str,
    alts1: list,
    start2: int,
    ref2: str,
    alts2: list,
):
    if start1 == start2:
        pos = start1
    elif start1 < start2:
        prefix = ref1[0:(start2 - start1)]
        ref2 = prefix + ref2
        alts2 = [prefix + x for x in alts2]
        pos = start1

    elif start1 > start2:
        prefix = ref2[0:(start1 - start2)]
        ref1 = prefix + ref1
        alts1 = [prefix + x for x in alts1]
        pos = start2
    else:
        raise Exception(start1, start2)

    return pos, ref1, alts1, ref2, alts2


def sync_suffix(
    end1: int,
    ref1: str,
    alts1: list,
    end2: int,
    ref2: str,
    alts2: list,
):

    if end1 == end2:
        pass
    elif end1 < end2:
        # -----
        # ---
        suffix = ref2[-(end2 - end1):]
        ref1 = ref1 + suffix
        alts1 = [x + suffix for x in alts1]
    elif end1 > end2:
        suffix = ref1[-(end1 - end2):]
        ref2 = ref2 + suffix
        alts2 = [x + suffix for x in alts2]
    else:
        raise Exception(end1, end2)

    return ref1, alts1, ref2, alts2
