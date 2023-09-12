import re

REGULAR_BASE = re.compile(r'^[ACGT]+$')

from .genome import Genome


class Variant():

    def __init__(
        self,
        chrom: str,
        pos: int,
        ref: str,
        alt: str,
        id_: str = None,
        into: str = None,
    ):

        self.id = id_
        self.chrom = chrom
        self.pos = pos
        self.ref = ref
        self.alt = alt
        self.info = info

        if not REGULAR_BASE.match(self.ref):
            print(
                f'[WARN]\tWEIRD REF: {self.chrom}:{self.pos} {self.ref}>{self.alt}'
            )

        for alt in self.alt.split(','):
            if not REGULAR_BASE.match(alt):
                print(
                    f'[WARN]\tWEIRD ALT: {self.chrom}:{self.pos} {self.ref}>{self.alt}'
                )

        self.region = get_region(
            self.pos,
            self.ref,
            self.alt,
        )

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
            info=self.info,
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
            info=self.info,
        )

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
            info=self.info,
        )

    def __eq__(self, other):
        if self is other:
            return True

        return (self.chrom == other.chrom and self.pos == other.pos
                and self.id == other.id and self.ref == other.ref
                and self.alt == other.alt)

    def __hash__(self):
        return hash(self.chrom, self.pos, self.id, self.ref, self.alt)

    def __str__(self):
        return f'{self.chrom}\t{self.pos}\t{self.id}\t{self.ref}\t{self.alt}'

    def __str__(self):
        return f'{self.chrom}\t{self.pos}\t{self.id}\t{self.ref}\t{self.alt}\t.\t.\t{self.info}'

    def __repr__(self):
        return f'{self.chrom}\t{self.pos}\t{self.id}\t{self.ref}\t{self.alt}\t.\t.\t{self.info}'


def get_region(pos, ref, alt):
    if is_snv(ref, alt):
        return Region(pos, pos)
    if is_ins(ref, alt):
        return Region(pos, pos + 1)
    if is_del(ref, alt):
        return Region(pos + 1, pos + len(ref) - 1)
    if is_ma(ref, alt):
        return Region(pos, pos)
    if is_mnp(ref, alt):
        return Region(pos, pos + len(ref) - 1)
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


def is_mnp(ref, alt):
    if Variant.is_snv(ref, alt):
        return False
    if Variant.is_ins(ref, alt):
        return False

    if Variant.is_del(ref, alt):
        return False

    if Variant.is_ma(ref, alt):
        return False

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
