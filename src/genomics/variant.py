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
        qual: str = None,
        filter_: str = None,
        info: str = None,
        format_: str = None,
        calls: str = None,
    ):

        self._chrom = chrom
        self._pos = int(pos)
        self._id = id_
        self._ref = ref.upper()
        self._alt = alt.upper()
        self._qual = qual
        self._filter = filter_
        self._info = info
        self._format = format_
        self._calls = calls

        if format_ is not None:
            assert calls is not None
        if calls is not None:
            assert format_ is not None

        self._region = get_region(
            self._chrom,
            self._pos,
            self._ref,
            self._alt,
        )

    def clone(self):
        return Variant(
            self.chrom,
            self.pos,
            self.ref,
            self.alt,
            self.id,
            self.qual,
            self.filter,
            self.info,
            self.format,
            self.calls,
        )

    @property
    def chrom(self):
        return self._chrom

    @property
    def pos(self):
        return self._pos

    @property
    def id(self):
        return self._id

    @id.setter
    def id(self, id_):
        self._id = id_

    @property
    def ref(self):
        return self._ref

    @property
    def alt(self):
        return self._alt

    @property
    def qual(self):
        return self._qual

    @property
    def filter(self):
        return self._filter

    @property
    def info(self):
        return self._info

    @property
    def format(self):
        return self._format

    @format.setter
    def format(self, format_):
        self._format = format_

    @property
    def calls(self):
        return self._calls

    @calls.setter
    def calls(self, calls):
        self._calls = calls

    @property
    def region(self):
        return self._region

    @property
    def alts(self) -> set:
        return set(self.alt.split(','))

    def is_overlapping(self, other):
        return self.region.is_overlapping(other.region)

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

        same_site = (self.chrom == other.chrom and self.pos == other.pos
                     and self.id == other.id and self.ref == other.ref
                     and self.alt == other.alt)

        if not same_site:
            return False

        if self.calls is None and other.calls is None:
            return True
        else:
            return self.calls == other.calls

    def __hash__(self):
        return hash(
            (self.chrom, self.pos, self.id, self.ref, self.alt, self.qual,
             self.filter, self.info, self.format, self.calls))

    def __str__(self):

        bag = list()

        bag.append(self.chrom)
        bag.append(str(self.pos))
        bag.append(self.id if self.id else '.')
        bag.append(self.ref)
        bag.append(self.alt)
        bag.append(self.qual if self.qual else '.')
        bag.append(self.filter if self.filter else '.')
        bag.append(self.info if self.info else '.')

        if self.format:
            bag.append(self.format)
        if self.calls:
            bag.append(self.calls)

        return '\t'.join(bag)

    def __repr__(self):
        return self.__str__()


def get_region(chrom, pos, ref, alt):

    if is_snv(ref, alt):
        return GenomicRegion(chrom, pos, pos)
    if is_ins(ref, alt):
        if len(ref) == 1:
            end = pos
        else:
            end = pos + len(ref) - 1
        return GenomicRegion(chrom, pos, end)
    if is_del(ref, alt):
        return GenomicRegion(chrom, pos, pos + len(ref) - 1)
    if is_ma(ref, alt):
        return GenomicRegion(chrom, pos, pos + len(ref) - 1)
    if is_mnv(ref, alt):
        return GenomicRegion(chrom, pos, pos + len(ref) - 1)
    raise Exception(f'{pos} {ref} {alt}')


def is_snv(ref, alt):

    if (len(ref) == len(alt) == 1) and not is_ma(ref, alt):
        return True

    return False


def is_ins(ref, alt):

    if alt.startswith(ref) and len(ref) < len(alt) and not is_ma(ref, alt):
        return True

    return False


def is_del(ref, alt):

    if ref.startswith(alt) and len(ref) > len(alt) and not is_ma(ref, alt):
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


def merge_alts(*alts):
    alts = set(*alts)
    alts.discard('.')
    if len(alts) == 0:
        return '.'
    return sorted(list(alts))


def merge_variant_id(*ids):
    id_ = sorted(
        list({id_
              for id_ in ids if (id_ is not None and id_ != '.')}))

    if len(id_) == 0:
        return None
    return ','.join(id_)


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


def transcode_gt(
    idx2allele: dict,
    allele2idx: dict,
    allele2allele: dict,
    calls: str,
):
    if not calls:
        return None

    calls = calls.split('\t')
    bag = list()
    for call in calls:
        bag.append(
            _transcode_gt(
                gt=call,
                idx2allele=idx2allele,
                allele2idx=allele2idx,
                allele2allele=allele2allele,
            ))

    return '\t'.join(bag)


def _transcode_gt(
    gt: str,
    idx2allele: dict,
    allele2allele: dict,
    allele2idx: dict,
):
    bag = list()

    delimiter = None
    if '/' in gt:
        delimiter = '/'
    elif '|' in gt:
        delimiter = '|'
    else:
        pass

    if delimiter:
        for idx in gt.split(delimiter):

            allele = idx2allele[idx]

            new_allele = allele2allele[allele]

            new_idx = allele2idx[new_allele]

            bag.append(new_idx)

        if '.' in bag:
            bag.remove('.')
            bag = sorted(bag)
            bag.append('.')
            new_gt = delimiter.join(bag)
        else:
            new_gt = delimiter.join(sorted(bag))

        return new_gt
    elif gt == '.':
        return '.'
    elif gt.isnumeric():
        idx = gt
        allele = idx2allele[idx]

        allele = allele2allele[allele]
        idx = allele2idx[allele]
        return idx
    else:
        raise Exception(gt)


def _load_allele2allele(old: list, new: list) -> dict:
    bag = dict()
    for x, y in zip(old, new):
        bag[x] = y
    bag['.'] = '.'

    return bag


def _load_allele2idx(ref: str, alt: str):

    allele_dict = {
        allele: str(index + 1)
        for index, allele in enumerate(
            [allele for allele in alt.split(',') if allele != '.'])
    }
    allele_dict[ref] = '0'
    allele_dict['.'] = '.'

    return allele_dict


def sync_alleles(
    v1: Variant,
    v2: Variant,
    merge_id: bool = False,
    merge_alt: bool = False,
):
    if not v1.is_overlapping(v2):
        return v1, v2

    if v1 == v2:
        v = Variant(
            chrom=v1.chrom,
            pos=v1.pos,
            id_=v1.id,
            ref=v1.ref,
            alt=v1.alt,
            qual=v1.qual,
            filter_=v1.filter,
            info=v1.info,
            calls=v1.calls,
        )

        return v, v

    sregion = v1.region
    oregion = v2.region

    new_pos, new_v1_ref, new_v1_alts, new_v2_ref, new_v2_alts = sync_prefix(
        start1=sregion.start,
        ref1=v1.ref,
        alts1=v1.alts,
        start2=oregion.start,
        ref2=v2.ref,
        alts2=v2.alts,
    )

    new_v1_ref, new_v1_alts, new_v2_ref, new_v2_alts = sync_suffix(
        end1=sregion.end,
        ref1=new_v1_ref,
        alts1=new_v1_alts,
        end2=oregion.end,
        ref2=new_v2_ref,
        alts2=new_v2_alts,
    )

    assert new_v1_ref == new_v2_ref, f'{new_v1_ref} != {new_v2_ref};{v1};{v2}'

    new_v1_pos = new_pos
    new_v2_pos = new_pos

    if merge_id:
        new_v1_id = merge_variant_id(v1.id, v2.id)
        new_v2_id = new_v1_id
    else:
        new_v1_id = v1.id
        new_v2_id = v2.id

    if merge_alt:
        new_v1_alt = ','.join(merge_alts([*new_v1_alts, *new_v2_alts]))
        new_v2_alt = new_v1_alt
    else:
        new_v1_alt = ','.join(sorted(new_v1_alts))
        new_v2_alt = ','.join(sorted(new_v2_alts))

    new_v1_format = v1.format
    new_v2_format = v2.format

    new_v1_calls = transcode_gt(
        idx2allele=_load_idx2allele(v1.ref, v1.alt),
        allele2idx=_load_allele2idx(new_v1_ref, new_v1_alt),
        allele2allele=_load_allele2allele(
            [v1.ref, *v1.alts],
            [new_v1_ref, *new_v1_alts],
        ),
        calls=v1.calls,
    )

    new_v2_calls = transcode_gt(
        idx2allele=_load_idx2allele(v2.ref, v2.alt),
        allele2idx=_load_allele2idx(new_v2_ref, new_v2_alt),
        allele2allele=_load_allele2allele(
            [v2.ref, *v2.alts],
            [new_v2_ref, *new_v2_alts],
        ),
        calls=v2.calls,
    )

    new_v1 = Variant(
        chrom=v1.chrom,
        pos=new_v1_pos,
        id_=new_v1_id,
        ref=new_v1_ref,
        alt=new_v1_alt,
        qual=v1.qual,
        filter_=v1.filter,
        info=v1.info,
        format_=new_v1_format,
        calls=new_v1_calls,
    )

    new_v2 = Variant(
        chrom=v2.chrom,
        pos=new_v2_pos,
        id_=new_v2_id,
        ref=new_v2_ref,
        alt=new_v2_alt,
        qual=v2.qual,
        filter_=v2.filter,
        info=v2.info,
        format_=new_v2_format,
        calls=new_v2_calls,
    )

    return new_v1, new_v2


def _load_idx2allele(ref: str, alt: str):
    allele_dict = {
        str(index + 1): allele
        for index, allele in enumerate(
            [allele for allele in alt.split(',') if allele != '.'])
    }

    allele_dict['0'] = ref
    allele_dict['.'] = '.'

    return allele_dict
