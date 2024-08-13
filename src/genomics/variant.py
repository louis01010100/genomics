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
    def calls(self) -> str:
        return self._calls

    @calls.setter
    def calls(self, calls: str):
        self._calls = calls

    @property
    def region(self):
        return self._region

    @property
    def alts(self) -> list:
        return self.alt.split(',')

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

    def align(self, ref_variant, genome):

        result = align(
                self.chrom, 
                self.pos, 
                self.ref, 
                self.alts, 
                ref_variant.pos, 
                ref_variant.ref, 
                ref_variant.alts,
                genome,
        )

        if result is None:
            return result

        new_chrom = result['chrom']
        new_pos = result['pos']
        new_ref = result['ref']
        new_alts = result['alts']
        merged_alts = merge_alts(ref_variant.alts, new_alts)

        new_calls = transcode_gt(
            idx2allele=_load_idx2allele(self.ref, self.alts),
            allele2allele=_load_allele2allele(
                [self.ref, *self.alts],
                [new_ref, *new_alts],
            ),
            allele2idx=_load_allele2idx(new_ref, merged_alts),
            calls=self.calls,
        )

        return Variant(
            chrom = new_chrom,
            pos = new_pos,
            ref = new_ref,
            alt = ','.join(merged_alts),
            id_ = self.id,
            qual = self.qual,
            filter_ = self.filter,
            info = self.info,
            format_ = self.format,
            calls = new_calls
        )



    def normalize(self, genome):
        result = normalize(self.chrom, self.pos, self.ref, self.alts, genome)

        return Variant(
            chrom = result['chrom'],
            pos = result['pos'],
            ref = result['ref'],
            alt = ','.join(result['alts']),
            id_ = self.id,
            qual = self.qual,
            filter_ = self.filter,
            info = self.info,
            format_ = self.format,
            calls = self.calls,
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

    def same_coordinate(self, other):

        return (self.chrom == other.chrom 
                and self.pos == other.pos 
                and self.ref == other.ref)

    def same_site(self, other, check_id = True):

        if check_id and self.id != other.id:
            return False

        return (self.chrom == other.chrom 
                and self.pos == other.pos 
                and self.ref == other.ref 
                and self.alt == other.alt)


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


def _load_allele2idx(ref: str, alts: list):

    allele_dict = {
        allele: str(index + 1)
        for index, allele in enumerate(
            [allele for allele in alts if allele != '.'])
    }
    allele_dict[ref] = '0'
    allele_dict['.'] = '.'

    return allele_dict


def align(chrom, pos, ref, alts, ref_pos, ref_ref, ref_alts, genome):

    if pos >= ref_pos + len(ref_ref):
        return None

    if pos > ref_pos:
        pos -= 1
        base = genome.slice(chrom, pos -1, pos)
        ref = base + ref
        alts = [base + alt for alt in alts]

        return align(chrom, pos, ref, alts, ref_pos, ref_ref, ref_alts, genome)

    elif pos < ref_pos:

        if all([alt[0] == ref[0] for alt in alts]):
            ref = ref[1:]
            alts = [alt[1:] for alt in alts]
        else:
            return None

        pos += 1
        base = genome.slice(chrom, pos -1, pos)
        ref = ref + base
        alts = [alt + base for alt in alts]
        return align(chrom, pos, ref, alts, ref_pos, ref_ref, ref_alts, genome)
    else:
        assert pos == ref_pos

        if ref == ref_ref:
            return {'chrom': chrom, 'pos': pos, 'ref': ref, 'alts': alts, 'status': 'PASS'}

        elif len(ref) < len(ref_ref):
            suffix = ref_ref[len(ref):]
            ref = ref + suffix
            alts = [alt + suffix for alt in alts]
            return align(chrom, pos, ref, alts, ref_pos, ref_ref, ref_alts, genome)
        else:
            return None

def rshift(chrom, pos, ref, alts, genome):

    # trim_left
    if all([alt[0] == ref[0] for alt in alts]):

        ref = ref[1:]
        alts = [alt[1:] for alt in alts]
        pos += 1

        if len(ref) == 0 or any([len(alt) == 0 for alt in alts]):
            base = genome.slice(chrom, pos -1, pos)

            if len(ref) == 0:
                ref = base

            bag = list()
            for alt in alts:
                if len(alt) == 0:
                    bag.append(base)
                else:
                    bag.append(alt)


        return rshift(chrom, pos, ref, alts, genome)

    # trim right
    
    if all([alt[-1] == ref[-1] for alt in alts]):

        pos += 1
        end = pos + len(ref)
        base = genome.slice(chrom, end -1, end)

        ref = ref[1:] + base
        alts = [alt[1:] + base for alt in alts]


    return { 'chrom': chrom , 'pos': pos, 'ref': ref, 'alts': alts }



    

def normalize(chrom, pos, ref, alts, genome):

    # trim right
    if all([alt[-1] == ref[-1] for alt in alts]):

        ref = ref[0:-1]
        alts = [alt[0:-1] for alt in alts]

        if len(ref) == 0 or any([len(alt) == 0 for alt in alts]):
            pos -= 1
            base = genome.slice(chrom, pos -1, pos)
            ref = base + ref
            alts = [base + alt for alt in alts]

        return normalize(chrom, pos, ref, alts, genome)

    
    # trim left
    if any([alt[0] != ref[0] for alt in alts]):

        return { 'chrom': chrom , 'pos': pos, 'ref': ref, 'alts': alts }

    if len(ref) == 1 or any([len(alt) == 1 for alt in alts]):
        return { 'chrom': chrom , 'pos': pos, 'ref': ref, 'alts': alts }

    ref = ref[1:]
    alts = [alt[1:] for alt in alts]
    pos += 1

    return normalize(chrom, pos, ref, alts, genome)



def merge_alts(ref_alts:list, alts:list) -> list:


    values = [*ref_alts]

    for alt in alts:
        if alt not in values:
            values.append(alt)
    return values



def _load_idx2allele(ref: str, alts: list):
    allele_dict = {
        str(index + 1): allele
        for index, allele in enumerate(
            [allele for allele in alts if allele != '.'])
    }

    allele_dict['0'] = ref
    allele_dict['.'] = '.'

    return allele_dict
