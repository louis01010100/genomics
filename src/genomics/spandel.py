
from .utils import copy_vcf_header, is_gzip, create_col2idx
import gzip
import logging

def _group_spandel(ifh, ofh):
    for line in ifh:
        if line.startswith('##'):
            continue
        col2idx = create_col2idx(line)
        break

    parent = None
    children = list()

    for line in ifh:
        line = line.strip()
        items = line.split('\t', 8)

        if '*' not in items[col2idx['ALT']]:
            if parent and children:
                
                record = SpanFamily(parent.split('\t'), children, col2idx).combine()
                ofh.write(f'{record}\n')
                parent = None
                children = list()
                ofh.write(f'{line}\n')
            elif parent: 
                ofh.write(f'{parent}\n')
                ofh.write(f'{line}\n')
                parent = None
            else:
                if len(items[col2idx['REF']]) > 1:
                    parent = line
                else:
                    ofh.write(f'{line}\n')

        else:
            children.append(line.split('\t'))

    if parent and children:
        record = SpanFamily(parent.split('\t'), children, col2idx).combine()
        ofh.write(f'{record}\n')
    elif parent: 
        ofh.write(f'{parent}\n')
    else:
        pass

def group_spandel(input_vcf_file, output_vcf_file):
    copy_vcf_header(input_vcf_file, output_vcf_file)

    if is_gzip(input_vcf_file):
        with gzip.open(input_vcf_file, 'rt') as ifh, output_vcf_file.open('at') as ofh:

            _group_spandel(ifh, ofh)
    else:
        with input_vcf_file.open('rt') as ifh, output_vcf_file.open('at') as ofh:
            _group_spandel(ifh, ofh)


    return output_vcf_file

def _expand_spandel(parent: list, children:list[list], col2idx:dict):
    dpos = int(parent[col2idx['POS']])
    dref = parent[col2idx['REF']]
    dalts = parent[col2idx['ALT']].split('\t')

    max_end = dpos + len(dref) -1
    max_ref = dref
    max_alts = dalts
    suffix = ''

    for child in children:
        tpos = int(child[col2idx['POS']])
        tref = child[col2idx['REF']]

        tend = tpos + len(tref) - 1

        if tend > max_end:
            suffix = tref[-(tend - max_end):]
            max_ref = max_ref + suffix
            max_alts = [max_alt + suffix for max_alt in max_alts]

            max_end = tpos + len(tref) -1

    parent_expanded = list()
    parent_expanded.extend(parent)

    parent_expanded[col2idx['ID']] = parent[col2idx['ID']]
    parent_expanded[col2idx['REF']] = max_ref
    parent_expanded[col2idx['ALT']] = ','.join(max_alts)
    return parent_expanded

class AlleleTranslator():

    def __init__(self):
        self._code2allele = {}
        self._allele2allele = {}
        self._allele2code = {}
        self._new_alts = None

    def __str__(self):
        return f'{self._code2allele}\n{self._allele2allele}\n{self._allele2code}\n'
    def __expr__(self):
        return str(self)

    def turn_code2allele(self, id_, code):
        return self._code2allele[id_][code]

    def turn_allele2code(self, allele):
        return self._allele2code[allele]

    def turn_allele2allele(self, id_, allele):
        return self._allele2allele[id_][allele]

    def update_code2allele(self, parent, children, col2idx):
        self._add_code2allele(parent[col2idx['ID']], 0, parent[col2idx['REF']])
        for code, allele in enumerate(parent[col2idx['ALT']].split(','), 1):
            self._add_code2allele(parent[col2idx['ID']], code, allele)

        for child in children:
            self._add_code2allele(child[col2idx['ID']], 0, child[col2idx['REF']])
            for code, allele in enumerate(child[col2idx['ALT']].split(','), 1):
                self._add_code2allele(child[col2idx['ID']], code, allele)

    def _add_code2allele(self, id_, code, allele):
        if id_ not in self._code2allele:
            self._code2allele[id_] = dict()

        self._code2allele[id_][code] = allele

    def update_allele2code(self, code, allele):
        self._allele2code[code] = allele

    def update_allele2allele(self, id_, old_allele, new_allele):
        if id_ not in self._allele2allele:
            self._allele2allele[id_] = dict()
        self._allele2allele[id_][old_allele] = new_allele

    @property
    def new_alts(self):
        return self._new_alts

    @new_alts.setter
    def new_alts(self, new_alts):
        self._new_alts = new_alts

    def translate(self, id_, gt):

        if '/' in gt:
            sep = '/'
        elif '|' in gt:
            sep = '|'
        else:
            assert len(gt) == 0
            sep = '/'

        new_codes = list()

        for code in gt.split(sep):
            if code == '.':
                new_codes.append('.')
            else:
                code = int(code)
                old_allele = self._code2allele[id_][code]
                if old_allele == '*':
                    return None
                new_allele = self._allele2allele[id_][old_allele]
                new_code = self._allele2code[new_allele]
                new_codes.append(new_code)
        return sep.join([str(x) for x in new_codes])

def set_id(record: list, col2idx:dict)-> list:
    if record[col2idx['ID']] == '.':
        record[col2idx['ID']] = ':'.join([record[0], record[1], record[3], record[4]])

    return record

def _new_allele_translator(
        parent: list[str], parent_expanded: list[str], 
        children:list[list[str]], col2idx:dict):
    at = AlleleTranslator()

    at.update_code2allele(parent, children, col2idx)

    dpos_expanded = int(parent_expanded[col2idx['POS']])
    dref_expanded = parent_expanded[col2idx['REF']]
    dalts_expanded = parent_expanded[col2idx['ALT']].split(',')

    at.update_allele2allele(
        parent[col2idx['ID']], 
        parent[col2idx['REF']], 
        dref_expanded,
    )

    for old_allele, new_allele in zip(parent[col2idx['ALT']].split(','), dalts_expanded):
        at.update_allele2allele(
            parent[col2idx['ID']], 
            old_allele,
            new_allele,
        )

    new_alts = dalts_expanded

    for child in children:
        tpos = int(child[col2idx['POS']])
        tref = child[col2idx['REF']]
        talts = child[col2idx['ALT']].split(',')
        prefix = dref_expanded[0:(tpos - dpos_expanded + 1 - 1)]
        suffix = dref_expanded[tpos - dpos_expanded + 1 -1 + len(tref) + 1 - 1:]
        at.update_allele2allele(child[col2idx['ID']], tref, dref_expanded)
        for talt in talts:
            if talt == '*':
                continue
            else:
                new_talt = prefix + talt + suffix
                at.update_allele2allele(child[col2idx['ID']], talt, new_talt)
            new_alts.append(new_talt)

    new_alts = sorted(new_alts)

    at.update_allele2code(dref_expanded, 0)
    for code, new_alt in enumerate(new_alts, 1):
        at.update_allele2code(new_alt, code)

    at.new_alts = new_alts

    return at



def __group_spandel(parent: list[str], children:list[list[str]], col2idx:dict) -> str:
    return SpanFamily(parent, children, col2idx).combine()


def update_calls(parent_expanded: list, children: list[list], col2idx, allele_translator)-> str:

    parent_expanded[col2idx['ALT']] = ','.join(allele_translator.new_alts)

    ref_allele = parent_expanded[col2idx['REF']]

    result = list()

    for idx in range(0, len(parent_expanded)):
        if idx < 9:
            result.append(parent_expanded[idx])
            continue

        alt_alleles = list()

        old_codes = parent_expanded[idx].split('/')

        cn = len(old_codes)

        n_alts_expected = 0

        for code in old_codes:
            if code == '.':
                continue
            if code == '0':
                continue
            n_alts_expected += 1

            code = int(code)

            allele_old = allele_translator.turn_code2allele(parent_expanded[col2idx['ID']], code)
            allele_new = allele_translator.turn_allele2allele(parent_expanded[col2idx['ID']], allele_old)

            alt_alleles.append(allele_new)

        conflict = False
        for child in children:
            n_stars = 0
            for code in child[idx].split('/'):
                if code == '.':
                    continue
                if code == '0':
                    continue

                code = int(code)
                allele_old = allele_translator.turn_code2allele(child[col2idx['ID']], code)

                if allele_old == '*':
                    n_stars += 1
                    continue

                allele_new = allele_translator.turn_allele2allele(child[col2idx['ID']], allele_old)

                alt_alleles.append(allele_new)

            if n_alts_expected != n_stars:
                conflict = True

        alleles = alt_alleles
        if conflict:
            call = './.'
            print('## multiallelic genotypes', parent_expanded, alleles)
        else:
            assert len(alleles) <= cn, f'CN:{cn}\t{alleles}'

            if len(alleles) < cn:
                for i in range(cn - len(alleles)):
                    alleles.append(ref_allele)
            codes = list()
            for allele in alleles:
                code = allele_translator.turn_allele2code(allele)
                codes.append(code)
            call = '/'.join([str(x) for x in sorted(codes)])

        result.append(call)

    return '\t'.join([str(x) for x in result])

def is_prefix(allele_new, alleles):
    for allele in alleles:
        if not allele.startswith(allele_new):
            return False
    return True
        


class SpanFamily:
    def __init__(self, parent, children, col2idx):
        self._raw_parent = parent
        self._raw_children = children
        self._col2idx = col2idx

    def combine(self):
        self._parent = set_id(self._raw_parent, self._col2idx)
        self._children = [set_id(child, self._col2idx) for child in self._raw_children]

        self._parent_expanded = _expand_spandel(self._parent, self._children, self._col2idx)

        at = _new_allele_translator(self._parent, self._parent_expanded, self._children, self._col2idx)

        grouped_record = update_calls(self._parent_expanded, self._children, self._col2idx, at)

        return grouped_record












        





