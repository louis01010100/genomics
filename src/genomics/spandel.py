
from .utils import copy_vcf_header, is_gzip, create_col2idx

def _group_spandel(ifh, ofh):
    for line in ifh:
        if line.startswith('##'):
            continue
        col2idx = create_col2idx(line)
        break

    deletion = None
    targets = list()

    for line in ifh:
        line = line.strip()
        items = line.split('\t', 8)

        if '*' not in items[col2idx['ALT']]:
            if deletion and targets:
                record = __group_spandel(deletion.split('\t'), targets, col2idx)
                ofh.write(f'{record}\n')
                deletion = None
                targets = list()
                ofh.write(f'{line}\n')
            elif deletion: 
                ofh.write(f'{deletion}\n')
                ofh.write(f'{line}\n')
                deletion = None
            else:
                if len(items[col2idx['REF']]) > 1:
                    deletion = line
                else:
                    ofh.write(f'{line}\n')

        else:
            targets.append(line.split('\t'))

    if deletion and targets:
        record = _group_spandel(deletion, targets, col2idx)
        ofh.write(f'{record}\n')
    elif deletion: 
        ofh.write(f'{deletion}\n')
    else:
        pass

def group_spandel(input_vcf_file, output_vcf_file):
    copy_vcf_header(input_vcf_file, output_vcf_file)

    if is_gzip(input_vcf_file):
        with gzip.open(input_vcf_file) as ifh, output_vcf_file.open('at') as ofh:

            _group_spandel(ifh, ofh)
    else:
        with input_vcf_file.open('rt') as ifh, output_vcf_file.open('at') as ofh:
            _group_spandel(ifh, ofh)


    return output_vcf_file

def _expand_spandel(deletion: list, targets:list[list], col2idx:dict)-> tuple[str, str]:
    dpos = int(deletion[col2idx['POS']])
    dref = deletion[col2idx['REF']]
    dalts = deletion[col2idx['ALT']].split('\t')

    max_end = dpos + len(dref) -1
    max_ref = dref
    max_alts = dalts
    suffix = ''

    for target in targets:
        tpos = int(target[col2idx['POS']])
        tref = target[col2idx['REF']]

        tend = tpos + len(tref) - 1

        if tend > max_end:
            suffix = tref[-(tend - max_end):]
            max_ref = max_ref + suffix
            max_alts = [max_alt + suffix for max_alt in max_alts]

            max_end = tpos + len(tref) -1

    expanded_deletion = list()
    expanded_deletion.extend(deletion)

    expanded_deletion[col2idx['REF']] = max_ref
    expanded_deletion[col2idx['ALT']] = ','.join(max_alts)
    return expanded_deletion

class AlleleTranslator():

    def __init__(self):
        self._code2allele = {}
        self._allele2allele = {}
        self._allele2code = {}

    def code2allele(self, id_, code):
        return self._code2allele[id_][code]

    def allele2code(self, allele):
        return self._allele2code[allele]

    def allele2allele(self, id_, allele):
        return self._allele2allele[id_][allele]

    def update_code2allele(self, deletion, targets, col2idx):
        self._add_code2allele(deletion[col2idx['ID']], 0, deletion[col2idx['REF']])
        for code, allele in enumerate(deletion[col2idx['ALT']].split(','), 1):
            self._add_code2allele(deletion[col2idx['ID']], code, allele)

        for target in targets:
            self._add_code2allele(target[col2idx['ID']], 0, target[col2idx['REF']])
            for code, allele in enumerate(target[col2idx['ALT']].split(','), 1):
                self._add_code2allele(target[col2idx['ID']], code, allele)

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

def __group_spandel(deletion: list[str], targets:list[list[str]], col2idx:dict) -> str:

    at = AlleleTranslator()

    deletion = set_id(deletion, col2idx)
    targets = [set_id(target, col2idx) for target in targets]

    at.update_code2allele(deletion, targets, col2idx)

    expanded_deletion = _expand_spandel(deletion, targets, col2idx)
    dpos_expanded = int(expanded_deletion[col2idx['POS']])
    dref_expanded = expanded_deletion[col2idx['REF']]
    dalts_expanded = expanded_deletion[col2idx['ALT']].split(',')

    at.update_allele2allele(
        deletion[col2idx['ID']], 
        deletion[col2idx['REF']], 
        dref_expanded,
    )

    for old_allele, new_allele in zip(deletion[col2idx['ALT']].split(','), dalts_expanded):
        at.update_allele2allele(
            deletion[col2idx['ID']], 
            old_allele,
            new_allele,
        )

    new_alts = dalts_expanded

    for target in targets:
        tpos = int(target[col2idx['POS']])
        tref = target[col2idx['REF']]
        talts = target[col2idx['ALT']].split(',')
        prefix = dref_expanded[0:(tpos - dpos_expanded + 1 - 1)]
        suffix = dref_expanded[tpos - dpos_expanded + 1 -1 + len(tref) + 1 - 1:]
        at.update_allele2allele(target[col2idx['ID']], tref, dref_expanded)
        for talt in talts:
            if talt == '*':
                continue
            else:
                new_talt = prefix + talt + suffix
                at.update_allele2allele(target[col2idx['ID']], talt, new_talt)
            new_alts.append(new_talt)

    new_alts = sorted(new_alts)

    at.update_allele2code(dref_expanded, 0)
    for code, new_alt in enumerate(new_alts, 1):
        at.update_allele2code(new_alt, code)

    expanded_deletion[col2idx['ALT']] = ','.join(new_alts)
    expanded_deletion[col2idx['ID']] = deletion[col2idx['ID']]

    # print(at.code2allele)
        # self.allele2allele = {}
    # print(at.allele2code)

    grouped_record = update_calls(expanded_deletion, targets, col2idx, at)

    return grouped_record


def update_calls(expanded_deletion: list, targets: list[list], col2idx, allele_translator)-> str:

    idx = 0

    result = list()

    for idx in range(0, len(expanded_deletion)):
        if idx < 10:
            result.append(expanded_deletion[idx])
            continue

        for gt in expanded_deletion[idx].split('\t'):

            alleles = list()

            for code in gt.split('/'):
                if code == '.':
                    continue
                code = int(code)
                allele_old = allele_translator.code2allele(expanded_deletion[col2idx['ID']], code)
                allele_new = allele_translator.allele2allele(expanded_deletion[col2idx['ID']], allele_old)
                alleles.append(allele_new)

            for target in targets:
                for code in target[idx].split('/'):
                    if code == '.':
                        continue
                    code = int(code)
                    allele_old = allele_translator.code2allele(target[col2idx['ID']], code)
                    if allele_old == '*':
                        continue
                    allele_new = allele_translator.allele2allele(target[col2idx['ID']], allele_old)
                    if allele_new not in alleles:
                        alleles.append(allele_new)
            if len(alleles) > 2:
                print('##', expanded_deletion[col2idx['ID']], alleles)
                call = './.'
            elif len(alleles) == 0:
                call = './.'
            else:
                codes = list()
                for allele in alleles:
                    code = allele_translator.allele2code(allele)
                    codes.append(code)
                call = '\t'.join(codes)
            result.append(call)
    return '\t'.join(result)


        











        





