#!/usr/bin/env python
from pathlib import Path
import pytest
from genomics.spandel import ( 
        group_spandel, 
        _expand_spandel, 
        update_calls, 
        AlleleTranslator, 
        _new_allele_translator, 
        _create_backbone, 
        SpanFamily, 
        Site,
)
import polars as pl
from genomics.vcf import Vcf

def test__expand_spandel():
    parent = ['chr22', '100','rs100','ATT','A','.','.','.']
    children = list()
    children.append(['chr22', '101','.','TTG','*,T','.','.','.'])

    col2idx = {'#CHROM':0, 'POS': 1, 'ID': 2, 'REF': 3, 'ALT': 4, 'QUAL': 5, 'FILTER': 6, 'INFO':7}

    record = _expand_spandel(parent, children, col2idx)

    assert list == type(record)
    assert 'chr22' == record[0]
    assert '100' == record[1]
    assert 'rs100' == record[2]
    assert 'ATTG' == record[3]
    assert 'AG' == record[4]

    parent = ['chr22', '100','rs100','AGATGAAAT','A','.','.','.']
    children = list()
    children.append(['chr22', '101','.','GATGAAATGATGAA','*','.','.','.'])
    children.append(['chr22', '102','.','A','ATGAAG,*','.','.','.'])

    col2idx = {'#CHROM':0, 'POS': 1, 'ID': 2, 'REF': 3, 'ALT': 4, 'QUAL': 5, 'FILTER': 6, 'INFO':7}

    record = _expand_spandel(parent, children, col2idx)
    
    assert list == type(record)
    assert 'chr22' == record[0]
    assert '100' == record[1]
    assert 'rs100' == record[2]
    assert 'AGATGAAATGATGAA' == record[3]
    assert 'AGATGAA' == record[4]



@pytest.mark.skip()
def test_SpanFamily_combine(tmp_path):

    parent = ['chr22', '10516150','.','GTA','G','.','.','.']
    children = list()
    children.append(['chr22', '10516151','.','T','*','.','.','.'])
    children.append(['chr22', '10516152','.','A','*,G','.','.','.'])

    col2idx = {'#CHROM':0, 'POS': 1, 'ID': 2, 'REF': 3, 'ALT': 4, 'QUAL': 5, 'FILTER': 6, 'INFO':7}

    record = SpanFamily(parent, children, col2idx).combine()

    assert str == type(record)
    assert 'chr22' == record.split('\t')[0]
    assert '10516150' == record.split('\t')[1]
    assert 'chr22:10516150:GTA:G' == record.split('\t')[2]
    assert 'GTA' == record.split('\t')[3]
    assert 'G,GTG' == record.split('\t')[4]

    parent = ['chr22', '100','.','ATT','A','.','.','.']
    children = list()
    children.append(['chr22', '101','.','TTG','*,T','.','.','.'])

    col2idx = {'#CHROM':0, 'POS': 1, 'ID': 2, 'REF': 3, 'ALT': 4, 'QUAL': 5, 'FILTER': 6, 'INFO':7}

    record = SpanFamily(parent, children, col2idx).combine()

    assert str == type(record)
    assert 'chr22' == record.split('\t')[0]
    assert '100' == record.split('\t')[1]
    assert 'chr22:100:ATT:A' == record.split('\t')[2]
    assert 'ATTG' == record.split('\t')[3]
    assert 'AG,AT' == record.split('\t')[4]

def test_update_calls():

    # 01234
    # GCCCCCACCC
    # G     
    #       A
    #       T

    col2idx = {'#CHROM': 0, 'POS': 1, 'ID': 2, 'REF': 3, 'ALT': 4}
    parent = ['chr1', 100, 'rs100', 'GCCCCACCC', 'G', '.', '.', '.', 'GT', '0/0']
    child1 = ['chr1', 105, 'rs101', 'A', 'T,*', '.', '.', '.', 'GT', '0/1']
    children = [child1]
    parent_expanded = ['chr1', 100, 'rs100', 'GCCCCACCC', 'G', '.', '.', '.', 'GT', '0/0']

    at = _new_allele_translator(parent, parent_expanded, children, col2idx)

    record = update_calls(parent, children, col2idx, at)

    assert str == type(record)
    assert 'chr1' == record.split('\t')[0]
    assert '100' == record.split('\t')[1]
    assert 'rs100' == record.split('\t')[2]
    assert 'GCCCCACCC' == record.split('\t')[3]
    assert 'G,GCCCCTCCC' == record.split('\t')[4]
    assert '0/2' == record.split('\t')[9]

    col2idx = {'#CHROM': 0, 'POS': 1, 'ID': 2, 'REF': 3, 'ALT': 4}
    parent = ['chr1', 100, 'rs100', 'GCCCCACCC', 'G', '.', '.', '.', 'GT', '0/1']
    child1 = ['chr1', 105, 'rs101', 'A', 'T,*', '.', '.', '.', 'GT', '1/2']
    children = [child1]
    parent_expanded = ['chr1', 100, 'rs100', 'GCCCCACCC', 'G', '.', '.', '.', 'GT', '0/1']

    at = _new_allele_translator(parent, parent_expanded, children, col2idx)

    record = update_calls(parent, children, col2idx, at)

    assert str == type(record)
    assert 'chr1' == record.split('\t')[0]
    assert '100' == record.split('\t')[1]
    assert 'rs100' == record.split('\t')[2]
    assert 'GCCCCACCC' == record.split('\t')[3]
    assert 'G,GCCCCTCCC' == record.split('\t')[4]
    assert '1/2' == record.split('\t')[9]

    col2idx = {'#CHROM': 0, 'POS': 1, 'ID': 2, 'REF': 3, 'ALT': 4}
    parent = ['chr1', 100, 'rs100', 'GCCCCACCC', 'G', '.', '.', '.', 'GT', '1/1']
    child1 = ['chr1', 105, 'rs101', 'A', 'T,*', '.', '.', '.', 'GT', '2/2']
    children = [child1]
    parent_expanded = ['chr1', 100, 'rs100', 'GCCCCACCC', 'G', '.', '.', '.', 'GT', '0/1']

    at = _new_allele_translator(parent, parent_expanded, children, col2idx)

    record = update_calls(parent, children, col2idx, at)

    assert str == type(record)
    assert 'chr1' == record.split('\t')[0]
    assert '100' == record.split('\t')[1]
    assert 'rs100' == record.split('\t')[2]
    assert 'GCCCCACCC' == record.split('\t')[3]
    assert 'G,GCCCCTCCC' == record.split('\t')[4]
    assert '1/1' == record.split('\t')[9]

    col2idx = {'#CHROM': 0, 'POS': 1, 'ID': 2, 'REF': 3, 'ALT': 4}
    parent = ['chr1', 100, 'rs100', 'GCCCCACCC', 'G', '.', '.', '.', 'GT', '0/1']
    child1 = ['chr1', 105, 'rs101', 'A', 'T,*', '.', '.', '.', 'GT', '0/2']
    children = [child1]
    parent_expanded = ['chr1', 100, 'rs100', 'GCCCCACCC', 'G', '.', '.', '.', 'GT', '0/1']

    at = _new_allele_translator(parent, parent_expanded, children, col2idx)

    record = update_calls(parent, children, col2idx, at)

    assert str == type(record)
    assert 'chr1' == record.split('\t')[0]
    assert '100' == record.split('\t')[1]
    assert 'rs100' == record.split('\t')[2]
    assert 'GCCCCACCC' == record.split('\t')[3]
    assert 'G,GCCCCTCCC' == record.split('\t')[4]
    assert '0/1' == record.split('\t')[9]

    # # 012
    # #
    # # CAA
    # # C
    # # CAAAA x
    # # CAAA  x
    # #
    # #  A    x
    # #  C
    # #
    # #   A   x
    # #   C
    #
    # col2idx = {'#CHROM': 0, 'POS': 1, 'ID': 2, 'REF': 3, 'ALT': 4}
    # parent = ['chr1', 100, 'rs100', 'CAA', 'C,CAAAA,CAAA', '.', '.', '.', 'GT', '2/3']
    # child1 = ['chr1', 101, 'rs101', 'A', '*,C', '.', '.', '.', 'GT', '0/0']
    # child2 = ['chr1', 102, 'rs102', 'A', '*,C', '.', '.', '.', 'GT', '0/0']
    # children = [child1, child2]
    # parent_expanded = ['chr1', 100, 'rs100', 'CAA', 'C,CAAAA,CAAA', '.', '.', '.', 'GT', '2/3']
    #
    # at = _new_allele_translator(parent, parent_expanded, children, col2idx)
    #
    # record = update_calls(parent, children, col2idx, at)
    #
    # assert str == type(record)
    # assert 'chr1' == record.split('\t')[0]
    # assert '100' == record.split('\t')[1]
    # assert 'rs100' == record.split('\t')[2]
    # assert 'CAA' == record.split('\t')[3]
    # assert 'C,CAAA,CAAAA,CAC,CCA' == record.split('\t')[4]
    # assert '2/3' == record.split('\t')[9]

    col2idx = {'#CHROM': 0, 'POS': 1, 'ID': 2, 'REF': 3, 'ALT': 4}
    parent = ['chr1', 100, 'rs100', 'CAA', 'C,CAAAA,CAAA', '.', '.', '.', 'GT', '1/2']
    child1 = ['chr1', 101, 'rs101', 'A', '*,C', '.', '.', '.', 'GT', '1/2']
    child2 = ['chr1', 102, 'rs102', 'A', '*,C', '.', '.', '.', 'GT', '1/2']
    children = [child1, child2]
    parent_expanded = ['chr1', 100, 'rs100', 'CAA', 'C,CAAAA,CAAA', '.', '.', '.', 'GT', '1/2']

    at = _new_allele_translator(parent, parent_expanded, children, col2idx)

    record = update_calls(parent, children, col2idx, at)

    assert str == type(record)
    assert 'chr1' == record.split('\t')[0]
    assert '100' == record.split('\t')[1]
    assert 'rs100' == record.split('\t')[2]
    assert 'CAA' == record.split('\t')[3]
    assert 'C,CAAA,CAAAA,CAC,CCA' == record.split('\t')[4]
    assert './.' == record.split('\t')[9]

    col2idx = {'#CHROM': 0, 'POS': 1, 'ID': 2, 'REF': 3, 'ALT': 4}
    parent = ['chr1', 100, 'rs100', 'CAA', 'C,CAAAA,CAAA', '.', '.', '.', 'GT', '0/0']
    child1 = ['chr1', 101, 'rs101', 'A', '*,C', '.', '.', '.', 'GT', '0/0']
    child2 = ['chr1', 102, 'rs102', 'A', '*,C', '.', '.', '.', 'GT', '0/0']
    children = [child1, child2]
    parent_expanded = ['chr1', 100, 'rs100', 'CAA', 'C,CAAAA,CAAA', '.', '.', '.', 'GT', '0/0']

    at = _new_allele_translator(parent, parent_expanded, children, col2idx)

    record = update_calls(parent, children, col2idx, at)

    assert str == type(record)
    assert 'chr1' == record.split('\t')[0]
    assert '100' == record.split('\t')[1]
    assert 'rs100' == record.split('\t')[2]
    assert 'CAA' == record.split('\t')[3]
    assert 'C,CAAA,CAAAA,CAC,CCA' == record.split('\t')[4]
    assert '0/0' == record.split('\t')[9]

    col2idx = {'#CHROM': 0, 'POS': 1, 'ID': 2, 'REF': 3, 'ALT': 4}
    parent = ['chrY', 100, 'rs100', 'CAA', 'C,CAAAA,CAAA', '.', '.', '.', 'GT', '1']
    child1 = ['chr1', 101, 'rs101', 'A', '*,C', '.', '.', '.', 'GT', '1']
    child2 = ['chr1', 102, 'rs102', 'A', '*,C', '.', '.', '.', 'GT', '1']
    children = [child1, child2]
    parent_expanded = ['chr1', 100, 'rs100', 'CAA', 'C,CAAAA,CAAA', '.', '.', '.', 'GT', '1']

    at = _new_allele_translator(parent, parent_expanded, children, col2idx)

    record = update_calls(parent, children, col2idx, at)

    assert str == type(record)
    assert 'chrY' == record.split('\t')[0]
    assert '100' == record.split('\t')[1]
    assert 'rs100' == record.split('\t')[2]
    assert 'CAA' == record.split('\t')[3]
    assert 'C,CAAA,CAAAA,CAC,CCA' == record.split('\t')[4]
    assert '1' == record.split('\t')[9]

# def test_group_spandel(tmp_path):
#
#     vcf_file = Path(__file__).parents[0] / 'fixture.vcf'
#     output_vcf_file = tmp_path / 'output.vcf'
#
#     # print(output_vcf_file)
#     output_vcf_file = group_spandel(vcf_file, output_vcf_file)
#     # print(Vcf(output_vcf_file, tmp_path / 'tmp').bgzip().to_df(site_only = True).select(pl.col(['chrom', 'pos', 'id', 'ref', 'alt'])))
#     #

def test__create_backbone():
    col2idx = {'#CHROM':0, 'POS': 1, 'ID': 2, 'REF': 3, 'ALT': 4, 'QUAL': 5, 'FILTER': 6, 'INFO':7}

    parent = ['chr22', '10516150','.','GTA','G','.','.','.']
    children = [ 
        ['chr22', '10516151','.','T','*','.','.','.'], 
        ['chr22', '10516152','.','A','*,G','.','.','.'],
    ]

    backbone = _create_backbone(parent, children, col2idx)
    assert 10516150 == backbone['pos']
    assert 10516152 == backbone['end']
    assert 'GTA' == backbone['seq']

    parent = ['chr22', '10516150','.','GTA','G','.','.','.']
    children = [ 
        ['chr22', '10516151','.','T','*','.','.','.'], 
        ['chr22', '10516152','.','AC','*,A','.','.','.'],
    ]

    backbone = _create_backbone(parent, children, col2idx)
    assert 10516150 == backbone['pos']
    assert 10516153 == backbone['end']
    assert 'GTAC' == backbone['seq']


def test_SpanFamily_expand():
    parent = ['chr22', '10516150','.','GTA','G','.','.','.']
    children = list()
    children.append(['chr22', '10516151','.','T','*','.','.','.'])
    children.append(['chr22', '10516152','.','A','*,G','.','.','.'])

    col2idx = {'#CHROM':0, 'POS': 1, 'ID': 2, 'REF': 3, 'ALT': 4, 'QUAL': 5, 'FILTER': 6, 'INFO':7}

    sf = SpanFamily(parent, children, col2idx)

    assert 0 == len(sf._expanded_sites)

    sf.expand()

    assert 3 == len(sf._expanded_sites)
    assert 'chr22' == sf._expanded_sites[0].chrom
    assert 10516150 == sf._expanded_sites[0].pos
    assert 'GTA' == sf._expanded_sites[0].ref
    assert 'G' == sf._expanded_sites[0].alt

    assert 'chr22' == sf._expanded_sites[1].chrom
    assert 10516150 == sf._expanded_sites[1].pos
    assert 'GTA' == sf._expanded_sites[1].ref
    assert '*' == sf._expanded_sites[1].alt

    assert 'chr22' == sf._expanded_sites[2].chrom
    assert 10516150 == sf._expanded_sites[2].pos
    assert 'GTA' == sf._expanded_sites[2].ref
    assert '*,GTG' == sf._expanded_sites[2].alt

def test_Site_expand():

    backbone = {'pos': 10516150, 'end': 10516152, 'seq': 'GTA'}

    record = Site('chr22', 10516150, 'GTA', 'G').expand(backbone)
    assert 'chr22' == record.chrom
    assert 10516150 == record.pos
    assert 'GTA' == record.ref
    assert 'G' == record.alt

    record = Site('chr22', 10516151, 'T', '*').expand(backbone)
    assert 'chr22' == record.chrom
    assert 10516150 == record.pos
    assert 'GTA' == record.ref
    assert '*' == record.alt

    record = Site('chr22', 10516152, 'A', 'G').expand(backbone)
    assert 'chr22' == record.chrom
    assert 10516150 == record.pos
    assert 'GTA' == record.ref
    assert 'GTG' == record.alt


    record = Site('chr1', 49515, 'G', '*').expand({'pos': 49514, 'end': 49515, 'seq': 'AG'})

    assert 'chr1' == record.chrom
    assert 49514 == record.pos
    assert 'AG' == record.ref
    assert '*' == record.alt

def test_Site_translate():

    backbone = {'pos': 10516150, 'end': 10516152, 'seq': 'GTA'}

    record = Site('chr22', 10516150, 'GTA', 'G').expand(backbone)
    assert 'GTA' == record.translate('GTA')
    assert 'G' == record.translate('G')

    record = Site('chr22', 10516151, 'T', '*').expand(backbone)
    assert 'GTA' == record.translate('T')
    assert '*' == record.translate('*')

    record = Site('chr22', 10516152, 'A', 'G').expand(backbone)
    assert 'GTA' == record.translate('A')
    assert 'GTG' == record.translate('G')

    record = Site('chr1', 49515, 'G', '*').expand({'pos': 49514, 'end': 49515, 'seq': 'AG'})

    assert 'AG' == record.translate('G')
    assert '*' == record.translate('*')

def test_Site_allele():

    assert 'GTA' == Site('chr22', 10516150, 'GTA', 'A,G').allele('0')
    assert 'A' == Site('chr22', 10516150, 'GTA', 'A,G').allele('1')
    assert 'G' == Site('chr22', 10516150, 'GTA', 'A,G').allele('2')
    assert None == Site('chr22', 10516150, 'GTA', 'A,G').allele('.')
    assert '*' == Site('chr22', 10516150, 'GTA', '*,G').allele('1')

