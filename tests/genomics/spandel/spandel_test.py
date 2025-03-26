#!/usr/bin/env python
from pathlib import Path
from genomics.spandel import ( group_spandel, _expand_spandel, update_calls, AlleleTranslator, _new_allele_translator, _expand, _create_backbone, SpanFamily)
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

    pos, end, ref = _create_backbone(parent, children, col2idx)
    assert 10516150 == pos
    assert 10516152 == end
    assert 'GTA' == ref

    parent = ['chr22', '10516150','.','GTA','G','.','.','.']
    children = [ 
        ['chr22', '10516151','.','T','*','.','.','.'], 
        ['chr22', '10516152','.','AC','*,A','.','.','.'],
    ]

    pos, end, ref = _create_backbone(parent, children, col2idx)
    assert 10516150 == pos
    assert 10516153 == end
    assert 'GTAC' == ref


def test__expand():
    col2idx = {'#CHROM':0, 'POS': 1, 'ID': 2, 'REF': 3, 'ALT': 4, 'QUAL': 5, 'FILTER': 6, 'INFO':7}

    snv = ['chr22', '10516150','.','GTA','G','.','.','.']

    record, map_ = _expand(snv, col2idx, 10516150, 10516152, 'GTA' )

    assert list == type(record)
    assert 'chr22' == record[0]
    assert 10516150 == record[1]
    assert 'chr22:10516150:GTA:G'[2]
    assert 'GTA' == record[3]
    assert 'G' == record[4]

    snv = ['chr22', '10516151','.','T','*','.','.','.']
    record, map_ = _expand(snv, col2idx, 10516150, 10516152, 'GTA' )

    assert list == type(record)
    assert 'chr22' == record[0]
    assert 10516150 == record[1]
    assert 'chr22:10516151:T:*'[2]
    assert 'GTA' == record[3]
    assert '*' == record[4]

    snv = ['chr22', '10516152','.','A','G','.','.','.']
    record, map_ = _expand(snv, col2idx, 10516150, 10516152, 'GTA' )

    assert list == type(record)
    assert 'chr22' == record[0]
    assert 10516150 == record[1]
    assert 'chr22:10516152:A:*'[2]
    assert 'GTA' == record[3]
    assert 'GTG' == record[4]


    snv = ['chr1', '49515','.','G','*','.','.','.']
    record, map_ = _expand(snv, col2idx, 49514, 49515, 'AG' )

    assert list == type(record)
    assert 'chr1' == record[0]
    assert 49514 == record[1]
    assert 'chr1::49515:G:*'[2]
    assert 'AG' == record[3]
    assert '*' == record[4]
