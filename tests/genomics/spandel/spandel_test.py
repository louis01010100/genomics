#!/usr/bin/env python
from pathlib import Path
from genomics.spandel import ( group_spandel, __group_spandel, _expand_spandel, update_calls, AlleleTranslator, _new_allele_translator)
import polars as pl
from genomics.vcf import Vcf

def test__expand_spandel():
    deletion = ['chr22', '100','.','ATT','A','.','.','.']
    targets = list()
    targets.append(['chr22', '101','.','TTG','*,T','.','.','.'])

    col2idx = {'#CHROM':0, 'POS': 1, 'ID': 2, 'REF': 3, 'ALT': 4, 'QUAL': 5, 'FILTER': 6, 'INFO':7}

    record = _expand_spandel(deletion, targets, col2idx)

    assert list == type(record)
    assert 'chr22' == record[0]
    assert '100' == record[1]
    assert '.' == record[2]
    assert 'ATTG' == record[3]
    assert 'AG' == record[4]

    deletion = ['chr22', '100','.','AGATGAAAT','A','.','.','.']
    targets = list()
    targets.append(['chr22', '101','.','GATGAAATGATGAA','*','.','.','.'])
    targets.append(['chr22', '102','.','A','ATGAAG,*','.','.','.'])

    col2idx = {'#CHROM':0, 'POS': 1, 'ID': 2, 'REF': 3, 'ALT': 4, 'QUAL': 5, 'FILTER': 6, 'INFO':7}

    record = _expand_spandel(deletion, targets, col2idx)
    
    assert list == type(record)
    assert 'chr22' == record[0]
    assert '100' == record[1]
    assert '.' == record[2]
    assert 'AGATGAAATGATGAA' == record[3]
    assert 'AGATGAA' == record[4]



def test___group_spandel(tmp_path):

    deletion = ['chr22', '10516150','.','GTA','G','.','.','.']
    targets = list()
    targets.append(['chr22', '10516151','.','T','*','.','.','.'])
    targets.append(['chr22', '10516152','.','A','*,G','.','.','.'])

    col2idx = {'#CHROM':0, 'POS': 1, 'ID': 2, 'REF': 3, 'ALT': 4, 'QUAL': 5, 'FILTER': 6, 'INFO':7}

    record = __group_spandel(deletion, targets, col2idx)

    assert str == type(record)
    assert 'chr22' == record.split('\t')[0]
    assert '10516150' == record.split('\t')[1]
    assert 'chr22:10516150:GTA:G' == record.split('\t')[2]
    assert 'GTA' == record.split('\t')[3]
    assert 'G,GTG' == record.split('\t')[4]

    deletion = ['chr22', '100','.','ATT','A','.','.','.']
    targets = list()
    targets.append(['chr22', '101','.','TTG','*,T','.','.','.'])

    col2idx = {'#CHROM':0, 'POS': 1, 'ID': 2, 'REF': 3, 'ALT': 4, 'QUAL': 5, 'FILTER': 6, 'INFO':7}

    record = __group_spandel(deletion, targets, col2idx)

    assert str == type(record)
    assert 'chr22' == record.split('\t')[0]
    assert '100' == record.split('\t')[1]
    assert 'chr22:100:ATT:A' == record.split('\t')[2]
    assert 'ATTG' == record.split('\t')[3]
    assert 'AG,AT' == record.split('\t')[4]

# def test_update_calls():
#
#     # 01234
#     # AAAAG   0
#     # A       1
#     # AAAAG   0
#     # AA      1
#
#     col2idx = {'#CHROM': 0, 'POS': 1, 'ID': 2, 'REF': 3, 'ALT': 4}
#     deletion = ['chr1', 100, 'rs100', 'AAAAG', 'A', '.', '.', '.', 'GT', '0/1']
#     target1 = ['chr1', 101, 'rs101', 'AAAG', 'A,*', '.', '.', '.', 'GT', '0/1']
#     targets = [target1]
#     expanded_deletion = ['chr1', 100, 'rs100', 'AAAAG', 'A', '.', '.', '.', 'GT', './.']
#
#     at = _new_allele_translator(deletion, expanded_deletion, targets, col2idx)
#
#     record = update_calls(deletion, targets, col2idx, at)
#
#     assert str == type(record)
#     assert 'chr1' == record.split('\t')[0]
#     assert '100' == record.split('\t')[1]
#     assert 'rs100' == record.split('\t')[2]
#     assert 'AAAAG' == record.split('\t')[3]
#     assert 'A,AA' == record.split('\t')[4]
#     assert './.' == record.split('\t')[9]
#
#
#     # 012
#     #
#     # CAA
#     # C
#     # CAAAA x
#     # CAAA  x
#     #
#     #  A    x
#     #  C
#     #
#     #   A   x
#     #   C
#
#     col2idx = {'#CHROM': 0, 'POS': 1, 'ID': 2, 'REF': 3, 'ALT': 4}
#     deletion = ['chr1', 100, 'rs100', 'CAA', 'C,CAAAA,CAAA', '.', '.', '.', 'GT', '2/3']
#     target1 = ['chr1', 101, 'rs101', 'A', '*,C', '.', '.', '.', 'GT', '0/0']
#     target2 = ['chr1', 102, 'rs102', 'A', '*,C', '.', '.', '.', 'GT', '0/0']
#     targets = [target1, target2]
#     expanded_deletion = ['chr1', 100, 'rs100', 'CAA', 'C,CAAAA,CAAA', '.', '.', '.', 'GT', '2/3']
#
#     at = _new_allele_translator(deletion, expanded_deletion, targets, col2idx)
#
#     record = update_calls(deletion, targets, col2idx, at)
#
#     assert str == type(record)
#     assert 'chr1' == record.split('\t')[0]
#     assert '100' == record.split('\t')[1]
#     assert 'rs100' == record.split('\t')[2]
#     assert 'CAA' == record.split('\t')[3]
#     assert 'C,CAAA,CAAAA,CAC,CCA' == record.split('\t')[4]
#     assert '2/3' == record.split('\t')[9]
#
#     col2idx = {'#CHROM': 0, 'POS': 1, 'ID': 2, 'REF': 3, 'ALT': 4}
#     deletion = ['chr1', 100, 'rs100', 'CAA', 'C,CAAAA,CAAA', '.', '.', '.', 'GT', '1/2']
#     target1 = ['chr1', 101, 'rs101', 'A', '*,C', '.', '.', '.', 'GT', '1/2']
#     target2 = ['chr1', 102, 'rs102', 'A', '*,C', '.', '.', '.', 'GT', '1/2']
#     targets = [target1, target2]
#     expanded_deletion = ['chr1', 100, 'rs100', 'CAA', 'C,CAAAA,CAAA', '.', '.', '.', 'GT', '1/2']
#
#     at = _new_allele_translator(deletion, expanded_deletion, targets, col2idx)
#
#     record = update_calls(deletion, targets, col2idx, at)
#
#     assert str == type(record)
#     assert 'chr1' == record.split('\t')[0]
#     assert '100' == record.split('\t')[1]
#     assert 'rs100' == record.split('\t')[2]
#     assert 'CAA' == record.split('\t')[3]
#     assert 'C,CAAA,CAAAA,CAC,CCA' == record.split('\t')[4]
#     assert './.' == record.split('\t')[9]



def test_group_spandel(tmp_path):

    vcf_file = Path(__file__).parents[0] / 'fixture.vcf'
    output_vcf_file = tmp_path / 'output.vcf'

    print(output_vcf_file)
    output_vcf_file = group_spandel(vcf_file, output_vcf_file)
    # print(Vcf(output_vcf_file, tmp_path / 'tmp').bgzip().to_df(site_only = True).select(pl.col(['chrom', 'pos', 'id', 'ref', 'alt'])))
    #
