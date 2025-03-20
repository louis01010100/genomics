#!/usr/bin/env python
from pathlib import Path
from genomics.spandel import ( group_spandel, __group_spandel, _expand_spandel, update_calls)
import polars as pl
from genomics.vcf import Vcf

# def test__expand_spandel():
#     deletion = ['chr22', '100','.','ATT','A','.','.','.']
#     targets = list()
#     targets.append(['chr22', '101','.','TTG','*,T','.','.','.'])
#
#     col2idx = {'#CHROM':0, 'POS': 1, 'ID': 2, 'REF': 3, 'ALT': 4, 'QUAL': 5, 'FILTER': 6, 'INFO':7}
#
#     record = _expand_spandel(deletion, targets, col2idx)
#
#     assert list == type(record)
#     assert 'chr22' == record[0]
#     assert '100' == record[1]
#     assert '.' == record[2]
#     assert 'ATTG' == record[3]
#     assert 'AG' == record[4]
#
#     deletion = ['chr22', '100','.','AGATGAAAT','A','.','.','.']
#     targets = list()
#     targets.append(['chr22', '101','.','GATGAAATGATGAA','*','.','.','.'])
#     targets.append(['chr22', '102','.','A','ATGAAG,*','.','.','.'])
#
#     col2idx = {'#CHROM':0, 'POS': 1, 'ID': 2, 'REF': 3, 'ALT': 4, 'QUAL': 5, 'FILTER': 6, 'INFO':7}
#
#     record = _expand_spandel(deletion, targets, col2idx)
#
#     assert list == type(record)
#     assert 'chr22' == record[0]
#     assert '100' == record[1]
#     assert '.' == record[2]
#     assert 'AGATGAAATGATGAA' == record[3]
#     assert 'AGATGAA' == record[4]
#
#
#
# def test___group_spandel(tmp_path):
#
#     deletion = ['chr22', '10516150','.','GTA','G','.','.','.']
#     targets = list()
#     targets.append(['chr22', '10516151','.','T','*','.','.','.'])
#     targets.append(['chr22', '10516152','.','A','*,G','.','.','.'])
#
#     col2idx = {'#CHROM':0, 'POS': 1, 'ID': 2, 'REF': 3, 'ALT': 4, 'QUAL': 5, 'FILTER': 6, 'INFO':7}
#
#     record = __group_spandel(deletion, targets, col2idx)
#
#     assert str == type(record)
#     assert 'chr22' == record.split('\t')[0]
#     assert '10516150' == record.split('\t')[1]
#     assert 'chr22:10516150:GTA:G' == record.split('\t')[2]
#     assert 'GTA' == record.split('\t')[3]
#     assert 'G,GTG' == record.split('\t')[4]
#
#     deletion = ['chr22', '100','.','ATT','A','.','.','.']
#     targets = list()
#     targets.append(['chr22', '101','.','TTG','*,T','.','.','.'])
#
#     col2idx = {'#CHROM':0, 'POS': 1, 'ID': 2, 'REF': 3, 'ALT': 4, 'QUAL': 5, 'FILTER': 6, 'INFO':7}
#
#     record = __group_spandel(deletion, targets, col2idx)
#
#     assert str == type(record)
#     assert 'chr22' == record.split('\t')[0]
#     assert '100' == record.split('\t')[1]
#     assert 'chr22:100:ATT:A' == record.split('\t')[2]
#     assert 'ATTG' == record.split('\t')[3]
#     assert 'AG,AT' == record.split('\t')[4]
#
# def test_update_calls():
#     col2idx = ['#CHROM': 0, 'POS': 1, 'REF': 2, 'ALT': 3,]
#     deletion = ['chr22', 100, 'A', 'AAT,AATAT', '.', '.', '.', 'GT', '0/2']
#     target1 = ['chr22', 101, 'AT', '*', '.', '.', '.', 'GT', '0/1']
#
#     expanded_deletion = ['chr22', 100, 'AT', 'AATT,AATATT', '.', '.', '.', 'GT', '0/2']
#
#     at.update_allele2allele(
#     at = AlleleTranslator()
#     at.update_code2allele(deletion, [target1], col2idx)
#
#     update_calls(deletion: list, targets: list[list], col2idx, allele_translator)-> str:

def test_group_spandel(tmp_path):

    vcf_file = Path(__file__).parents[0] / 'fixture.vcf'
    output_vcf_file = tmp_path / 'output.vcf'

    print(output_vcf_file)
    output_vcf_file = group_spandel(vcf_file, output_vcf_file)
    # print(Vcf(output_vcf_file, tmp_path / 'tmp').bgzip().to_df(site_only = True).select(pl.col(['chrom', 'pos', 'id', 'ref', 'alt'])))
    #
