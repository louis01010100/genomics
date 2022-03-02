#!/usr/bin/env python
import pandas as pd
from vcf import Vcf
from vcf import bgzip
from vcf import index
from pathlib import Path
import shutil


def test_bgzip(tmp_path):
    fixture_dir = Path(__file__).parents[0] / 'fixture'
    output_filepath = bgzip(fixture_dir / 'sample.vcf', tmp_path)
    assert 'sample.vcf.bgz' == output_filepath.name
    output_filepath = bgzip(fixture_dir / 'sample.vcf.gz', tmp_path)
    assert 'sample.vcf.bgz' == output_filepath.name


def test_index(tmp_path):
    src_data_file = Path(__file__).parents[0] / 'fixture/sample.vcf.bgz'
    dst_data_file = tmp_path / src_data_file.name

    shutil.copy2(src_data_file, dst_data_file)

    index(dst_data_file, tmp_path)
    index_file = dst_data_file.with_suffix('.bgz.csi')
    assert index_file.exists()


def test_to_tsv(tmp_path):
    vcf_file = Path(__file__).parents[0] / 'fixture/sample.vcf.bgz'
    vcf = Vcf(vcf_file, tmp_path)
    output_file = vcf.to_tsv(
        '[%CHROM\t%POS\t%ID\t%REF\t%ALT\t%SAMPLE\t%GT\n]')

    data = pd.read_csv(output_file, header=0, sep='\t', dtype='str')

    record1 = data.iloc[0, :]

    assert record1['CHROM'] == 'chr21'
    assert record1['POS'] == '5030578'
    assert record1['ID'] == '.'
    assert record1['REF'] == 'C'
    assert record1['ALT'] == 'T'
    assert record1['SAMPLE'] == 'HG00403'
    assert record1['GT'] == '0|0'
