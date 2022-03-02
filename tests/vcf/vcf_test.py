#!/usr/bin/env python
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
