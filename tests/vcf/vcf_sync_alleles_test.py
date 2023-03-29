
import filecmp
from pathlib import Path

from vcf.vcf import Vcf, sync_alleles


def test_sync_alleles_forward(tmp_path: Path):

    fixture_dir = Path(__file__).parents[0] / 'fixture_sync_alleles' 
    vcf_file_1 = fixture_dir / 'one.vcf'
    vcf_file_2 = fixture_dir / 'two.vcf'

    tmp_dir_before = tmp_path / 'before'
    tmp_dir_after = tmp_path / 'after'

    vcf_file_1 = Vcf(vcf_file_1, tmp_dir_before).bgzip()
    vcf_file_2 = Vcf(vcf_file_2, tmp_dir_before).bgzip()

    vcf_file_1_df = vcf_file_1.to_df()
    vcf_file_2_df = vcf_file_2.to_df()

    assert {'A'} == set(vcf_file_1_df['ref'])
    assert {'AGGAGTC'} == set(vcf_file_2_df['ref'])

    sync_alleles(vcf_file_1.filepath, vcf_file_2.filepath, tmp_path / 'after')

    vcf1_observed = Vcf(tmp_dir_after / 'one-sync.vcf.bgz', tmp_dir_after).to_df()
    vcf2_observed = Vcf(tmp_dir_after / 'two-sync.vcf.bgz', tmp_dir_after).to_df()

    assert {'AGGAGTC'} == set(vcf1_observed['ref'])
    assert {'AGGAGTC'} == set(vcf2_observed['ref'])


def test_sync_alleles_reverse(tmp_path: Path):

    fixture_dir = Path(__file__).parents[0] / 'fixture_sync_alleles' 
    vcf_file_1 = fixture_dir / 'one.vcf'
    vcf_file_2 = fixture_dir / 'two.vcf'

    tmp_dir_before = tmp_path / 'before'
    tmp_dir_after = tmp_path / 'after'

    vcf_file_1 = Vcf(vcf_file_1, tmp_dir_before).bgzip()
    vcf_file_2 = Vcf(vcf_file_2, tmp_dir_before).bgzip()

    vcf_file_1_df = vcf_file_1.to_df()
    vcf_file_2_df = vcf_file_2.to_df()

    assert {'A'} == set(vcf_file_1_df['ref'])
    assert {'AGGAGTC'} == set(vcf_file_2_df['ref'])

    sync_alleles(vcf_file_2.filepath, vcf_file_1.filepath, tmp_path / 'after')

    vcf1_observed = Vcf(tmp_dir_after / 'one-sync.vcf.bgz', tmp_dir_after).to_df()
    vcf2_observed = Vcf(tmp_dir_after / 'two-sync.vcf.bgz', tmp_dir_after).to_df()

    assert {'AGGAGTC'} == set(vcf1_observed['ref'])
    assert {'AGGAGTC'} == set(vcf2_observed['ref'])
