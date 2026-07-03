"""End-to-end test for `genomics cram-depth` over a small real-derived fixture.

The fixture in ``cram_depth_e2e_fixture/`` was sliced from real 1000G GRCh38
CRAMs so that four samples (2 male, 2 female) exercise every edge case in a
few thousand positions:

  * chr1  -> autosome, single mean over all samples           (autosomes-depth.tsv)
  * chrM  -> mitochondria, treated like an autosome           (autosomes-depth.tsv)
  * chrX  -> non-PAR; female (diploid) depth > male (haploid) (sex-depth.tsv)
  * chrY  -> MSY; male present, female ~0 with zero-depth rows (sex-depth.tsv)
  * chrUn_...decoy -> excluded from both outputs

Regenerate the fixture with ``cram_depth_e2e.sh`` (or the branch's build_fixture.sh).

Run:  wt/main/.venv/bin/python -m pytest cram_depth_e2e_test.py -v
"""
import csv
import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).parent
FIXTURE = ROOT / 'cram_depth_e2e_fixture'
SAMPLES = ['HG00403', 'HG00406', 'HG00404', 'HG00419']


def _resolve_src():
    """Locate the `genomics` package source under test, regardless of where this
    e2e folder lives: walk up to the repo container (holds `.bare` and `wt/`),
    then prefer the feature branch worktree, falling back to main."""
    override = os.environ.get('GENOMICS_SRC')
    if override:
        return override
    for container in (ROOT, *ROOT.parents):
        if (container / '.bare').exists() and (container / 'wt').is_dir():
            for name in ('feat-cram-depth', 'main'):
                candidate = container / 'wt' / name / 'src'
                if (candidate / 'genomics').is_dir():
                    return str(candidate)
            break
    return None


SRC = _resolve_src()

pytestmark = pytest.mark.skipif(
    shutil.which('samtools') is None or SRC is None,
    reason='samtools not on PATH or genomics source not found',
)


def _read_tsv(path):
    with path.open() as fh:
        return list(csv.DictReader(fh, delimiter='\t'))


def _run_cli(work_dir, n_threads):
    work_dir.mkdir(parents=True, exist_ok=True)
    crams_file = work_dir / 'crams.tsv'
    with crams_file.open('w') as fh:
        fh.write('sample\tcram\n')
        for s in SAMPLES:
            fh.write(f'{s}\t{FIXTURE / f"{s}.cram"}\n')

    out = work_dir / 'depth'
    env = dict(os.environ)
    env['PYTHONPATH'] = SRC
    subprocess.run(
        [sys.executable, '-m', 'genomics', 'cram-depth',
         '--crams-file', str(crams_file),
         '--genders-file', str(FIXTURE / 'genders.tsv'),
         '--genome-file', str(FIXTURE / 'miniref.fa'),
         '--output-dir', str(out),
         '--n-threads', str(n_threads)],
        check=True, env=env, capture_output=True, text=True,
    )
    return out


@pytest.fixture(scope='module')
def outputs(tmp_path_factory):
    return _run_cli(tmp_path_factory.mktemp('cram_depth_out'), n_threads=4)


def test_only_expected_files(outputs):
    assert sorted(p.name for p in outputs.iterdir()) == [
        'autosomes-depth.tsv', 'depth.log', 'sex-depth.tsv',
    ]
    assert not (outputs / 'depths').exists()


def test_autosomes_schema_and_contigs(outputs):
    rows = _read_tsv(outputs / 'autosomes-depth.tsv')
    assert list(rows[0].keys()) == ['chrom', 'pos', 'depth_mean', 'n_samples']
    by_contig = {}
    for r in rows:
        by_contig[r['chrom']] = by_contig.get(r['chrom'], 0) + 1
    # chr1 + chrM only; chrM (mitochondria) is treated like an autosome
    assert set(by_contig) == {'chr1', 'chrM'}
    assert by_contig['chr1'] == 2000 and by_contig['chrM'] == 2000
    assert all(r['n_samples'] == '4' for r in rows)


def test_sex_schema_and_contigs(outputs):
    rows = _read_tsv(outputs / 'sex-depth.tsv')
    assert list(rows[0].keys()) == [
        'chrom', 'pos', 'mean_male', 'n_male', 'mean_female', 'n_female',
    ]
    assert {r['chrom'] for r in rows} == {'chrX', 'chrY'}
    assert all(r['n_male'] == '2' and r['n_female'] == '2' for r in rows)


def test_chrM_only_in_autosomes(outputs):
    auto = {r['chrom'] for r in _read_tsv(outputs / 'autosomes-depth.tsv')}
    sex = {r['chrom'] for r in _read_tsv(outputs / 'sex-depth.tsv')}
    assert 'chrM' in auto and 'chrM' not in sex
    assert not ({'chrX', 'chrY'} & auto)


def test_chrX_female_diploid_gt_male_haploid(outputs):
    rows = [r for r in _read_tsv(outputs / 'sex-depth.tsv') if r['chrom'] == 'chrX']
    male = sum(float(r['mean_male']) for r in rows) / len(rows)
    female = sum(float(r['mean_female']) for r in rows) / len(rows)
    assert female > 1.5 * male > 0


def test_chrY_male_present_female_near_zero(outputs):
    rows = [r for r in _read_tsv(outputs / 'sex-depth.tsv') if r['chrom'] == 'chrY']
    male = sum(float(r['mean_male']) for r in rows) / len(rows)
    female = sum(float(r['mean_female']) for r in rows) / len(rows)
    zero_female = sum(1 for r in rows if float(r['mean_female']) == 0.0)
    assert male > 5.0
    assert female < male / 5
    assert zero_female > 0


def test_decoy_contig_excluded(outputs):
    for name in ('autosomes-depth.tsv', 'sex-depth.tsv'):
        assert 'decoy' not in (outputs / name).read_text()


def test_output_is_thread_invariant(tmp_path):
    single = _run_cli(tmp_path / 'n1', n_threads=1)
    multi = _run_cli(tmp_path / 'n4', n_threads=4)
    for name in ('autosomes-depth.tsv', 'sex-depth.tsv'):
        assert (single / name).read_bytes() == (multi / name).read_bytes()
