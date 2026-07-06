#!/usr/bin/env python
import logging
import subprocess
from pathlib import Path

import pytest

from genomics.samples import export_samples

# Two source VCFs carrying QUAL/FILTER/INFO and an extra FORMAT (DP) to prove
# stripping. Each source is individually coordinate-sorted, but the sources
# INTERLEAVE across positions so the per-sample concat must sort to order them:
#   source1: samples A, B  (chr1:100, chr1:300)
#   source2: samples B, C  (chr1:200, chr2:100)
# B spans both sources (its records interleave: 100,300 | 200,100@chr2);
# A only source1; C only source2.

SOURCE1 = """\
##fileformat=VCFv4.2
##contig=<ID=chr1>
##contig=<ID=chr2>
##FILTER=<ID=PASS,Description="passed">
##INFO=<ID=AC,Number=A,Type=Integer,Description="ac">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="depth">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA\tB
chr1\t100\trs1\tC\tT\t60\tPASS\tAC=2\tGT:DP\t0/0:25\t0/1:15
chr1\t300\trs3\tA\tG\t50\tPASS\tAC=1\tGT:DP\t0/1:30\t1/1:20
"""

SOURCE2 = """\
##fileformat=VCFv4.2
##contig=<ID=chr1>
##contig=<ID=chr2>
##FILTER=<ID=PASS,Description="passed">
##INFO=<ID=AC,Number=A,Type=Integer,Description="ac">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="depth">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tB\tC
chr1\t200\trs2\tG\tA\t40\tPASS\tAC=1\tGT:DP\t0/1:22\t1/1:18
chr2\t100\trs4\tT\tC\t30\tPASS\tAC=1\tGT:DP\t0/0:10\t0/1:12
"""


def _vcf(out, s):
    return out / 'samples' / f'{s}.vcf.bgz'


def _bcftools(*args):
    return subprocess.run(['bcftools', *args], capture_output=True, text=True,
                          check=True).stdout


def _build_sources(root: Path):
    (root / 'source1.vcf').write_text(SOURCE1)
    (root / 'source2.vcf').write_text(SOURCE2)
    return [root / 'source1.vcf', root / 'source2.vcf']


def _samples_file(root: Path, names):
    p = root / 'samples.txt'
    p.write_text('sample_name\n' + '\n'.join(names) + '\n')
    return p


def test_export_samples_worked_example(tmp_path):
    sources = _build_sources(tmp_path)
    # 'B' duplicated -> must collapse to one output.
    samples_file = _samples_file(tmp_path, ['A', 'B', 'B', 'C'])
    out = tmp_path / 'out'

    export_samples(vcf_files=sources, samples_file=samples_file,
                   output_dir=out, n_threads=1)

    # Flat samples/ folder; duplicate B collapsed. Top level is only samples/
    # (tmp/ removed on success).
    assert sorted(p.name for p in out.iterdir()) == ['samples']
    sdir = out / 'samples'
    for s in ('A', 'B', 'C'):
        assert (sdir / f'{s}.vcf.bgz').exists()
        assert list(sdir.glob(f'{s}.vcf.bgz.csi')) or \
            list(sdir.glob(f'{s}.vcf.bgz.tbi')), 'index missing'
    # samples/ holds only flat per-sample VCFs + indexes -- no subdirectories.
    assert all(p.is_file() for p in sdir.iterdir())

    # Single sample column each.
    assert _bcftools('query', '-l', str(_vcf(out, 'B'))).split() == ['B']
    assert _bcftools('query', '-l', str(_vcf(out, 'A'))).split() == ['A']

    # B = concat across both sources, sorted by CHROM,POS.
    b_rows = _bcftools('query', '-f', '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n',
                       str(_vcf(out, 'B'))).strip().split('\n')
    assert b_rows == [
        'chr1\t100\tC\tT\t0/1',
        'chr1\t200\tG\tA\t0/1',
        'chr1\t300\tA\tG\t1/1',
        'chr2\t100\tT\tC\t0/0',
    ]

    # A only from source1; C only from source2.
    a_rows = _bcftools('query', '-f', '%CHROM\t%POS[\t%GT]\n',
                       str(_vcf(out, 'A'))).strip().split('\n')
    assert a_rows == ['chr1\t100\t0/0', 'chr1\t300\t0/1']
    c_rows = _bcftools('query', '-f', '%CHROM\t%POS[\t%GT]\n',
                       str(_vcf(out, 'C'))).strip().split('\n')
    assert c_rows == ['chr1\t200\t1/1', 'chr2\t100\t0/1']

    # GT-only: FORMAT column is exactly GT; QUAL/FILTER/INFO stripped to '.'.
    full = _bcftools('view', '-H', str(_vcf(out, 'B'))).strip().split('\n')
    for line in full:
        f = line.split('\t')
        assert f[5] == '.', f'QUAL not stripped: {line}'
        assert f[6] == '.', f'FILTER not stripped: {line}'
        assert f[7] == '.', f'INFO not stripped: {line}'
        assert f[8] == 'GT', f'FORMAT not GT-only: {line}'
    hdr = _bcftools('view', '-h', str(_vcf(out, 'B')))
    assert '##FORMAT=<ID=DP' not in hdr


def test_warn_and_continue_on_absent(tmp_path, caplog):
    sources = _build_sources(tmp_path)
    samples_file = _samples_file(tmp_path, ['A', 'Z'])  # Z is in neither source
    out = tmp_path / 'out'

    with caplog.at_level(logging.WARNING):
        export_samples(vcf_files=sources, samples_file=samples_file,
                       output_dir=out, n_threads=1)

    # Present sample A is produced; absent Z is skipped with a warning; success.
    assert (out / 'samples' / 'A.vcf.bgz').exists()
    assert not (out / 'samples' / 'Z.vcf.bgz').exists()
    assert sorted(p.name for p in out.iterdir()) == ['samples']
    assert [p.name for p in (out / 'samples').glob('*.vcf.bgz')] == ['A.vcf.bgz']
    assert 'Z' in caplog.text


def test_samples_file_requires_sample_name_header(tmp_path):
    sources = _build_sources(tmp_path)
    # A file without the 'sample_name' header is rejected.
    headerless = tmp_path / 'nohdr.txt'
    headerless.write_text('A\nB\n')
    out = tmp_path / 'out'

    with pytest.raises(Exception):
        export_samples(vcf_files=sources, samples_file=headerless,
                       output_dir=out, n_threads=1)


def test_error_when_all_absent(tmp_path):
    sources = _build_sources(tmp_path)
    samples_file = _samples_file(tmp_path, ['Y', 'Z'])  # both absent from all
    out = tmp_path / 'out'

    with pytest.raises(Exception):
        export_samples(vcf_files=sources, samples_file=samples_file,
                       output_dir=out, n_threads=1)

    # No output published.
    assert not list(out.glob('samples/*.vcf.bgz')) if out.exists() else True


def test_deterministic_across_threads(tmp_path):
    sources = _build_sources(tmp_path)
    samples_file = _samples_file(tmp_path, ['A', 'B', 'C'])

    out1 = tmp_path / 'out1'
    out2 = tmp_path / 'out2'
    export_samples(vcf_files=sources, samples_file=samples_file,
                   output_dir=out1, n_threads=1)
    export_samples(vcf_files=sources, samples_file=samples_file,
                   output_dir=out2, n_threads=2)

    for s in ('A', 'B', 'C'):
        r1 = _bcftools('view', '-H', str(_vcf(out1, s)))
        r2 = _bcftools('view', '-H', str(_vcf(out2, s)))
        assert r1 == r2
