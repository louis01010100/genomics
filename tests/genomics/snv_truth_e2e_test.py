"""End-to-end test for `genomics snv-truth` (backbone-driven exact-match rework).

A self-contained fixture is generated at runtime (no committed binaries): an all-'A'
mini reference (chr1 + chrM), a two-sample input VCF (one male, one female), a backbone
VCF with SNP families, two `cram-depth`-style depth tables, and genders/samples files.

Exercised through the real CLI:
  * exact match, biallelic  (MALE1 chr1:100 A>C matches one family member)
  * exact match, multiallelic ALT-set (MALE1 chr1:300 A>G,T matches the 2-ALT member)
  * matched member kept, siblings dropped
  * no-match fill: ONE record per family member (FEM1 fills the 3- and 2-ALT families)
  * homref vs nocall by autosome depth, and nocall on a site missing from the table
  * chrM haploid fill genotype
  * per-sample outputs: <sample>.truth.vcf.bgz (+.csi) and gzip <sample>.truth.tsv.gz
Sex-chromosome ploidy (haploid chrX/chrY, female-chrY nocall) is covered separately by
the fixploidy unit test at real GRCh38 coordinates.
"""
import gzip
import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

HERE = Path(__file__).parent
SRC = HERE.parents[1] / 'src'

pytestmark = pytest.mark.skipif(
    shutil.which('bcftools') is None or shutil.which('samtools') is None
    or shutil.which('bgzip') is None or shutil.which('tabix') is None,
    reason='bcftools/samtools/bgzip/tabix not on PATH',
)

INPUT_VCF = (
    '##fileformat=VCFv4.2\n'
    '##contig=<ID=chr1,length=400>\n'
    '##contig=<ID=chrM,length=100>\n'
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
    '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tMALE1\tFEM1\n'
    'chr1\t100\t.\tA\tC\t.\t.\t.\tGT\t0/1\t0/0\n'      # MALE1 het A/C -> matches F1 member C
    'chr1\t150\t.\tA\tG\t.\t.\t.\tGT\t0/0\t./.\n'      # MALE1 obs homref; FEM1 obs nocall (AX-9)
    'chr1\t300\t.\tA\tG,T\t.\t.\t.\tGT\t1/2\t0/0\n'    # MALE1 G/T -> matches multiallelic member
)

BACKBONE_VCF = (
    '##fileformat=VCFv4.2\n'
    '##contig=<ID=chr1,length=400>\n'
    '##contig=<ID=chrM,length=100>\n'
    '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
    # 3-ALT family: full multiallelic + 3 splits, same ID repeated
    'chr1\t100\tAX-1,AX-2,AX-3\tA\tC,G,T\t.\t.\t.\n'
    'chr1\t100\tAX-1,AX-2,AX-3\tA\tC\t.\t.\t.\n'
    'chr1\t100\tAX-1,AX-2,AX-3\tA\tG\t.\t.\t.\n'
    'chr1\t100\tAX-1,AX-2,AX-3\tA\tT\t.\t.\t.\n'
    # 2-ALT family: full multiallelic + 2 splits
    'chr1\t300\tAX-2,AX-4\tA\tG,T\t.\t.\t.\n'
    'chr1\t300\tAX-2,AX-4\tA\tG\t.\t.\t.\n'
    'chr1\t300\tAX-2,AX-4\tA\tT\t.\t.\t.\n'
    # single-ALT families for the fill branch
    'chr1\t150\tAX-9\tA\tG\t.\t.\t.\n'      # depth 20 >= 2 -> homref
    'chr1\t200\tAX-10\tA\tT\t.\t.\t.\n'     # depth 1 < 2   -> nocall
    'chr1\t250\tAX-50\tA\tC\t.\t.\t.\n'     # absent from table -> nocall
    'chrM\t50\tAX-40\tA\tG\t.\t.\t.\n'      # chrM depth 900 -> haploid homref
)

AUTOSOMES_DEPTH = (
    'chrom\tpos\tdepth_mean\tn_samples\n'
    'chr1\t100\t25.0\t2\n'
    'chr1\t150\t20.0\t2\n'
    'chr1\t200\t1.0\t2\n'
    'chr1\t300\t25.0\t2\n'
    'chrM\t50\t900.0\t2\n'
)

SEX_DEPTH = (
    'chrom\tpos\tmean_male\tn_male\tmean_female\tn_female\n'
    'chrX\t3000000\t10.0\t1\t15.0\t1\n'
    'chrY\t3000000\t8.0\t1\t0.0\t1\n'
)

GENDERS = 'sample\tgender\nMALE1\tmale\nFEM1\tfemale\n'
SAMPLES = 'sample\nMALE1\nFEM1\n'


def _index_depths(tsv_file: Path, text: str):
    """bgzip + tabix-index a (chrom, pos) depth table like `cram-depth` produces."""
    tsv_file.write_text(text)
    bgz_file = Path(f'{tsv_file}.bgz')
    subprocess.run(f'bgzip -f -c {tsv_file} > {bgz_file}', shell=True, check=True)
    subprocess.run(
        ['tabix', '-f', '-s', '1', '-b', '2', '-e', '2', '-S', '1', str(bgz_file)],
        check=True,
    )


def _build_fixture(d: Path):
    d.mkdir(parents=True, exist_ok=True)
    # all-'A' mini reference so every SNV REF ('A') matches the reference
    genome = d / 'genome.fa'
    with genome.open('w') as fh:
        for contig, length in (('chr1', 400), ('chrM', 100)):
            fh.write(f'>{contig}\n')
            seq = 'A' * length
            for i in range(0, length, 60):
                fh.write(seq[i:i + 60] + '\n')
    subprocess.run(['samtools', 'faidx', str(genome)], check=True)

    (d / 'backbone.vcf').write_text(BACKBONE_VCF)
    _index_depths(d / 'autosomes-depth.tsv', AUTOSOMES_DEPTH)
    _index_depths(d / 'sex-depth.tsv', SEX_DEPTH)
    (d / 'genders.tsv').write_text(GENDERS)
    (d / 'samples.tsv').write_text(SAMPLES)

    plain = d / 'input.vcf'
    plain.write_text(INPUT_VCF)
    subprocess.run(['bgzip', '-f', str(plain)], check=True)   # -> input.vcf.gz
    (d / 'input.vcf.gz').rename(d / 'input.vcf.bgz')          # Vcf expects .vcf.bgz
    subprocess.run(['bcftools', 'index', '-f', str(d / 'input.vcf.bgz')], check=True)
    return genome


def _run(work: Path):
    fixture = work / 'fixture'
    genome = _build_fixture(fixture)
    out = work / 'out'
    env = dict(os.environ)
    env['PYTHONPATH'] = str(SRC)
    subprocess.run(
        [sys.executable, '-m', 'genomics', 'snv-truth',
         '--coordinates-file', str(fixture / 'backbone.vcf'),
         '--samples-file', str(fixture / 'samples.tsv'),
         '--genders-file', str(fixture / 'genders.tsv'),
         '--autosomes-depths-file', str(fixture / 'autosomes-depth.tsv.bgz'),
         '--sex-depths-file', str(fixture / 'sex-depth.tsv.bgz'),
         '--genome-file', str(genome),
         '--output-dir', str(out),
         str(fixture / 'input.vcf.bgz')],
        check=True, env=env, capture_output=True, text=True,
    )
    return out


@pytest.fixture(scope='module')
def outputs(tmp_path_factory):
    return _run(tmp_path_factory.mktemp('snv_truth'))


def _read_tsv(path: Path):
    with gzip.open(path, 'rt') as fh:
        rows = [line.rstrip('\n').split('\t') for line in fh]
    header, data = rows[0], rows[1:]
    return header, [dict(zip(header, r)) for r in data]


def test_output_files_exist(outputs):
    for s in ('MALE1', 'FEM1'):
        assert (outputs / f'{s}.truth.vcf.bgz').exists()
        assert (outputs / f'{s}.truth.vcf.bgz.csi').exists()
        assert (outputs / f'{s}.truth.tsv.gz').exists()


def test_tsv_columns(outputs):
    header, _ = _read_tsv(outputs / 'MALE1.truth.tsv.gz')
    assert header == ['chrom', 'pos', 'id', 'ref', 'alt', 'tgt']


def test_no_legacy_artifacts(outputs):
    names = {p.name for p in outputs.iterdir()}
    assert not any('snv_profile' in n for n in names)
    assert not any('truth_synced' in n for n in names)
    assert not (outputs / 'truth.vcf.bgz').exists()  # no --merge-vcf combined output


def test_male_matched_members_and_fills(outputs):
    _, rows = _read_tsv(outputs / 'MALE1.truth.tsv.gz')
    by = {(r['id'], r['pos'], r['alt']): r for r in rows}

    # biallelic match: matched member kept with sample GT, siblings dropped
    fam1 = [r for r in rows if r['id'] == 'AX-1,AX-2,AX-3']
    assert len(fam1) == 1
    assert fam1[0]['alt'] == 'C' and fam1[0]['tgt'] == 'A/C'

    # multiallelic ALT-set match
    fam2 = [r for r in rows if r['id'] == 'AX-2,AX-4']
    assert len(fam2) == 1
    assert set(fam2[0]['alt'].split(',')) == {'G', 'T'} and fam2[0]['tgt'] == 'G/T'

    # autosome fills: homref (depth>=2), nocall (depth<2), nocall (missing depth)
    assert by[('AX-9', '150', 'G')]['tgt'] == 'A/A'
    assert by[('AX-10', '200', 'T')]['tgt'] == './.'
    assert by[('AX-50', '250', 'C')]['tgt'] == './.'
    # chrM haploid homref
    assert by[('AX-40', '50', 'G')]['tgt'] == 'A'

    assert len(rows) == 6


def test_female_fill_one_record_per_member(outputs):
    _, rows = _read_tsv(outputs / 'FEM1.truth.tsv.gz')

    # no sample record -> every member of each multi-member family is emitted
    fam1 = [r for r in rows if r['id'] == 'AX-1,AX-2,AX-3']
    assert len(fam1) == 4
    assert all(r['tgt'] == 'A/A' for r in fam1)          # chr1:100 depth 25 -> homref
    assert {r['alt'] for r in fam1} == {'C,G,T', 'C', 'G', 'T'}

    fam2 = [r for r in rows if r['id'] == 'AX-2,AX-4']
    assert len(fam2) == 3
    assert all(r['tgt'] == 'A/A' for r in fam2)

    by = {r['id']: r for r in rows}
    # observed ./. at chr1:150 is honored as nocall, overriding the depth homref
    assert by['AX-9']['tgt'] == './.'
    assert by['AX-10']['tgt'] == './.'
    assert by['AX-50']['tgt'] == './.'
    assert by['AX-40']['tgt'] == 'A'                     # chrM haploid homref

    assert len(rows) == 11


def test_truth_vcf_sorted_and_indexed(outputs):
    # index is usable and the VCF is coordinate-sorted
    out = subprocess.run(
        ['bcftools', 'query', '-f', '%CHROM\t%POS\n', str(outputs / 'MALE1.truth.vcf.bgz')],
        check=True, capture_output=True, text=True,
    ).stdout.strip().split('\n')
    positions = [(c, int(p)) for c, p in (line.split('\t') for line in out)]
    assert positions == sorted(positions, key=lambda x: (x[0], x[1]))
