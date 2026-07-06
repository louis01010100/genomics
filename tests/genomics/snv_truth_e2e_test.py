"""End-to-end test for `genomics snv-truth` (snv-family-driven exact-match rework).

A self-contained fixture is generated at runtime (no committed binaries): an all-'A'
mini reference (chr1 + chrM), a two-sample input VCF (one male, one female), a snv-family
gzip TSV of SNP families (id = `FM-<n>`), two `cram-depth`-style depth tables, and
genders/samples files.

Exercised through the real CLI:
  * exact match, biallelic  (MALE1 chr1:100 A>C matches one family member)
  * exact match, multiallelic ALT-set (MALE1 chr1:300 A>G,T matches the 2-ALT member)
  * matched member kept, siblings dropped
  * no-match fill: ONE record per family member (FEM1 fills the 3- and 2-ALT families)
  * homref vs nocall by autosome depth, and nocall on a site missing from the table
  * chrM haploid fill genotype
  * per-sample VCF outputs: <sample>.truth.vcf.bgz (+.csi)
  * one combined gzip TSV for all samples: truth.tsv.gz (sample_name-tagged rows)
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
    'chrM\t10\t.\tA\tG\t.\t.\t.\tGT\t0/1\t0/1\n'       # non-backbone chrM record (keeps the
)                                                      # per-chrom chrM input non-empty)

# snv-family TSV (as `axiom snv-family` emits): ID column is the family id `FM-<n>`,
# probesets in `psid`. Same families/records as the old backbone fixture, keyed by fmid:
#   FM-1 = chr1:100 (AX-1,AX-2,AX-3)   FM-2 = chr1:150 (AX-9)   FM-3 = chr1:200 (AX-10)
#   FM-4 = chr1:250 (AX-50)   FM-5 = chr1:300 (AX-2,AX-4)   FM-6 = chrM:50 (AX-40)
SNV_FAMILY_TSV = (
    'chrom\tpos\tref\talt\tfmid\tmsid\tasid\tpsid\n'
    # 3-ALT family: full multiallelic + 3 splits, one fmid
    'chr1\t100\tA\tC,G,T\tFM-1\tNA\tAffx-1,Affx-2,Affx-3\tAX-1,AX-2,AX-3\n'
    'chr1\t100\tA\tC\tFM-1\tNA\tAffx-1,Affx-2,Affx-3\tAX-1,AX-2,AX-3\n'
    'chr1\t100\tA\tG\tFM-1\tNA\tAffx-1,Affx-2,Affx-3\tAX-1,AX-2,AX-3\n'
    'chr1\t100\tA\tT\tFM-1\tNA\tAffx-1,Affx-2,Affx-3\tAX-1,AX-2,AX-3\n'
    # 2-ALT family: full multiallelic + 2 splits
    'chr1\t300\tA\tG,T\tFM-5\tNA\tAffx-2,Affx-4\tAX-2,AX-4\n'
    'chr1\t300\tA\tG\tFM-5\tNA\tAffx-2,Affx-4\tAX-2,AX-4\n'
    'chr1\t300\tA\tT\tFM-5\tNA\tAffx-2,Affx-4\tAX-2,AX-4\n'
    # single-ALT families for the fill branch
    'chr1\t150\tA\tG\tFM-2\tNA\tAffx-9\tAX-9\n'       # depth 20 >= 2 -> homref
    'chr1\t200\tA\tT\tFM-3\tNA\tAffx-10\tAX-10\n'     # depth 1 < 2   -> nocall
    'chr1\t250\tA\tC\tFM-4\tNA\tAffx-50\tAX-50\n'     # absent from table -> nocall
    'chrM\t50\tA\tG\tFM-6\tNA\tAffx-40\tAX-40\n'      # chrM depth 900 -> haploid homref
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

    with gzip.open(d / 'snv-family.tsv.gz', 'wt') as fh:
        fh.write(SNV_FAMILY_TSV)
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


def _run_snv_truth(fixture: Path, genome: Path, out: Path, vcf_args, n_threads=1, check=True):
    env = dict(os.environ)
    env['PYTHONPATH'] = str(SRC)
    return subprocess.run(
        [sys.executable, '-m', 'genomics', 'snv-truth',
         '--snv-family-file', str(fixture / 'snv-family.tsv.gz'),
         '--samples-file', str(fixture / 'samples.tsv'),
         '--genders-file', str(fixture / 'genders.tsv'),
         '--autosomes-depths-file', str(fixture / 'autosomes-depth.tsv.bgz'),
         '--sex-depths-file', str(fixture / 'sex-depth.tsv.bgz'),
         '--genome-file', str(genome),
         '--output-dir', str(out),
         '--n-threads', str(n_threads)] + [str(v) for v in vcf_args],
        check=check, env=env, capture_output=True, text=True,
    )


def _run(work: Path):
    fixture = work / 'fixture'
    genome = _build_fixture(fixture)
    out = work / 'out'
    # inputs must be one chromosome per file
    _run_snv_truth(fixture, genome, out, _split_input_by_chrom(fixture))
    return out


def _split_input_by_chrom(fixture: Path):
    """Split the single multi-sample input into per-chromosome, indexed VCFs (each
    still multi-sample) — a by-chromosome split of the same data."""
    src = fixture / 'input.vcf.bgz'
    paths = []
    for chrom in ('chr1', 'chrM'):
        p = fixture / f'{chrom}.vcf.bgz'
        subprocess.run(f'bcftools view -r {chrom} -O z -o {p} {src}', shell=True, check=True)
        subprocess.run(['bcftools', 'index', '-f', str(p)], check=True)
        paths.append(p)
    return paths


@pytest.fixture(scope='module')
def outputs(tmp_path_factory):
    return _run(tmp_path_factory.mktemp('snv_truth'))


def _read_tsv(path: Path):
    with gzip.open(path, 'rt') as fh:
        rows = [line.rstrip('\n').split('\t') for line in fh]
    header, data = rows[0], rows[1:]
    return header, [dict(zip(header, r)) for r in data]


def _rows_for(outputs, sample):
    _, rows = _read_tsv(outputs / 'truth.tsv.gz')
    return [r for r in rows if r['sample_name'] == sample]


def test_output_files_exist(outputs):
    assert (outputs / 'truth.tsv.gz').exists()            # one combined tsv
    for s in ('MALE1', 'FEM1'):
        assert (outputs / f'{s}.truth.vcf.bgz').exists()
        assert (outputs / f'{s}.truth.vcf.bgz.csi').exists()
        assert not (outputs / f'{s}.truth.tsv.gz').exists()  # per-sample tsvs removed


def test_tsv_columns(outputs):
    header, _ = _read_tsv(outputs / 'truth.tsv.gz')
    assert header == ['sample_name', 'chrom', 'pos', 'id', 'ref', 'alt', 'tgt']


def test_combined_tsv_groups_samples(outputs):
    _, rows = _read_tsv(outputs / 'truth.tsv.gz')
    names = [r['sample_name'] for r in rows]
    assert set(names) == {'MALE1', 'FEM1'}
    # samples emitted in sorted order, grouped in contiguous blocks
    assert names == sorted(names)
    assert len(rows) == 6 + 11                            # MALE1 + FEM1 rows


def test_no_legacy_artifacts(outputs):
    names = {p.name for p in outputs.iterdir()}
    assert not any('snv_profile' in n for n in names)
    assert not any('truth_synced' in n for n in names)
    assert not (outputs / 'truth.vcf.bgz').exists()  # no --merge-vcf combined output


def test_male_matched_members_and_fills(outputs):
    rows = _rows_for(outputs, 'MALE1')
    by = {(r['id'], r['pos'], r['alt']): r for r in rows}

    # biallelic match: matched member kept with sample GT, siblings dropped
    fam1 = [r for r in rows if r['id'] == 'FM-1']
    assert len(fam1) == 1
    assert fam1[0]['alt'] == 'C' and fam1[0]['tgt'] == 'A/C'

    # multiallelic ALT-set match
    fam2 = [r for r in rows if r['id'] == 'FM-5']
    assert len(fam2) == 1
    assert set(fam2[0]['alt'].split(',')) == {'G', 'T'} and fam2[0]['tgt'] == 'G/T'

    # autosome fills: homref (depth>=2), nocall (depth<2), nocall (missing depth)
    assert by[('FM-2', '150', 'G')]['tgt'] == 'A/A'
    assert by[('FM-3', '200', 'T')]['tgt'] == './.'
    assert by[('FM-4', '250', 'C')]['tgt'] == './.'
    # chrM haploid homref
    assert by[('FM-6', '50', 'G')]['tgt'] == 'A'

    assert len(rows) == 6


def test_female_fill_one_record_per_member(outputs):
    rows = _rows_for(outputs, 'FEM1')

    # no sample record -> every member of each multi-member family is emitted
    fam1 = [r for r in rows if r['id'] == 'FM-1']
    assert len(fam1) == 4
    assert all(r['tgt'] == 'A/A' for r in fam1)          # chr1:100 depth 25 -> homref
    assert {r['alt'] for r in fam1} == {'C,G,T', 'C', 'G', 'T'}

    fam2 = [r for r in rows if r['id'] == 'FM-5']
    assert len(fam2) == 3
    assert all(r['tgt'] == 'A/A' for r in fam2)

    by = {r['id']: r for r in rows}
    # observed ./. at chr1:150 is honored as nocall, overriding the depth homref
    assert by['FM-2']['tgt'] == './.'
    assert by['FM-3']['tgt'] == './.'
    assert by['FM-4']['tgt'] == './.'
    assert by['FM-6']['tgt'] == 'A'                      # chrM haploid homref

    assert len(rows) == 11


def test_truth_vcf_sorted_and_indexed(outputs):
    # index is usable and the VCF is coordinate-sorted
    out = subprocess.run(
        ['bcftools', 'query', '-f', '%CHROM\t%POS\n', str(outputs / 'MALE1.truth.vcf.bgz')],
        check=True, capture_output=True, text=True,
    ).stdout.strip().split('\n')
    positions = [(c, int(p)) for c, p in (line.split('\t') for line in out)]
    assert positions == sorted(positions, key=lambda x: (x[0], x[1]))


def test_determinism_input_order_and_threads(tmp_path):
    # Same per-chromosome inputs in different order and thread counts -> byte-identical
    # output (AC5). n_threads > sample count also exercises per-sample parallelism (AC4).
    fixture = tmp_path / 'fixture'
    genome = _build_fixture(fixture)
    split_vcfs = _split_input_by_chrom(fixture)

    a = tmp_path / 'a'
    _run_snv_truth(fixture, genome, a, split_vcfs, n_threads=1)
    b = tmp_path / 'b'
    _run_snv_truth(fixture, genome, b, list(reversed(split_vcfs)), n_threads=8)

    assert gzip.open(a / 'truth.tsv.gz', 'rt').read() == gzip.open(b / 'truth.tsv.gz', 'rt').read()
    for sample in ('MALE1', 'FEM1'):
        assert (b / f'{sample}.truth.vcf.bgz').exists()
        assert (b / f'{sample}.truth.vcf.bgz.csi').exists()


def test_extract_once_per_input(tmp_path):
    # Each input (chromosome) is extracted exactly once: one extraction dir per input
    # file, not one per (sample, chromosome) (AC1).
    fixture = tmp_path / 'fixture'
    genome = _build_fixture(fixture)
    split_vcfs = _split_input_by_chrom(fixture)   # two inputs: chr1, chrM
    out = tmp_path / 'out'
    _run_snv_truth(fixture, genome, out, split_vcfs)
    chrom_dirs = sorted(p.name for p in (out / 'tmp' / 'extract').iterdir() if p.is_dir())
    assert chrom_dirs == ['chr1', 'chrM']         # one per input, not 4 (per sample x chrom)


def test_backbone_file_flag_removed(tmp_path):
    # the old --backbone-file option is gone; the CLI rejects it
    fixture = tmp_path / 'fixture'
    genome = _build_fixture(fixture)
    result = subprocess.run(
        [sys.executable, '-m', 'genomics', 'snv-truth',
         '--backbone-file', str(fixture / 'snv-family.tsv.gz'),
         '--samples-file', str(fixture / 'samples.tsv'),
         '--genders-file', str(fixture / 'genders.tsv'),
         '--autosomes-depths-file', str(fixture / 'autosomes-depth.tsv.bgz'),
         '--sex-depths-file', str(fixture / 'sex-depth.tsv.bgz'),
         '--genome-file', str(genome),
         '--output-dir', str(tmp_path / 'out')],
        env={**os.environ, 'PYTHONPATH': str(SRC)}, capture_output=True, text=True,
    )
    # fails: --backbone-file is not recognized and does not satisfy the new required flag
    assert result.returncode != 0
    assert '--snv-family-file' in result.stderr


def test_multichromosome_input_fails(tmp_path):
    # An input spanning more than one chromosome is rejected up front, no output (AC2).
    fixture = tmp_path / 'fixture'
    genome = _build_fixture(fixture)   # builds the indexed multi-chrom input.vcf.bgz
    out = tmp_path / 'out'
    result = _run_snv_truth(fixture, genome, out, [fixture / 'input.vcf.bgz'], check=False)
    assert result.returncode != 0
    assert 'exactly one chromosome' in (result.stderr + result.stdout)
    assert not (out / 'truth.tsv.gz').exists()


def test_unindexed_input_fails(tmp_path):
    # An input lacking an index aborts with an actionable error, before any output.
    fixture = tmp_path / 'fixture'
    genome = _build_fixture(fixture)
    split_vcfs = _split_input_by_chrom(fixture)

    noidx = fixture / 'noidx.vcf.bgz'
    shutil.copy(split_vcfs[0], noidx)   # copy a per-chrom input WITHOUT its index

    out = tmp_path / 'out'
    result = _run_snv_truth(fixture, genome, out, [noidx, split_vcfs[1]], check=False)
    assert result.returncode != 0
    assert 'not indexed' in (result.stderr + result.stdout)
    assert not (out / 'truth.tsv.gz').exists()


def test_unindexed_depth_table_fails(tmp_path):
    # A depth table lacking its tabix index aborts fail-fast (before any output),
    # since depth is validated up front even though it is fetched lazily.
    fixture = tmp_path / 'fixture'
    genome = _build_fixture(fixture)
    (fixture / 'autosomes-depth.tsv.bgz.tbi').unlink()   # drop the tabix index

    out = tmp_path / 'out'
    result = _run_snv_truth(fixture, genome, out, _split_input_by_chrom(fixture), check=False)
    assert result.returncode != 0
    assert 'not indexed' in (result.stderr + result.stdout)
    assert not (out / 'truth.tsv.gz').exists()
