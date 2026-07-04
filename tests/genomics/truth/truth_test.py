"""Unit tests for the reworked `snv-truth` pure helpers in genomics.truth.

These target the backbone-driven exact-match rework: reading the backbone VCF as
coordinates, grouping into SNP families by ID, loading the two `cram-depth` tables,
the contig-class/gender-aware depth selection, the exact-match key, and the fill
homref/nocall decision.
"""
import gzip
import math
import shutil
from pathlib import Path

import pytest

from genomics.truth import (
    load_coordinates,
    group_families,
    load_autosomes_depths,
    load_sex_depths,
    depth_for,
    match_key,
    is_homref,
    fixploidy_lines,
    GRCH38_PLOIDY,
)

BACKBONE = (
    '##fileformat=VCFv4.2\n'
    '##contig=<ID=chr1>\n'
    '##contig=<ID=chrX>\n'
    '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
    # a 3-ALT family: full multiallelic + its splits, same ID string repeated
    'chr1\t930165\tAX-1,AX-2,AX-3\tA\tC,G,T\t.\t.\t.\n'
    'chr1\t930165\tAX-1,AX-2,AX-3\tA\tc\t.\t.\t.\n'
    'chr1\t930165\tAX-1,AX-2,AX-3\tA\tG\t.\t.\t.\n'
    'chr1\t930165\tAX-1,AX-2,AX-3\tA\tT\t.\t.\t.\n'
    # a lone single-ALT family
    'chr1\t11080559\tAX-9\tAAAAAAC\tA\t.\t.\t.\n'
)


def _write(path: Path, text: str, gz: bool):
    if gz:
        with gzip.open(path, 'wt') as fh:
            fh.write(text)
    else:
        path.write_text(text)
    return path


@pytest.mark.parametrize('gz,name', [(True, 'bb.vcf.bgz'), (False, 'bb.vcf')])
def test_load_coordinates_reads_vcf(tmp_path, gz, name):
    f = _write(tmp_path / name, BACKBONE, gz)
    rows = load_coordinates(f)
    assert len(rows) == 5
    r = rows[0]
    assert r['chrom'] == 'chr1'
    assert r['pos'] == 930165 and isinstance(r['pos'], int)
    assert r['id'] == 'AX-1,AX-2,AX-3'
    assert r['ref'] == 'A'
    assert r['alt'] == 'C,G,T'
    # lowercase alt is uppercased
    assert rows[1]['alt'] == 'C'


def test_group_families_by_id_first_seen_order(tmp_path):
    rows = load_coordinates(_write(tmp_path / 'bb.vcf', BACKBONE, False))
    fams = group_families(rows)
    assert list(fams.keys()) == ['AX-1,AX-2,AX-3', 'AX-9']
    assert len(fams['AX-1,AX-2,AX-3']) == 4
    assert len(fams['AX-9']) == 1


def test_load_autosomes_depths(tmp_path):
    f = tmp_path / 'autosomes-depth.tsv'
    f.write_text(
        'chrom\tpos\tdepth_mean\tn_samples\n'
        'chr1\t100\t12.5\t4\n'
        'chrM\t5\t900.0\t4\n'
    )
    d = load_autosomes_depths(f)
    assert d[('chr1', 100)] == pytest.approx(12.5)
    assert d[('chrM', 5)] == pytest.approx(900.0)


def test_load_sex_depths_blank_maps_to_none(tmp_path):
    f = tmp_path / 'sex-depth.tsv'
    f.write_text(
        'chrom\tpos\tmean_male\tn_male\tmean_female\tn_female\n'
        'chrX\t100\t5.0\t2\t11.0\t2\n'
        'chrY\t200\t8.0\t2\t\t0\n'   # no females in cohort -> blank
    )
    d = load_sex_depths(f)
    assert d[('chrX', 100)]['male'] == pytest.approx(5.0)
    assert d[('chrX', 100)]['female'] == pytest.approx(11.0)
    assert d[('chrY', 200)]['male'] == pytest.approx(8.0)
    assert d[('chrY', 200)]['female'] is None


@pytest.fixture
def depth_tables(tmp_path):
    a = tmp_path / 'autosomes-depth.tsv'
    a.write_text(
        'chrom\tpos\tdepth_mean\tn_samples\n'
        'chr1\t100\t12.5\t4\n'
        'chrM\t5\t900.0\t4\n'
    )
    s = tmp_path / 'sex-depth.tsv'
    s.write_text(
        'chrom\tpos\tmean_male\tn_male\tmean_female\tn_female\n'
        'chrX\t100\t5.0\t2\t11.0\t2\n'
        'chrY\t200\t8.0\t2\t\t0\n'
    )
    return load_autosomes_depths(a), load_sex_depths(s)


def test_depth_for_autosome_and_mito(depth_tables):
    auto, sex = depth_tables
    assert depth_for('chr1', 100, 'male', auto, sex) == pytest.approx(12.5)
    assert depth_for('chrM', 5, 'female', auto, sex) == pytest.approx(900.0)


def test_depth_for_sex_by_gender(depth_tables):
    auto, sex = depth_tables
    assert depth_for('chrX', 100, 'male', auto, sex) == pytest.approx(5.0)
    assert depth_for('chrX', 100, 'female', auto, sex) == pytest.approx(11.0)
    assert depth_for('chrY', 200, 'male', auto, sex) == pytest.approx(8.0)
    # blank female column -> None (below threshold -> nocall)
    assert depth_for('chrY', 200, 'female', auto, sex) is None


def test_depth_for_missing_key_is_none(depth_tables):
    auto, sex = depth_tables
    assert depth_for('chr1', 999, 'male', auto, sex) is None
    assert depth_for('chrX', 999, 'female', auto, sex) is None


def test_match_key_altset_order_and_case_invariant():
    assert match_key('chr1', 930165, 'A', 'C,G') == match_key('chr1', 930165, 'A', 'g,c')
    assert match_key('chr1', '930165', 'A', 'C') == match_key('chr1', 930165, 'a', 'c')
    assert match_key('chr1', 930165, 'A', 'C,G') != match_key('chr1', 930165, 'A', 'C')
    assert match_key('chr1', 930165, 'A', 'C') != match_key('chr1', 930165, 'A', 'T')


def test_is_homref_threshold():
    assert is_homref(4.0, 2) is True
    assert is_homref(2.0, 2) is True
    assert is_homref(1.9, 2) is False
    assert is_homref(None, 2) is False
    assert is_homref(math.nan, 2) is False


def test_fill_depth_is_per_member_position():
    """Members of one family can sit at different positions (axiom backbone
    re-normalizes/left-aligns expanded records), so the homref/nocall decision must
    be looked up per member, not once for the whole family."""
    from genomics.truth import call_families
    families = {
        'AX-1,AX-2': [
            {'chrom': 'chr1', 'pos': 100, 'id': 'AX-1,AX-2', 'ref': 'A', 'alt': 'C'},
            {'chrom': 'chr1', 'pos': 500, 'id': 'AX-1,AX-2', 'ref': 'A', 'alt': 'G'},
        ]
    }
    auto = {('chr1', 100): 10.0, ('chr1', 500): 1.0}   # pos 100 homref, pos 500 nocall
    matched, fills = call_families(families, {}, 'male', auto, {}, min_depth=2)
    assert matched == []
    gt = {line.split('\t')[1]: line.split('\t')[9] for line in fills}
    assert gt['100'] == '0/0'
    assert gt['500'] == './.'


HEADER = (
    '##fileformat=VCFv4.2\n'
    '##contig=<ID=chr1,length=250000000>\n'
    '##contig=<ID=chrX,length=156040895>\n'
    '##contig=<ID=chrY,length=57227415>\n'
    '##contig=<ID=chrM,length=16569>\n'
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
    '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n'
)


@pytest.mark.skipif(shutil.which('bcftools') is None, reason='bcftools not on PATH')
def test_fixploidy_grch38_gender_aware(tmp_path):
    """Fill records get GRCh38 gender-aware ploidy at real non-PAR coordinates via
    the module's GRCH38_PLOIDY: male sex-chrom/chrM haploid, female chrY nocall."""
    header_file = tmp_path / 'header.vcf'
    header_file.write_text(HEADER)
    ploidy_file = tmp_path / 'grch38.ploidy'
    ploidy_file.write_text(GRCH38_PLOIDY)

    fills = [
        'chr1\t100\ta\tA\tC\t.\t.\t.\tGT\t0/0',       # autosome -> diploid
        'chrX\t3000000\tb\tA\tC\t.\t.\t.\tGT\t0/0',   # non-PAR chrX
        'chrY\t3000000\tc\tA\tC\t.\t.\t.\tGT\t0/0',   # MSY chrY
        'chrM\t50\td\tA\tC\t.\t.\t.\tGT\t0/0',        # mito
        'chrY\t3000001\te\tA\tC\t.\t.\t.\tGT\t./.',   # nocall on MSY
    ]

    (tmp_path / 'm').mkdir()
    (tmp_path / 'f').mkdir()
    male = fixploidy_lines(fills, header_file, 'S1', 'male', ploidy_file, tmp_path / 'm')
    female = fixploidy_lines(fills, header_file, 'S1', 'female', ploidy_file, tmp_path / 'f')

    def gt(lines, i):
        return lines[i].split('\t')[9]

    # male: autosome diploid; non-PAR chrX/chrY and chrM haploid; MSY nocall -> '.'
    assert [gt(male, i) for i in range(5)] == ['0/0', '0', '0', '0', '.']
    # female: chrX diploid (default 2), chrY absent (ploidy 0 -> nocall), chrM haploid
    assert [gt(female, i) for i in range(5)] == ['0/0', '0/0', '.', '0', '.']
