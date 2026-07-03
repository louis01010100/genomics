from genomics.depth import merge
from genomics.depth import (
    classify_contig,
    validate_genders,
    mean_single,
    mean_per_gender,
)
import numpy as np
import polars as pl
import pytest


def test_classify_contig():
    for i in range(1, 23):
        assert classify_contig(f'chr{i}') == 'single'
    assert classify_contig('chrMT') == 'single'
    assert classify_contig('chrX') == 'sex'
    assert classify_contig('chrY') == 'sex'
    # non-target contigs excluded from both outputs
    assert classify_contig('chrM') == 'exclude'
    assert classify_contig('chr1_KI270706v1_random') == 'exclude'
    assert classify_contig('chrUn_GL000195v1') == 'exclude'
    assert classify_contig('HLA-A*01:01:01:01') == 'exclude'


def test_validate_genders_ok():
    # all samples have a valid gender -> no error
    validate_genders({'s1': 'a.cram', 's2': 'b.cram'}, {'s1': 'male', 's2': 'female'})


def test_validate_genders_missing():
    with pytest.raises(ValueError):
        validate_genders({'s1': 'a.cram'}, {})


def test_validate_genders_invalid():
    with pytest.raises(ValueError):
        validate_genders({'s1': 'a.cram'}, {'s1': 'unknown'})


def test_mean_single():
    v1 = np.array([0, 10, 20], dtype=np.int64)
    v2 = np.array([2, 12, 24], dtype=np.int64)
    mean, n = mean_single([v1, v2])
    assert n == 2
    assert list(mean) == [1.0, 11.0, 22.0]


def test_mean_per_gender():
    items = [
        ('male', np.array([10, 0], dtype=np.int64)),
        ('male', np.array([20, 4], dtype=np.int64)),
        ('female', np.array([30, 8], dtype=np.int64)),
    ]
    mm, nm, mf, nf = mean_per_gender(items)
    assert nm == 2
    assert nf == 1
    assert list(mm) == [15.0, 2.0]
    assert list(mf) == [30.0, 8.0]


def test_merge():
    chrom2coordinate = dict()
    chrom2coordinate['chr1'] = list()
    chrom2coordinate['chr1'].append(99)
    chrom2coordinate['chr1'].append(100)
    chrom2coordinate['chr1'].append(110)
    chrom2coordinate['chr1'].append(200)

    chrom2depth = dict()
    chrom2depth['chr1'] = list()
    chrom2depth['chr1'].append({
        'chrom': 'chr1',
        'start': 100,
        'end': 110,
        'depth': 10
    })
    chrom2depth['chr1'].append({
        'chrom': 'chr1',
        'start': 190,
        'end': 200,
        'depth': 20
    })

    result = pl.from_pandas(merge(chrom2coordinate, chrom2depth))

    assert len(result) == 3

    record0 = result.row(0, named=True)
    record1 = result.row(1, named=True)
    record2 = result.row(2, named=True)

    assert record0['chrom'] == 'chr1'
    assert record0['pos'] == 100
    assert record0['depth'] == 10

    assert record1['chrom'] == 'chr1'
    assert record1['pos'] == 110
    assert record1['depth'] == 10

    assert record2['chrom'] == 'chr1'
    assert record2['pos'] == 200
    assert record2['depth'] == 20
