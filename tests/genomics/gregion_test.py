from pathlib import Path
import polars as pl
from importlib import resources
import pytest

from genomics.gregion import GenomicRegion
from genomics.gregion import create_database

def test_intersects():
    assert GenomicRegion('chr1', 100, 200).intersects( GenomicRegion('chr2', 100, 200)) is None
    assert GenomicRegion('chr1', 100, 200) == GenomicRegion('chr1', 100, 200).intersects( GenomicRegion('chr1', 100, 200))

    assert GenomicRegion('chr1', 150, 200) == GenomicRegion('chr1', 150, 200).intersects( GenomicRegion('chr1', 100, 200))
    assert GenomicRegion('chr1', 150, 200) == GenomicRegion('chr1', 150, 200).intersects( GenomicRegion('chr1', 100, 201))

    assert GenomicRegion('chr14', 21385192,22829820) == GenomicRegion('chr14', 18967499, 23423741).intersects( GenomicRegion('chr14', 21385192, 22829820))



# def test_GenomicRegionDatabase_complex_regions():
#     complex_regions_file = (
#         resources.files("genomics") / "resources" / "complex_regions-hg38.tsv"
#     )
#
#     db = create_database(pl.read_csv(complex_regions_file, has_header = True, separator = '\t').to_dicts())
#
#     matches = db.find_overlap({'chrom': 'chr6', 'start': 32478823, 'end': 32605693})
#
#     print('####')
#     print(matches)
#     print('####')
# {'sample_name': '20230615_GT1_ARRAY_2247_G02_TPM1.CEL', 'chrom': 'chr6', 'start': 32478823, 'end': 32605693, 'cn_state': '3', 'n_markers': 555, 'barcode': '5509754473359122124058', 'region_idx': 736}
# [{'chrom': 'chr6', 'start': 28510119, 'end': 29972223, 'category': 'MHC;HYPERVARIABLE', 'note': ''}]
#
#     data = pl.from_dict({
#         'chrom': ['chr6'],
#         'start': [28510119,],
#         'end': [29972223],
#         'name': ['r1', 'r2', 'r3'],
#     }).to_dicts()
#
#     regions = create_database(data)


def test_GenomicRegionDatabase():
    data = pl.from_dict({
        'chrom': ['chr1', 'chr1', 'chr2'],
        'start': [10000, 30000, 30000],
        'end': [20000, 40000, 40000],
        'name': ['r1', 'r2', 'r3'],
    }).to_dicts()

    regions = create_database(data)

    assert {'chr1', 'chr2'} == regions.chroms

    print(regions.find_overlap({'chrom': 'chr2', 'start': 40001, 'end': 40002}))
    assert len(regions.find_overlap({'chrom': 'chr2', 'start': 40001, 'end': 40002})) == 0
    assert len(regions.find_overlap({'chrom': 'chr3', 'start': 40001, 'end': 40002})) == 0

    result =  regions.find_overlap({'chrom':'chr1', 'start': 10000, 'end': 10001})
    assert len(result) == 1
    assert result[0]['chrom'] == 'chr1'
    assert result[0]['start'] == 10000
    assert result[0]['end'] == 20000
    assert result[0]['name'] == 'r1'

    result =  regions.find_overlap({'chrom': 'chr1', 'start': 20000, 'end': 30001})
    assert len(result) == 1
    assert result[0]['chrom'] == 'chr1'
    assert result[0]['start'] == 30000
    assert result[0]['end'] == 40000
    assert result[0]['name'] == 'r2'





def test_overlaps():
    assert not GenomicRegion('chr1', 100, 200).overlaps( GenomicRegion('chr2', 100, 200))
    assert not GenomicRegion('chr1', 100, 200).overlaps( GenomicRegion('chr1', 200, 300))
    assert GenomicRegion('chr1', 100, 101).overlaps(GenomicRegion('chr1', 100, 101))
    assert GenomicRegion('chr1', 100, 200).overlaps(GenomicRegion('chr1', 150, 160))
    assert GenomicRegion('chr1', 100, 200).overlaps(GenomicRegion('chr1', 150, 250))
    assert not GenomicRegion('chr1', 100, 200).overlaps(GenomicRegion('chr1', 50, 100))
    assert GenomicRegion('chr1', 100, 200).overlaps(GenomicRegion('chr1', 50, 101))
    assert not GenomicRegion('chr1', 100, 200).overlaps(GenomicRegion('chr1', 201, 300))


def test_merge():
    assert GenomicRegion('chr1', 100, 160) == \
        GenomicRegion('chr1', 100, 150).merge(GenomicRegion('chr1', 110, 160))

    with pytest.raises(Exception):
        GenomicRegion('chr1', 100, 150).merge(GenomicRegion('chr2', 110, 160))


def test_get_reciprocal_overlap():
    x = GenomicRegion('chr1', 100, 200)
    y = GenomicRegion('chr1', 100, 200)
    assert 1.0 == x.get_reciprocal_overlap(y)['ratio']

    x = GenomicRegion('chr1', 150, 200)
    y = GenomicRegion('chr1', 100, 200)
    assert 1.0 == x.get_reciprocal_overlap(y)['ratio']

    x = GenomicRegion('chr1', 150, 250)
    y = GenomicRegion('chr1', 100, 200)
    assert 0.5 == x.get_reciprocal_overlap(y)['ratio']

    x = GenomicRegion('chr1', 150, 200)
    y = GenomicRegion('chr1', 100, 150)
    assert 0.0 == x.get_reciprocal_overlap(y)['ratio']
