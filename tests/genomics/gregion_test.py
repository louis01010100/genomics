from pathlib import Path
import polars as pl
import pytest

from genomics.gregion import GenomicRegion, create_genomic_regions


def test_GenomicRegions():
    data = pl.from_dict({
        'chrom': ['chr1', 'chr1', 'chr2'],
        'start': [10000, 30000, 30000],
        'end': [20000, 40000, 40000],
        'name': ['r1', 'r2', 'r3'],
    })

    regions = create_genomic_regions(data)

    assert {'chr1', 'chr2'} == regions.chroms

    assert len(regions.find_overlap('chr2', 40001, 40002)) == 0
    assert len(regions.find_overlap('chr3', 40001, 40002)) == 0

    result =  regions.find_overlap('chr1', 10000, 10001)
    assert len(result) == 1
    assert result[0] == GenomicRegion('chr1', 10000, 20000, 'r1')

    result =  regions.find_overlap('chr1', 20000, 30001)
    assert len(result) == 1
    assert result[0] == GenomicRegion('chr1', 30000, 40000, 'r2')


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
