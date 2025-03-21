import pytest

from genomics.gregion import GenomicRegion


def test_overlaps():
    assert not GenomicRegion('chr1', 100, 200).overlaps(
        GenomicRegion('chr2', 100, 200))
    assert not GenomicRegion('chr1', 100, 200).overlaps(
        GenomicRegion('chr1', 201, 300))
    assert GenomicRegion('chr1', 100,
                         100).overlaps(GenomicRegion('chr1', 100, 100))
    assert GenomicRegion('chr1', 100,
                         200).overlaps(GenomicRegion('chr1', 150, 160))
    assert GenomicRegion('chr1', 100,
                         200).overlaps(GenomicRegion('chr1', 150, 250))


def test_merge():
    assert GenomicRegion('chr1', 100,
                         160) == GenomicRegion('chr1', 100, 150).merge(
                             GenomicRegion('chr1', 110, 160))
    with pytest.raises(Exception):
        GenomicRegion('chr1', 100, 150).merge(GenomicRegion('chr2', 110, 160))
