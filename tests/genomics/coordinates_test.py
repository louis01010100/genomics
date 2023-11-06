from genomics.coordinate import merge
from genomics.variant import Variant


def test_merge():

    records = [{
        'chrom': 'chr1',
        'pos': 784860,
        'id': 'AX-32104085',
        'ref': 'T',
        'alt': 'C'
    }]

    result = merge(records)

    assert result == Variant(
        chrom='chr1',
        pos=784860,
        id_='AX-32104085',
        ref='T',
        alt='C',
    )
