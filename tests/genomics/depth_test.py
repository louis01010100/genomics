from genomics.depth import merge


def test_merge():
    chrom2coordinate = dict()
    chrom2coordinate['chr1'] = list()
    chrom2coordinate['chr1'].append({'chrom': 'chr1', 'pos': 99})
    chrom2coordinate['chr1'].append({'chrom': 'chr1', 'pos': 100})
    chrom2coordinate['chr1'].append({'chrom': 'chr1', 'pos': 110})
    chrom2coordinate['chr1'].append({'chrom': 'chr1', 'pos': 200})

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

    result = pl.from_pandasd(merge(chrom2coordinate, chrom2depth))

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
