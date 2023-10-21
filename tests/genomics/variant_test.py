from genomics.variant import is_vcf, Variant, sync_prefix, sync_suffix
from genomics.gregion import GenomicRegion


def test_region():

    assert Variant(
        'chr1',
        100,
        'AC',
        'A',
    ).region == GenomicRegion('chr1', 100, 101)

    assert Variant(
        'chr1',
        101,
        'C',
        'G',
    ).region == GenomicRegion('chr1', 101, 101)

    assert Variant(
        'chr1',
        100,
        'CT',
        'AT,C',
    ).region == GenomicRegion('chr1', 100, 101)

    assert Variant(
        'chr1',
        100,
        'T',
        'TG',
    ).region == GenomicRegion('chr1', 100, 100)


def test_is_vcf():
    assert is_vcf(100, 100, 'A', 'C')
    assert is_vcf(100, 100, 'AC', 'A')
    assert is_vcf(100, 100, 'A', 'AC')
    assert is_vcf(100, 100, 'AT', 'CG')
    assert is_vcf(100, 100, 'AT', 'ACCG')

    assert not is_vcf(100, 100, 'A', '-')
    assert not is_vcf(100, 100, '-', 'A')


def test_sync_prefix():
    pos, ref1, alts1, ref2, alts2 = sync_prefix(
        start1=100,
        ref1='AC',
        alts1=['A'],
        start2=101,
        ref2='C',
        alts2=['G'],
    )

    assert pos == 100
    assert ref1 == 'AC'
    assert alts1 == ['A']
    assert ref2 == 'AC'
    assert alts2 == ['AG']

    pos, ref1, alts1, ref2, alts2 = sync_prefix(
        start1=101,
        ref1='C',
        alts1=['G'],
        start2=100,
        ref2='AC',
        alts2=['A'],
    )

    assert pos == 100
    assert ref1 == 'AC'
    assert alts1 == ['AG']
    assert ref2 == 'AC'
    assert alts2 == ['A']


def test_sync_suffix():
    # 100-102 AAA CCC
    # 100-101 AA AC
    ref1, alts1, ref2, alts2 = sync_suffix(
        end1=102,
        ref1='TAA',
        alts1=['CCC'],
        end2=101,
        ref2='TA',
        alts2=['GG'],
    )

    assert ref1 == 'TAA'
    assert alts1 == ['CCC']
    assert ref2 == 'TAA'
    assert alts2 == ['GGA']

    ref1, alts1, ref2, alts2 = sync_suffix(
        end1=100,
        ref1='C',
        alts1=['T'],
        end2=102,
        ref2='CTT',
        alts2=['C'],
    )

    assert ref1 == 'CTT'
    assert alts1 == ['TTT']
    assert ref2 == 'CTT'
    assert alts2 == ['C']

    ref1, alts1, ref2, alts2 = sync_suffix(
        end1=100,
        ref1='C',
        alts1=['A'],
        end2=101,
        ref2='CT',
        alts2=['AT,C'],
    )

    assert ref1 == 'CT'
    assert alts1 == ['AT']
    assert ref2 == 'CT'
    assert alts2 == ['AT,C']


def test_merge__identical():
    v1 = Variant('chr1', 100, 'A', 'C')
    v2 = Variant('chr1', 100, 'A', 'C')

    result = v1.merge(v2)
    assert result == Variant('chr1', 100, 'A', 'C')


def test_merge__diff_alt():
    v1 = Variant('chr1', 100, 'A', 'C')
    v2 = Variant('chr1', 100, 'A', 'G')
    v3 = Variant('chr1', 100, 'A', 'AT,G')

    assert v1.merge(v2) == Variant('chr1', 100, 'A', 'C,G')
    assert v1.merge(v2).merge(v3) == Variant('chr1', 100, 'A', 'AT,C,G')


def test_merge__diff_pos():
    v1 = Variant('chr1', 100, 'AC', 'A')
    v2 = Variant('chr1', 101, 'C', 'G')

    assert v1.merge(v2) == Variant('chr1', 100, 'AC', 'A,AG')

    v1 = Variant('chr15', 48474581, 'C', 'T')
    v2 = Variant('chr15', 48474581, 'CTT', 'C')

    assert v1.merge(v2) == Variant('chr15', 48474581, 'CTT', 'C,TTT')

    v1 = Variant('chr3', 37008823, 'C', 'A')
    v2 = Variant('chr3', 37008823, 'CT', 'AT,C')

    assert v1.merge(v2) == Variant('chr3', 37008823, 'CT', 'AT,C')
