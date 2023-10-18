from genomics.variant import is_vcf, Variant


def test_is_vcf():
    assert is_vcf(100, 100, 'A', 'C')
    assert is_vcf(100, 100, 'AC', 'A')
    assert is_vcf(100, 100, 'A', 'AC')
    assert is_vcf(100, 100, 'AT', 'CG')
    assert is_vcf(100, 100, 'AT', 'ACCG')

    assert not is_vcf(100, 100, 'A', '-')
    assert not is_vcf(100, 100, '-', 'A')


def test_merge__identical():
    v1 = Variant('chr1', 100, 'A', 'C')
    v2 = Variant('chr1', 100, 'A', 'C')

    result = v1.merge(v2)
    assert result == Variant('chr1', 100, 'A', 'C')


def test_merge__diff_alt():
    v1 = Variant('chr1', 100, 'A', 'C')
    v2 = Variant('chr1', 100, 'A', 'G')

    result = v1.merge(v2)
    assert result.chrom == Variant('chr1', 100, 'A', 'C,G')
