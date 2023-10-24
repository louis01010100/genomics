from genomics.variant import is_vcf, Variant, sync_prefix, sync_suffix, _load_allele2idx, _load_idx2allele, _transcode_gt
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

    assert Variant(
        'chr13',
        32332376,
        'GTA',
        'GTAG,GTG,TT',
    ).region == GenomicRegion('chr13', 32332376, 32332378)


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


def test_sync_alleles__identical():
    v1 = Variant('chr1', 100, 'A', 'C')
    v2 = Variant('chr1', 100, 'A', 'C')

    result = v1.sync_alleles(v2)
    assert result == Variant('chr1', 100, 'A', 'C')


def test_sync_alleles__diff_alt():
    v1 = Variant('chr1', 100, 'A', 'C')
    v2 = Variant('chr1', 100, 'A', 'G')
    v3 = Variant('chr1', 100, 'A', 'AT,G')

    assert v1.sync_alleles(
        v2,
        site_only=True,
    ) == Variant('chr1', 100, 'A', 'C,G')

    assert v1.sync_alleles(v2, site_only=True).sync_alleles(
        v3, site_only=True) == Variant('chr1', 100, 'A', 'AT,C,G')

    v1 = Variant('chr15', 48474581, 'C', 'T')
    v2 = Variant('chr15', 48474581, 'CTT', 'C')

    assert v1.sync_alleles(
        v2,
        site_only=True,
    ) == Variant('chr15', 48474581, 'CTT', 'C,TTT')

    v1 = Variant('chr3', 37008823, 'C', 'A')
    v2 = Variant('chr3', 37008823, 'CT', 'AT,C')

    assert v1.sync_alleles(
        v2,
        site_only=True,
    ) == Variant('chr3', 37008823, 'CT', 'AT,C')


def test_sync_alleles__diff_pos():
    v1 = Variant('chr1', 100, 'AC', 'A')
    v2 = Variant('chr1', 101, 'C', 'G')

    assert v1.sync_alleles(
        v2,
        site_only=True,
    ) == Variant('chr1', 100, 'AC', 'A,AG')


def test_sync_alleles__with_gt():
    v1 = Variant('chr1', 100, 'AC', 'A')
    v2 = Variant(
        chrom='chr1',
        pos=101,
        id_='rs123',
        ref='C',
        alt='G',
        format_='GT',
        calls='\t'.join(['0/0', '0/1', '1|1', '0', '1', '.']),
    )

    assert v1.sync_alleles(v2) == Variant(
        chrom='chr1',
        pos=100,
        id_='rs123',
        ref='AC',
        alt='A,AG',
        format_='GT',
        calls='\t'.join(['0/0', '0/2', '2|2', '0', '2', '.']),
    )


def test__load_allele2idx():
    assert _load_allele2idx('A', 'C') == {'A': '0', 'C': '1', '.': '.'}
    assert _load_allele2idx('C', 'CT,G') == {
        'C': '0',
        'CT': '1',
        'G': '2',
        '.': '.'
    }


def test__load_idx2allele():
    assert _load_idx2allele('A', 'C') == {'0': 'A', '1': 'C', '.': '.'}
    assert _load_idx2allele('C', 'CT,G') == {
        '0': 'C',
        '1': 'CT',
        '2': 'G',
        '.': '.'
    }


def test__transcode_gt():

    idx2allele = {'0': 'A', '1': 'G', '.': '.'}
    allele2idx = {'A': '0', 'C': '1', 'G': '2', '.': '.'}

    assert _transcode_gt('0/0', idx2allele, allele2idx) == '0/0'
    assert _transcode_gt('0/1', idx2allele, allele2idx) == '0/2'
    assert _transcode_gt('.', idx2allele, allele2idx) == '.'
    assert _transcode_gt('0', idx2allele, allele2idx) == '0'
    assert _transcode_gt('./0', idx2allele, allele2idx) == '0/.'
    assert _transcode_gt('0/.', idx2allele, allele2idx) == '0/.'


def test___str__():

    expected = '\t'.join([
        'chr1', '100', 'rs123', 'AC', 'A,AG', '.', '.', '.', 'GT', '0/0',
        '0/2', '2|2', '0', '2', '.'
    ])
    actual = str(
        Variant(
            chrom='chr1',
            pos=100,
            id_='rs123',
            ref='AC',
            alt='A,AG',
            qual='.',
            filter_='.',
            info='.',
            format_='GT',
            calls='\t'.join(['0/0', '0/2', '2|2', '0', '2', '.']),
        ))

    assert expected == actual
