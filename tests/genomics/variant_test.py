from genomics.variant import is_vcf, Variant, sync_prefix, sync_suffix, _load_allele2idx, _load_idx2allele, transcode_gt, _transcode_gt, sync_alleles
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

    assert Variant(
        'chr13',
        48037782,
        'AGGAGTC',
        'AGGAGTCGGAGTC',
    ).region == GenomicRegion('chr13', 48037782, 48037788)


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


def test_sync_alleles__merge_id():
    v1 = Variant('chr13', 48037782, 'A', 'AGGAGTC', id_='AX-593989106')
    v2 = Variant('chr13', 48037782, 'A', 'AGGAGTC', id_='AX-593989107')
    v3 = Variant('chr13', 48037782, 'AGGAGTC', 'A', id_='AX-314679669')
    v4 = Variant('chr13', 48037782, 'AGGAGTC', 'A', id_='AX-593989108')
    v5 = Variant('chr13',
                 48037782,
                 'AGGAGTC',
                 'A,AGGAGTCGGAGTC',
                 id_='AX-593989104')
    v6 = Variant('chr13',
                 48037782,
                 'AGGAGTC',
                 'A,AGGAGTCGGAGTC',
                 id_='AX-593989105')

    result = v1
    for v in [v2, v3, v4, v5, v6]:
        result, _ = sync_alleles(result, v, merge_id=True)

    assert result == Variant(
        'chr13',
        48037782,
        'AGGAGTC',
        'AGGAGTCGGAGTC',
        id_=
        'AX-314679669,AX-593989106,AX-593989107,AX-593989108,AX-593989104,AX-593989105',
    )


def test_sync_alleles__merge_alt():
    v1 = Variant('chr1', 100, 'A', 'C')
    v2 = Variant('chr1', 100, 'A', 'G')

    new_v1, new_v2 = sync_alleles(v1, v2, merge_alt=False)
    assert new_v1 == Variant('chr1', 100, 'A', 'C')
    assert new_v2 == Variant('chr1', 100, 'A', 'G')

    new_v1, new_v2 = sync_alleles(v1, v2, merge_alt=True)
    assert new_v1 == Variant('chr1', 100, 'A', 'C,G')
    assert new_v2 == Variant('chr1', 100, 'A', 'C,G')

    v1 = Variant('chr15', 48474581, 'C', 'T')
    v2 = Variant('chr15', 48474581, 'CTT', 'C')
    new_v1, new_v2 = sync_alleles(v1, v2, merge_alt=False)
    assert new_v1 == Variant('chr15', 48474581, 'CTT', 'TTT')
    assert new_v2 == Variant('chr15', 48474581, 'CTT', 'C')

    new_v1, new_v2 = sync_alleles(v1, v2, merge_alt=True)
    assert new_v1 == Variant('chr15', 48474581, 'CTT', 'C,TTT')
    assert new_v2 == Variant('chr15', 48474581, 'CTT', 'C,TTT')


def test_sync_alleles__identical():
    v1 = Variant('chr1', 100, 'A', 'C')
    v2 = Variant('chr1', 100, 'A', 'C')

    new_v1, new_v2 = sync_alleles(v1, v2)
    assert new_v1 == Variant('chr1', 100, 'A', 'C')
    assert new_v2 == Variant('chr1', 100, 'A', 'C')


def test_sync_alleles__diff_alt():
    v1 = Variant('chr1', 100, 'A', 'C')
    v2 = Variant('chr1', 100, 'A', 'G')
    v3 = Variant('chr1', 100, 'A', 'AT,G')

    new_v1, new_v2 = sync_alleles(v1, v2)
    assert new_v1 == Variant('chr1', 100, 'A', 'C')
    assert new_v2 == Variant('chr1', 100, 'A', 'G')

    v1 = Variant('chr15', 48474581, 'C', 'T')
    v2 = Variant('chr15', 48474581, 'CTT', 'C')

    new_v1, new_v2 = sync_alleles(v1, v2)
    assert new_v1 == Variant('chr15', 48474581, 'CTT', 'TTT')
    assert new_v2 == Variant('chr15', 48474581, 'CTT', 'C')

    v1 = Variant('chr3', 37008823, 'C', 'A')
    v2 = Variant('chr3', 37008823, 'CT', 'AT,C')

    new_v1, new_v2 = sync_alleles(v1, v2)
    assert new_v1 == Variant('chr3', 37008823, 'CT', 'AT')
    assert new_v2 == Variant('chr3', 37008823, 'CT', 'AT,C')

    v1 = Variant('chr13', 48037782, 'A', 'AGGAGTC', id_='AX-593989106')
    v2 = Variant('chr13', 48037782, 'A', 'AGGAGTC', id_='AX-593989107')
    v3 = Variant('chr13', 48037782, 'AGGAGTC', 'A', id_='AX-314679669')
    v4 = Variant('chr13', 48037782, 'AGGAGTC', 'A', id_='AX-593989108')
    v5 = Variant('chr13',
                 48037782,
                 'AGGAGTC',
                 'A,AGGAGTCGGAGTC',
                 id_='AX-593989104')
    v6 = Variant('chr13',
                 48037782,
                 'AGGAGTC',
                 'A,AGGAGTCGGAGTC',
                 id_='AX-593989105')

    result = v1
    for v in [v2, v3, v4, v5, v6]:
        result, _ = sync_alleles(result, v, merge_alt=True, merge_id=True)

    assert result == Variant(
        'chr13',
        48037782,
        'AGGAGTC',
        'A,AGGAGTCGGAGTC',
        id_=
        'AX-314679669,AX-593989106,AX-593989107,AX-593989108,AX-593989104,AX-593989105',
    )


def test_sync_alleles__diff_pos():
    v1 = Variant('chr1', 100, 'AC', 'A')
    v2 = Variant('chr1', 101, 'C', 'G')

    new_v1, new_v2 = sync_alleles(v1, v2)
    assert new_v1 == Variant('chr1', 100, 'AC', 'A')
    assert new_v2 == Variant('chr1', 100, 'AC', 'AG')

    new_v1, new_v2 = sync_alleles(
        v1,
        v2,
        merge_alt=True,
    )
    assert new_v1 == Variant('chr1', 100, 'AC', 'A,AG')
    assert new_v2 == Variant('chr1', 100, 'AC', 'A,AG')


def test_sync_alleles__with_gt():

    v1 = Variant(
        chrom='chr1',
        pos=100,
        id_='rs123',
        ref='AC',
        alt='A',
        format_='GT',
        calls='\t'.join(['0/1', '1/1', '0|0', '0', '1', '.']),
    )

    v2 = Variant(
        chrom='chr1',
        pos=101,
        id_='rs456',
        ref='C',
        alt='G',
        format_='GT',
        calls='\t'.join(['0/0', '0/1', '1|1', '0', '1', '.']),
    )

    new_v1, new_v2 = sync_alleles(v1, v2)

    assert new_v1 == Variant(
        chrom='chr1',
        pos=100,
        id_='rs123',
        ref='AC',
        alt='A',
        format_='GT',
        calls='\t'.join(['0/1', '1/1', '0|0', '0', '1', '.']),
    )

    assert new_v2 == Variant(
        chrom='chr1',
        pos=100,
        id_='rs456',
        ref='AC',
        alt='AG',
        format_='GT',
        calls='\t'.join(['0/0', '0/1', '1|1', '0', '1', '.']),
    )

    new_v1, new_v2 = sync_alleles(v1, v2, merge_alt=True)

    assert new_v1 == Variant(
        chrom='chr1',
        pos=100,
        id_='rs123',
        ref='AC',
        alt='A,AG',
        format_='GT',
        calls='\t'.join(['0/1', '1/1', '0|0', '0', '1', '.']),
    )

    assert new_v2 == Variant(
        chrom='chr1',
        pos=100,
        id_='rs456',
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
    allele2allele = {'A': 'A', 'G': 'G', '.': '.'}

    assert _transcode_gt('0/0', idx2allele, allele2allele,
                         allele2idx) == '0/0'
    assert _transcode_gt('0/1', idx2allele, allele2allele,
                         allele2idx) == '0/2'
    assert _transcode_gt('.', idx2allele, allele2allele, allele2idx) == '.'
    assert _transcode_gt('0', idx2allele, allele2allele, allele2idx) == '0'
    assert _transcode_gt('./0', idx2allele, allele2allele,
                         allele2idx) == '0/.'
    assert _transcode_gt('0/.', idx2allele, allele2allele,
                         allele2idx) == '0/.'


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


def test_transcode_gt():

    calls = transcode_gt(
        idx2allele={
            '0': 'C',
            '1': 'G'
        },
        allele2idx={
            'AC': '0',
            'A': '1',
            'AG': '2',
        },
        allele2allele={
            'C': 'AC',
            'G': 'AG'
        },
        calls='\t'.join(['0/0', '0/1', '1|1', '0', '1', '.']),
    )

    assert calls == '0/0\t0/2\t2|2\t0\t2\t.'
