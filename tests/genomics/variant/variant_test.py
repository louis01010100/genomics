from genomics.variant import is_vcf, Variant, _load_allele2idx, _load_idx2allele, transcode_gt, _transcode_gt, normalize, align
from genomics.gregion import GenomicRegion
from genomics.genome import Genome
from pathlib import Path


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



def test__load_allele2idx():
    assert _load_allele2idx('A', ['C']) == {'A': '0', 'C': '1', '.': '.'}
    assert _load_allele2idx('C', ['CT','G']) == {
        'C': '0',
        'CT': '1',
        'G': '2',
        '.': '.'
    }


def test__load_idx2allele():
    assert _load_idx2allele('A', ['C']) == {'0': 'A', '1': 'C', '.': '.'}
    assert _load_idx2allele('C', ['CT','G']) == {
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

def test_align():
    genome_file = Path(__file__).parents[0] / 'seq.fa.bgz'
    genome = Genome(genome_file)

    v = align(
        chrom = 'chr1', 
        pos = 4, ref = 'TA', alts = ['T'], 
        ref_pos = 4, ref_ref = 'TA', ref_alts = ['T'], 
        genome = genome
    )

    assert 'chr1' == v['chrom']
    assert 4 == v['pos']
    assert 'TA' == v['ref']
    assert ['T'] == v['alts']

    v = align(
        chrom = 'chr1', 
        pos = 4, ref = 'TA', alts = ['T','TAA'], 
        ref_pos = 6, ref_ref = 'AA', ref_alts = ['A','AAA'], 
        genome = genome
    )

    assert 'chr1' == v['chrom']
    assert 6 == v['pos']
    assert 'AA' == v['ref']
    assert ['A','AAA'] == v['alts']


    try:
        v = align(
            chrom = 'chr1', 
            pos = 4, ref = 'TA', alts = 'AG', 
            ref_pos = 6, ref_ref = 'AA', ref_alts = 'A,AAA', 
            genome = genome
        )
        assert False
    except Exception as e:
        assert True



def test_normalize(tmp_path):
    genome_file = Path(__file__).parents[0] / 'seq.fa.bgz'
    genome = Genome(genome_file)


    result = Variant(chrom = 'chr1', pos = 2, ref = 'C', alt = 'G').normalize(genome)
    assert 'chr1' == result.chrom
    assert 2 == result.pos
    assert 'C' == result.ref
    assert 'G' == result.alt

    result = Variant(chrom = 'chr1', pos = 5, ref = 'A', alt = 'AA').normalize(genome)

    assert 'chr1' == result.chrom
    assert 4 == result.pos
    assert 'T' == result.ref
    assert 'TA' == result.alt


    result = Variant(chrom = 'chr1', pos = 6, ref = 'AA', alt = 'A').normalize(genome)

    assert 'chr1' == result.chrom
    assert 4 == result.pos
    assert 'TA' == result.ref
    assert 'T' == result.alt

    result = Variant(chrom = 'chr1', pos = 6, ref = 'AA', alt = 'A,AAA').normalize(genome)

    assert 'chr1' == result.chrom
    assert 4 == result.pos
    assert 'TA' == result.ref
    assert 'T,TAA' == result.alt

    result = Variant(chrom = 'chr1', pos = 6, ref = 'AA', alt = 'AAA,A').normalize(genome)

    assert 'chr1' == result.chrom
    assert 4 == result.pos
    assert 'TA' == result.ref
    assert 'TAA,T' == result.alt

def test_variant_align():

    genome_file = Path(__file__).parents[0] / 'seq.fa.bgz'
    genome = Genome(genome_file)

    ref_variant = Variant(chrom = 'chr1', pos = 6, ref = 'AA', alt = 'A,AAA')
    v = Variant(chrom = 'chr1', pos = 4, ref = 'TA', alt = 'TAA', format_ = 'GT', calls = '0/0\t0/1\t1/1')
    v1 = v.align(ref_variant, genome)

    assert 'chr1' == v1.chrom 
    assert 6 == v1.pos 
    assert 'AA' == v1.ref 
    assert 'A,AAA' == v1.alt
    assert '0/0\t0/2\t2/2' == v1.calls

