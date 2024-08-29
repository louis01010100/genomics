from genomics.variant import is_vcf, Variant, _load_allele2idx, _load_idx2allele, transcode_gt, _transcode_gt, normalize, align, denormalize, normalize_chrom_name, sync
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
    assert is_vcf( 'A', 'C')
    assert is_vcf( 'AC', 'A')
    assert is_vcf( 'A', 'AC')
    assert is_vcf( 'AT', 'CG')
    assert is_vcf( 'AT', 'ACCG')

    assert not is_vcf( 'A', '-')
    assert not is_vcf( '-', 'A')


def test_to_vcf():
    genome_file = Path(__file__).parents[0] / 'seq.fa'
    genome = Genome(genome_file)

    v = Variant(chrom = 'chr1', pos = 2, ref = 'C', alt = 'G')
    assert v == v.to_vcf(genome)

    assert Variant(chrom = 'chr1', pos = 2, ref = 'C', alt = 'CG') \
            == Variant(chrom = 'chr1', pos = 2, ref = '-', alt = 'G').to_vcf(genome)

    assert Variant(chrom = 'chr1', pos = 1, ref = 'AC', alt = 'A') \
            == Variant(chrom = 'chr1', pos = 2, ref = 'C', alt = '-').to_vcf(genome)

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
    genome_file = Path(__file__).parents[0] / 'seq.fa'
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

    v = align(
        chrom = 'chr2', 
        pos = 12, 
        ref = 'T', 
        alts = ['A'], 
        ref_pos = 10, 
        ref_ref = 'TTTCA', 
        ref_alts = ['GAATGATC','TTACA','TTTA','TTTTA'],
        genome = genome,
    )

    assert 'chr2' == v['chrom']
    assert 10 == v['pos']
    assert 'TTTCA' == v['ref']
    assert ['TTACA']  == v['alts']


    v = align(
        chrom = 'chr1', 
        pos = 4, ref = 'TA', alts = 'AG', 
        ref_pos = 6, ref_ref = 'AA', ref_alts = 'A,AAA', 
        genome = genome
    )
    assert None == v


def test_denormalize(tmp_path):
    genome_file = Path(__file__).parents[0] / 'seq.fa'
    genome = Genome(genome_file)

    # 123456789012
    # ACGTAAAAAAAT
    # A
    # T
    result = Variant('chr1', 1, 'A', 'T').denormalize(genome)
    assert 'chr1' == result.chrom
    assert 1 == result.pos
    assert 'A' == result.ref
    assert 'T' == result.alt

    # 123456789012
    # ACGTAAAAAAAT
    # A
    # T
    result = Variant('chr1', 1, 'AC', 'TG').denormalize(genome)
    assert 'chr1' == result.chrom
    assert 1 == result.pos
    assert 'AC' == result.ref
    assert 'TG' == result.alt

    # 123456789012
    # ACGTAAAAAAAT
    #    TA 
    #    T
    #          AA 
    #          A

    result = Variant('chr1', 4, 'TA', 'T').denormalize(genome)
    assert 'chr1' == result.chrom
    assert 10 == result.pos
    assert 'AA' == result.ref
    assert 'A' == result.alt

    # 123456789012
    # ACGTAAAAAAAT
    #    T 
    #    TA
    #           A 
    #           AA
    result = Variant('chr1', 4, 'T', 'TA').denormalize(genome)
    assert 'chr1' == result.chrom
    assert 11 == result.pos
    assert 'A' == result.ref
    assert 'AA' == result.alt

    # 123456789012
    # ACGTAAAAAAAT
    #    T 
    #    TA
    #    TAA
    #
    #           A
    #           AA
    #           AAA

    result = Variant('chr1', 4, 'T', 'TA,TAA').denormalize(genome)
    assert 'chr1' == result.chrom
    assert 11 == result.pos
    assert 'A' == result.ref
    assert 'AA,AAA' == result.alt


    # 12345
    # GTTGG
    # GT
    # G
    #  TT
    #  T
    result = Variant('chr3', 1, 'GT', 'G').denormalize(genome)
    assert 'chr3' == result.chrom
    assert  2== result.pos
    assert 'TT' == result.ref
    assert 'T' == result.alt
            

    # 1234567890
    # ATTCTTCTCG
    # ATTC
    # A
    #     TTCT
    #     T

    result = Variant('chr4', 1, 'ATTC', 'A').denormalize(genome)
    assert 'chr4' == result.chrom
    assert 'TTCT' == result.ref
    assert 'T' == result.alt
    assert  5== result.pos


def test_normalize(tmp_path):
    genome_file = Path(__file__).parents[0] / 'seq.fa'
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

    genome_file = Path(__file__).parents[0] / 'seq.fa'
    genome = Genome(genome_file)

    ref_variant = Variant(chrom = 'chr1', pos = 6, ref = 'AA', alt = 'A,AAA')
    v = Variant(chrom = 'chr1', pos = 4, ref = 'TA', alt = 'TAA', format_ = 'GT', calls = '0/0\t0/1\t1/1')
    v1 = v.align(ref_variant, genome)

    assert 'chr1' == v1.chrom 
    assert 6 == v1.pos 
    assert 'AA' == v1.ref 
    assert 'A,AAA' == v1.alt
    assert '0/0\t0/2\t2/2' == v1.calls

# chr13   32394895        AX-599399085    TTTCA   GAATGATC,TTACA,TTTA,TTTTA
    ref_variant = Variant(chrom = 'chr2', pos = 10, ref = 'TTTCA', alt = 'GAATGATC,TTACA,TTTA,TTTTA')

    v = Variant(chrom = 'chr2', pos = 12, ref = 'T', alt = 'A', format_ = 'GT', calls = '1/1\t0/1\t0/0')

    v1 = v.align(ref_variant, genome)

    assert 'chr2' == v1.chrom 
    assert 10 == v1.pos 
    assert 'TTTCA' == v1.ref 
    assert 'GAATGATC,TTACA,TTTA,TTTTA' == v1.alt
    assert '2/2\t0/2\t0/0' == v1.calls

def test_normalize_chrom_name():
    assert 'chr1' == normalize_chrom_name(1)
    assert 'chrX' == normalize_chrom_name('x')
    assert 'chrY' == normalize_chrom_name('Y')
    assert 'chrM' == normalize_chrom_name('MT')
    assert 'chr1' == normalize_chrom_name('chr1')


def test_sync():
    genome_file = Path(__file__).parents[0] / 'seq.fa'
    genome = Genome(genome_file)

    vx = Variant(chrom = 'chr1', pos = 3, ref = 'G', alt = 'C')
    vy = Variant(chrom = 'chr1', pos = 3, ref = 'G', alt = 'T')
    nvx, nvy = sync(vx,vy, genome)
    assert nvx.pos == nvy.pos == 3
    assert nvx.ref == nvy.ref == 'G'
    assert nvx.alt == 'C'
    assert nvy.alt == 'T'

    vx = Variant(chrom = 'chr1', pos = 3, ref = 'GC', alt = 'CC')
    vy = Variant(chrom = 'chr1', pos = 3, ref = 'G', alt = 'T')
    nvx, nvy = sync(vx,vy, genome)
    assert nvx.pos == nvy.pos == 3
    assert nvx.ref == nvy.ref == 'GC'
    assert nvx.alt == 'CC'
    assert nvy.alt == 'TC'

    vx = Variant(chrom = 'chr1', pos = 3, ref = 'G', alt = 'C')
    vy = Variant(chrom = 'chr1', pos = 3, ref = 'GC', alt = 'TT')
    nvx, nvy = sync(vx,vy, genome)
    assert 3 == nvx.pos == nvy.pos
    assert 'GC' == nvx.ref == nvy.ref
    assert 'CC' == nvx.alt
    assert 'TT' == nvy.alt

    vx = Variant(chrom = 'chr1', pos = 3, ref = 'G', alt = 'C')
    vy = Variant(chrom = 'chr1', pos = 4, ref = 'T', alt = 'G')
    nvx, nvy = sync(vx,vy, genome)
    assert 3 == nvx.pos == nvy.pos
    assert 'GT' == nvx.ref
    assert 'GT' == nvy.ref
    assert 'CT' == nvx.alt
    assert 'GG' == nvy.alt

    vx = Variant(chrom = 'chr1', pos = 3, ref = 'G', alt = 'C')
    vy = Variant(chrom = 'chr1', pos = 2, ref = 'C', alt = 'A')
    nvx, nvy = sync(vx,vy, genome)
    assert 2 == nvx.pos == nvy.pos
    assert 'CG' == nvx.ref == nvy.ref
    assert 'CC' == nvx.alt
    assert 'AG' == nvy.alt

    vx = Variant(chrom = 'chr1', pos = 3, ref = 'G', alt = 'A')
    vy = Variant(chrom = 'chr1', pos = 2, ref = 'CG', alt = 'C')
    nvx, nvy = sync(vx,vy, genome)
    assert 2 == nvx.pos == nvy.pos
    assert 'CG' == nvx.ref == nvy.ref
    assert 'CA' == nvx.alt
    assert 'C' == nvy.alt

    vx = Variant(chrom = 'chr2', pos = 3, ref = 'A', alt = 'AA,AAC')
    vy = Variant(chrom = 'chr2', pos = 2, ref = 'A', alt = 'C')
    nvx, nvy = sync(vx,vy, genome)
    assert 2 == nvx.pos == nvy.pos
    assert 'AA' == nvx.ref == nvy.ref
    assert 'AAA,AAAC' == nvx.alt
    assert 'CA' == nvy.alt
