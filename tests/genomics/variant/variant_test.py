from genomics.variant import is_vcf, Variant, normalize, denormalize, normalize_chrom_name, sync, get_max_region
from genomics.gregion import GenomicRegion
from genomics.genome import Genome
from pathlib import Path


def test_region():

    assert Variant(
        'chr1',
        10013,
        'TA',
        'T',
    ).region == GenomicRegion('chr1', 10012, 10014)

    assert Variant(
        'chr1',
        100,
        'AC',
        'A',
    ).region == GenomicRegion('chr1', 99, 101)

    assert Variant(
        'chr1',
        101,
        'C',
        'G',
    ).region == GenomicRegion('chr1', 100, 101)

    assert Variant(
        'chr1',
        100,
        'CT',
        'AT,C',
    ).region == GenomicRegion('chr1', 99, 101)

    assert Variant(
        'chr1',
        100,
        'T',
        'TG',
    ).region == GenomicRegion('chr1', 99, 100)

    assert Variant(
        'chr13',
        32332376,
        'GTA',
        'GTAG,GTG,TT',
    ).region == GenomicRegion('chr13', 32332375, 32332378)

    assert Variant(
        'chr13',
        48037782,
        'AGGAGTC',
        'AGGAGTCGGAGTC',
    ).region == GenomicRegion('chr13', 48037781, 48037788)


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
            == Variant(chrom = 'chr1', pos = 2, ref = '-', alt = 'G').to_vcf(genome.chromosome('chr1'))

    assert Variant(chrom = 'chr1', pos = 1, ref = 'AC', alt = 'A') \
            == Variant(chrom = 'chr1', pos = 2, ref = 'C', alt = '-').to_vcf(genome.chromosome('chr1'))

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



def test_get_max_region(tmp_path):
    # 123456789012
    # TCAGAGAAA
    #  CAG 
    #  C
    #     AGA
    #     A
    genome_file = Path(__file__).parents[0] / 'seq.fa'
    genome = Genome(genome_file)
    result = Variant('chr5', 2, 'CAG', 'C')

    region = get_max_region(result, genome.chromosome('chr5'))

    assert 1 == region.start 
    assert 7 == region.end

    # 12345678901
    # CCCAAGACGTT
    #   CAAGA
    #   CGT
    genome_file = Path(__file__).parents[0] / 'seq.fa'
    genome = Genome(genome_file)
    result = Variant('chr6', 3, 'CAAGA', 'CGT')

    region = get_max_region(result, genome.chromosome('chr5'))

    assert 2 == region.start 
    assert 7 == region.end

def test_denormalize(tmp_path):
    genome_file = Path(__file__).parents[0] / 'seq.fa'
    genome = Genome(genome_file)

    # 123456789012
    # ACGTAAAAAAAT
    # A
    # T
    result = Variant('chr1', 1, 'A', 'T').denormalize(genome.chromosome('chr1'))
    assert 'chr1' == result.chrom
    assert 1 == result.pos
    assert 'A' == result.ref
    assert 'T' == result.alt

    # 123456789012
    # ACGTAAAAAAAT
    # A
    # T
    result = Variant('chr1', 1, 'AC', 'TG').denormalize(genome.chromosome('chr1'))
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

    result = Variant('chr1', 4, 'TA', 'T').denormalize(genome.chromosome('chr1'))
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
    result = Variant('chr1', 4, 'T', 'TA').denormalize(genome.chromosome('chr1'))
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

    result = Variant('chr1', 4, 'T', 'TA,TAA').denormalize(genome.chromosome('chr1'))
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
    result = Variant('chr3', 1, 'GT', 'G').denormalize(genome.chromosome('chr3'))
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

    result = Variant('chr4', 1, 'ATTC', 'A').denormalize(genome.chromosome('chr4'))
    assert 'chr4' == result.chrom
    assert 'TTCT' == result.ref
    assert 'T' == result.alt
    assert  5== result.pos


def test_normalize(tmp_path):
    genome_file = Path(__file__).parents[0] / 'seq.fa'
    genome = Genome(genome_file)

    result = Variant(chrom = 'chr1', pos = 2, ref = 'C', alt = 'G').normalize(genome.chromosome('chr1'))
    assert 'chr1' == result.chrom
    assert 2 == result.pos
    assert 'C' == result.ref
    assert 'G' == result.alt

    result = Variant(chrom = 'chr1', pos = 5, ref = 'A', alt = 'AA').normalize(genome.chromosome('chr1'))

    assert 'chr1' == result.chrom
    assert 4 == result.pos
    assert 'T' == result.ref
    assert 'TA' == result.alt


    result = Variant(chrom = 'chr1', pos = 6, ref = 'AA', alt = 'A').normalize(genome.chromosome('chr1'))

    assert 'chr1' == result.chrom
    assert 4 == result.pos
    assert 'TA' == result.ref
    assert 'T' == result.alt

    result = Variant(chrom = 'chr1', pos = 6, ref = 'AA', alt = 'A,AAA').normalize(genome.chromosome('chr1'))

    assert 'chr1' == result.chrom
    assert 4 == result.pos
    assert 'TA' == result.ref
    assert 'T,TAA' == result.alt

    result = Variant(chrom = 'chr1', pos = 6, ref = 'AA', alt = 'AAA,A').normalize(genome.chromosome('chr1'))

    assert 'chr1' == result.chrom
    assert 4 == result.pos
    assert 'TA' == result.ref
    assert 'TAA,T' == result.alt


def test_normalize_chrom_name():
    assert 'chr1' == normalize_chrom_name(1)
    assert 'chrX' == normalize_chrom_name('x')
    assert 'chrY' == normalize_chrom_name('Y')
    assert 'chrM' == normalize_chrom_name('MT')
    assert 'chr1' == normalize_chrom_name('chr1')

def test_expand():
    genome_file = Path(__file__).parents[0] / 'seq.fa'
    genome = Genome(genome_file)
    v = Variant(chrom = 'chr2', pos = 3, ref = 'A', alt = 'C').expand(genome.chromosome('chr2'))
    assert 'chr2' == v.chrom
    assert 3 == v.pos
    assert 'A' == v.ref
    assert 'C' == v.alt

    v = Variant(chrom = 'chr2', pos = 3, ref = 'A', alt = 'CC').expand(genome.chromosome('chr2'))
    assert 'chr2' == v.chrom
    assert 3 == v.pos
    assert 'A' == v.ref
    assert 'CC' == v.alt

    v = Variant(chrom = 'chr2', pos = 3, ref = 'A', alt = 'AA').expand(genome.chromosome('chr2'))
    assert 'chr2' == v.chrom
    assert 2 == v.pos
    assert 'GAAAA' == v.ref
    assert 'GAAAAA' == v.alt

    genome_file = Path(__file__).parents[0] / 'seq.fa'
    genome = Genome(genome_file)
    v = Variant(chrom = 'chr2', pos = 4, ref = 'A', alt = 'AA,AAC')

    v = v.expand(genome.chromosome('chr2'))
    assert 'chr2' == v.chrom
    assert 2 == v.pos
    assert 'GAAAA' == v.ref
    assert 'GAAAAA,GAAACAA' == v.alt


def test_sync():
    genome_file = Path(__file__).parents[0] / 'seq.fa'
    genome = Genome(genome_file)

    # vx = Variant(chrom = 'chr1', pos = 3, ref = 'G', alt = 'C')
    # vy = Variant(chrom = 'chr1', pos = 3, ref = 'G', alt = 'T')
    # nvx, nvy = sync(vx,vy, genome.chromosome('chr1'))
    # assert nvx.pos == nvy.pos == 3
    # assert nvx.ref == nvy.ref == 'G'
    # assert nvx.alt == 'C'
    # assert nvy.alt == 'T'
    #
    # vx = Variant(chrom = 'chr1', pos = 3, ref = 'GT', alt = 'CC')
    # vy = Variant(chrom = 'chr1', pos = 3, ref = 'G', alt = 'T')
    # nvx, nvy = sync(vx,vy, genome.chromosome('chr1'))
    # assert 3 == nvx.pos == nvy.pos
    # assert 'GT' == nvx.ref == nvy.ref
    # assert 'CC' == nvx.alt
    # assert 'TT' == nvy.alt
    #
    #
    # vx = Variant(chrom = 'chr1', pos = 3, ref = 'G', alt = 'C')
    # vy = Variant(chrom = 'chr1', pos = 3, ref = 'GT', alt = 'AA')
    # nvx, nvy = sync(vx,vy, genome.chromosome('chr1'))
    # assert 3 == nvx.pos == nvy.pos
    # assert 'GT' == nvx.ref == nvy.ref
    # assert 'CT' == nvx.alt
    # assert 'AA' == nvy.alt

    vx = Variant(chrom = 'chr1', pos = 3, ref = 'G', alt = 'C')
    vy = Variant(chrom = 'chr1', pos = 4, ref = 'T', alt = 'G')
    nvx, nvy = sync(vx,vy, genome.chromosome('chr1'))
    assert 3 == nvx.pos == nvy.pos
    assert 'GT' == nvx.ref
    assert 'GT' == nvy.ref
    assert 'CT' == nvx.alt
    assert 'GG' == nvy.alt

    # vx = Variant(chrom = 'chr1', pos = 3, ref = 'G', alt = 'C')
    # vy = Variant(chrom = 'chr1', pos = 2, ref = 'C', alt = 'A')
    # nvx, nvy = sync(vx,vy, genome.chromosome('chr1'))
    # assert 2 == nvx.pos == nvy.pos
    # assert 'CG' == nvx.ref == nvy.ref
    # assert 'CC' == nvx.alt
    # assert 'AG' == nvy.alt
    #
    # vx = Variant(chrom = 'chr1', pos = 3, ref = 'G', alt = 'A')
    # vy = Variant(chrom = 'chr1', pos = 2, ref = 'CG', alt = 'C')
    # nvx, nvy = sync(vx,vy, genome.chromosome('chr1'))
    # assert 2 == nvx.pos == nvy.pos
    # assert 'CG' == nvx.ref == nvy.ref
    # assert 'CA' == nvx.alt
    # assert 'C' == nvy.alt
    #
    # vx = Variant(chrom = 'chr2', pos = 4, ref = 'A', alt = 'AA,AAC')
    # vy = Variant(chrom = 'chr2', pos = 3, ref = 'A', alt = 'C')
    #
    # nvx, nvy = sync(vx,vy, genome.chromosome('chr2'))
    # assert 2 == nvx.pos == nvy.pos
    # assert 'GAAAA' == nvx.ref == nvy.ref
    # assert 'GAAAAA,GAAACAA' == nvx.alt
    # assert 'GCAAA' == nvy.alt
