from pathlib import Path

from genomics.utils import _AllelePairs, vcf2dict


def test__AllelePairs():

    ap = _AllelePairs()
    ap.add_allele_pair('A', 'AGGAGTC')
    ap.add_allele_pair('AGGAGTC', 'A')
    ap.add_allele_pair('AGGAGTC', 'A,GGAGTCGGAGTC')

    assert {
        'ref': 'AGGAGTC',
        'alt': 'AGGAGTCGGAGTC'
    } == ap.get_updated_allele_pair('A', 'AGGAGTC')
    assert {
        'ref': 'AGGAGTC',
        'alt': 'A'
    } == ap.get_updated_allele_pair('AGGAGTC', 'A')
    assert {
        'ref': 'AGGAGTC',
        'alt': 'A,GGAGTCGGAGTC'
    } == ap.get_updated_allele_pair('AGGAGTC', 'A,GGAGTCGGAGTC')


def test_vcf2dict():
    fixture = Path(__file__).parents[0] / 'fixture.vcf'

    result = vcf2dict(fixture, fixture)

    allele_pairs = result[('chr13', '48037782')]
    assert {
        'ref': 'AGGAGTC',
        'alt': 'A'
    } == allele_pairs.get_updated_allele_pair('AGGAGTC', 'A')
    assert {
        'ref': 'AGGAGTC',
        'alt': 'A,AGGAGTCGGAGTC'
    } == allele_pairs.get_updated_allele_pair('AGGAGTC', 'A,AGGAGTCGGAGTC')
    assert {
        'ref': 'AGGAGTC',
        'alt': 'AGGAGTCGGAGTC'
    } == allele_pairs.get_updated_allele_pair('A', 'AGGAGTC')
