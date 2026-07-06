#!/usr/bin/env python
from pathlib import Path

import pandas as pd
import pytest
from genomics.gregion import GenomicRegion
from genomics.variant import Variant
from genomics.vcf import Vcf, concat, filter_variants, split_rtrim
import polars as pl


def test_meta(tmp_path):
    vcf_file = Path(__file__).parents[0] / 'fixture/sample.vcf'
    vcf = Vcf(vcf_file, tmp_path)
    expected = (
        '##fileformat=VCFv4.1\n'
        '##FILTER=<ID=PASS,Description="All filters passed">\n'
        '##contig=<ID=chr21>\n'
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Phased Genotype">')
    assert expected == vcf.meta


def test_contigs(tmp_path):
    vcf_file = Path(__file__).parents[0] / 'fixture/sample.vcf'
    vcf = Vcf(vcf_file, tmp_path).bgzip().index()

    contig_names = vcf.contigs

    assert {'chr21'} == contig_names


def test_header(tmp_path):
    vcf_file = Path(__file__).parents[0] / 'fixture/header.vcf'
    vcf = Vcf(vcf_file, tmp_path)
    expected = '#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	HG00403	HG00404	HG00405'
    assert expected == vcf.header


def test_bgzip(tmp_path):
    vcf_file = Path(__file__).parents[0] / 'fixture/sample.vcf'
    vcf = Vcf(vcf_file, tmp_path)
    vcf = vcf.bgzip()
    assert 'sample.vcf.bgz' == vcf.filepath.name
    assert vcf.filepath.exists()


def test_index(tmp_path):
    vcf_file = Path(__file__).parents[0] / 'fixture/sample.vcf.bgz'
    vcf = Vcf(vcf_file, tmp_path)
    vcf.index()

    index_file = vcf_file.with_suffix('.bgz.csi')
    assert index_file.exists()

    index_file.unlink()


def test_to_tsv(tmp_path):
    vcf_file = Path(__file__).parents[0] / 'fixture/sample.vcf.bgz'
    vcf = Vcf(vcf_file, tmp_path)
    output_file = vcf.to_tsv(
        '[%CHROM\t%POS\t%ID\t%REF\t%ALT\t%SAMPLE\t%GT\n]')

    data = pd.read_csv(output_file, header=0, sep='\t', dtype='str')

    record1 = data.iloc[0, :]

    assert record1['chrom'] == 'chr21'
    assert record1['pos'] == '5030578'
    assert record1['id'] == '.'
    assert record1['ref'] == 'C'
    assert record1['alt'] == 'T'
    assert record1['sample'] == 'HG00403'
    assert record1['gt'] == '0|0'


# def test_annotate_tag(tmp_path):
#     target_data = [
#         "##fileformat=VCFv4.1",
#         "##contig=<ID=chr21>",
#         "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO",
#         "chr21	1000000	AX-100	G	C	.	.	.",
#     ]
#     annotation_data = [
#         "##fileformat=VCFv4.1",
#         "##contig=<ID=chr21>",
#         "##INFO=<ID=AF,Number=A,Type=Float>",
#         "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO",
#         "chr21	1000000	.	G	C	.	.	AF=0.5",
#     ]
#
#     data_dir = tmp_path / 'data'
#     data_dir.mkdir()
#
#     target_file = data_dir / 'target.vcf'
#     annotation_file = data_dir / 'annot.vcf'
#
#     _to_vcf(target_data, target_file)
#     _to_vcf(annotation_data, annotation_file)
#
#     target = Vcf(target_file, tmp_path).bgzip()
#     target.index()
#     annotation = Vcf(annotation_file, tmp_path).bgzip()
#     annotation.index()
#
#     result = target.annotate(annotation, 'INFO/AF')
#
#     record = result.to_df('%CHROM\t%POS\t%ID\t%REF\t%ALT\t%INFO\n').iloc[0, :]
#
#     assert 'AX-100' == record['id']
#     assert 'AF=0.5' == record['info']
#
#
# def test_annotate__match_ba(tmp_path):
#     target_data = [
#         "##fileformat=VCFv4.1",
#         "##contig=<ID=chr21>",
#         "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO",
#         "chr21	1000000	AX-100	A	C	.	.	.",
#     ]
#     annotation_data = [
#         "##fileformat=VCFv4.1",
#         "##contig=<ID=chr21>",
#         "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO",
#         "chr21	1000000	rs100	A	C	.	.	.",
#     ]
#
#     data_dir = tmp_path / 'data'
#     data_dir.mkdir()
#
#     target_file = data_dir / 'target.vcf'
#     annotation_file = data_dir / 'annot.vcf'
#
#     _to_vcf(target_data, target_file)
#     _to_vcf(annotation_data, annotation_file)
#
#     target = Vcf(target_file, tmp_path).bgzip()
#     target.index()
#     annotation = Vcf(annotation_file, tmp_path).bgzip()
#     annotation.index()
#
#     result = target.annotate(annotation, 'ID')
#
#     record = result.to_df('%CHROM\t%POS\t%ID\t%REF\t%ALT\t%INFO\n').iloc[0, :]
#     assert 'rs100' == record['id']
#
#
# def test_annotate__match_ma(tmp_path):
#     target_data = [
#         "##fileformat=VCFv4.1",
#         "##contig=<ID=chr21>",
#         "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO",
#         "chr21	1000000	AX-100	A	C,G	.	.	.",
#     ]
#     annotation_data = [
#         "##fileformat=VCFv4.1",
#         "##contig=<ID=chr21>",
#         "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO",
#         "chr21	1000000	rs100	A	C,G	.	.	.",
#     ]
#
#     data_dir = tmp_path / 'data'
#     data_dir.mkdir()
#
#     target_file = data_dir / 'target.vcf'
#     annotation_file = data_dir / 'annot.vcf'
#
#     _to_vcf(target_data, target_file)
#     _to_vcf(annotation_data, annotation_file)
#
#     target = Vcf(target_file, tmp_path).bgzip()
#     target.index()
#     annotation = Vcf(annotation_file, tmp_path).bgzip()
#     annotation.index()
#
#     result = target.annotate(annotation, 'ID')
#
#     record = result.to_df('%CHROM\t%POS\t%ID\t%REF\t%ALT\t%INFO\n').iloc[0, :]
#     assert 'rs100' == record['id']
#
#
# def test_annotate__ma_children(tmp_path):
#     target_data = [
#         "##fileformat=VCFv4.1",
#         "##contig=<ID=chr21>",
#         "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO",
#         "chr21	1000000	AX-100	A	C,G	.	.	.",
#         "chr21	1000000	AX-101	A	C	.	.	.",
#         "chr21	1000000	AX-102	A	G	.	.	.",
#     ]
#     annotation_data = [
#         "##fileformat=VCFv4.1",
#         "##contig=<ID=chr21>",
#         "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO",
#         "chr21	1000000	rs100	A	C,G	.	.	.",
#     ]
#
#     data_dir = tmp_path / 'data'
#     data_dir.mkdir()
#
#     target_file = data_dir / 'target.vcf'
#     annotation_file = data_dir / 'annot.vcf'
#
#     _to_vcf(target_data, target_file)
#     _to_vcf(annotation_data, annotation_file)
#
#     target = Vcf(target_file, tmp_path).bgzip()
#     target.index()
#     annotation = Vcf(annotation_file, tmp_path).bgzip()
#     annotation.index()
#
#     result = target.annotate(annotation, 'ID')
#
#     df = result.to_df('%CHROM\t%POS\t%ID\t%REF\t%ALT\t%INFO\n')
#
#     assert 'rs100' == df.iloc[0, :]['id']
#     assert 'rs100' == df.iloc[1, :]['id']
#     assert 'rs100' == df.iloc[2, :]['id']
#
#
# def test_annotate__subset(tmp_path):
#     target_data = [
#         "##fileformat=VCFv4.1",
#         "##contig=<ID=chr21>",
#         "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO",
#         "chr21	1000000	AX-100	A	C,G	.	.	.",
#     ]
#     annotation_data = [
#         "##fileformat=VCFv4.1",
#         "##contig=<ID=chr21>",
#         "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO",
#         "chr21	1000000	rs100	A	C	.	.	.",
#     ]
#
#     data_dir = tmp_path / 'data'
#     data_dir.mkdir()
#
#     target_file = data_dir / 'target.vcf'
#     annotation_file = data_dir / 'annot.vcf'
#
#     _to_vcf(target_data, target_file)
#     _to_vcf(annotation_data, annotation_file)
#
#     target = Vcf(target_file, tmp_path).bgzip()
#     target.index()
#     annotation = Vcf(annotation_file, tmp_path).bgzip()
#     annotation.index()
#
#     result = target.annotate(annotation, 'ID')
#
#     df = result.to_df('%CHROM\t%POS\t%ID\t%REF\t%ALT\t%INFO\n')
#     assert 'rs100' == df.iloc[0, :]['id']
#
#
# def test_annotate__superset(tmp_path):
#     target_data = [
#         "##fileformat=VCFv4.1",
#         "##contig=<ID=chr21>",
#         "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO",
#         "chr21	1000000	AX-100	A	G	.	.	.",
#     ]
#     annotation_data = [
#         "##fileformat=VCFv4.1",
#         "##contig=<ID=chr21>",
#         "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO",
#         "chr21	1000000	rs100	A	C,G	.	.	.",
#     ]
#
#     data_dir = tmp_path / 'data'
#     data_dir.mkdir()
#
#     target_file = data_dir / 'target.vcf'
#     annotation_file = data_dir / 'annot.vcf'
#
#     _to_vcf(target_data, target_file)
#     _to_vcf(annotation_data, annotation_file)
#
#     target = Vcf(target_file, tmp_path).bgzip()
#     target.index()
#     annotation = Vcf(annotation_file, tmp_path).bgzip()
#     annotation.index()
#
#     result = target.annotate(annotation, 'ID')
#
#     df = result.to_df('%CHROM\t%POS\t%ID\t%REF\t%ALT\t%INFO\n')
#     assert 'rs100' == df.iloc[0, :]['id']
#
#
# def test_annotate__span(tmp_path):
#     target_data = [
#         "##fileformat=VCFv4.1",
#         "##contig=<ID=chr21>",
#         "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO",
#         "chr21	1000000	.	A	C,*	.	.	.",
#         "chr21	1000010	.	A	*,C	.	.	.",
#     ]
#     annotation_data = [
#         "##fileformat=VCFv4.1",
#         "##contig=<ID=chr21>",
#         "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO",
#         "chr21	1000000	AX-100	A	C	.	.	.",
#         "chr21	1000010	AX-200	A	C	.	.	.",
#     ]
#
#     data_dir = tmp_path / 'data'
#     data_dir.mkdir()
#
#     target_file = data_dir / 'target.vcf'
#     annotation_file = data_dir / 'annot.vcf'
#
#     _to_vcf(target_data, target_file)
#     _to_vcf(annotation_data, annotation_file)
#
#     target = Vcf(target_file, tmp_path).bgzip()
#     target.index()
#     annotation = Vcf(annotation_file, tmp_path).bgzip()
#     annotation.index()
#
#     result = target.annotate(annotation, 'ID')
#     df = result.to_df('%CHROM\t%POS\t%ID\t%REF\t%ALT\t%INFO\n')
#
#     assert 'AX-100' == df.iloc[0, :]['id']
#     assert 'AX-200' == df.iloc[1, :]['id']
#
#
# def test_annotate__dup_target(tmp_path):
#     target_data = [
#         "##fileformat=VCFv4.1",
#         "##contig=<ID=chr21>",
#         "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO",
#         "chr21	1000000	AX-100	G	C	.	.	.",
#         "chr21	1000000	AX-100	G	C	.	.	.",
#     ]
#     annotation_data = [
#         "##fileformat=VCFv4.1",
#         "##contig=<ID=chr21>",
#         "##INFO=<ID=AF,Number=A,Type=Float>",
#         "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO",
#         "chr21	1000000	.	G	C	.	.	AF=0.5",
#     ]
#
#     data_dir = tmp_path / 'data'
#     data_dir.mkdir()
#
#     target_file = data_dir / 'target.vcf'
#     annotation_file = data_dir / 'annot.vcf'
#
#     _to_vcf(target_data, target_file)
#     _to_vcf(annotation_data, annotation_file)
#
#     target = Vcf(target_file, tmp_path).bgzip()
#     target.index()
#     annotation = Vcf(annotation_file, tmp_path).bgzip()
#     annotation.index()
#
#     result = target.annotate(annotation, 'INFO/AF')
#     df = result.to_df('%CHROM\t%POS\t%ID\t%REF\t%ALT\t%INFO\n')
#
#     assert 'AX-100' == df.iloc[0, :]['id']
#     assert 'AF=0.5' == df.iloc[0, :]['info']
#     assert 'AX-100' == df.iloc[1, :]['id']
#     assert 'AF=0.5' == df.iloc[1, :]['info']
#
#
# def test_annotate__dup_annot(tmp_path):
#     target_data = [
#         "##fileformat=VCFv4.1",
#         "##contig=<ID=chr21>",
#         "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO",
#         "chr21	1000000	.	T	G	.	.	.",
#     ]
#     annotation_data = [
#         "##fileformat=VCFv4.1",
#         "##contig=<ID=chr21>",
#         "##INFO=<ID=AF,Number=A,Type=Float>",
#         "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO",
#         "chr21	1000000	AX-100	T	G	.	.	.",
#         "chr21	1000000	AX-101	T	G	.	.	.",
#     ]
#
#     data_dir = tmp_path / 'data'
#     data_dir.mkdir()
#
#     target_file = data_dir / 'target.vcf'
#     annotation_file = data_dir / 'annot.vcf'
#
#     _to_vcf(target_data, target_file)
#     _to_vcf(annotation_data, annotation_file)
#
#     target = Vcf(target_file, tmp_path).bgzip()
#     target.index()
#     annotation = Vcf(annotation_file, tmp_path).bgzip()
#     annotation.index()
#
#     result = target.annotate(annotation, 'ID')
#
#     df = result.to_df('%CHROM\t%POS\t%ID\t%REF\t%ALT\t%INFO\n')
#     assert 'AX-100' == df.iloc[0, :]['id']
#
#
# def test_annotate__dup_annot_ma(tmp_path):
#     target_data = [
#         "##fileformat=VCFv4.1",
#         "##contig=<ID=chr21>",
#         "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO",
#         "chr21	1000000	.	A	C	.	.	.",
#         "chr21	1000000	.	A	G	.	.	.",
#     ]
#     annotation_data = [
#         "##fileformat=VCFv4.1",
#         "##contig=<ID=chr21>",
#         "##INFO=<ID=AF,Number=A,Type=Float>",
#         "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO",
#         "chr21	1000000	AX-100	A	C	.	.	.",
#         "chr21	1000000	AX-101	A	C	.	.	.",
#         "chr21	1000000	AX-200	A	G	.	.	.",
#         "chr21	1000000	AX-201	A	G	.	.	.",
#     ]
#
#     data_dir = tmp_path / 'data'
#     data_dir.mkdir()
#
#     target_file = data_dir / 'target.vcf'
#     annotation_file = data_dir / 'annot.vcf'
#
#     _to_vcf(target_data, target_file)
#     _to_vcf(annotation_data, annotation_file)
#
#     target = Vcf(target_file, tmp_path).bgzip().index()
#     annotation = Vcf(annotation_file, tmp_path).bgzip().index()
#
#     result = target.annotate(annotation, 'ID')
#
#     df = result.to_df('%CHROM\t%POS\t%ID\t%REF\t%ALT\t%INFO\n')
#     assert 'AX-100' == df.iloc[0, :]['id']
#     assert 'AX-200' == df.iloc[1, :]['id']

# def test_annotate__dup_target_ma(tmp_path):
#
#     data_dir = tmp_path / 'data'
#     data_dir.mkdir()
#
#     target_data = [
#         "##fileformat=VCFv4.1",
#         "##contig=<ID=chr21>",
#         "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO",
#         "chr21	1000100	AX-100	T	C,G	.	.	.",
#         "chr21	1000100	AX-101	T	C,G	.	.	.",
#         "chr21	1000100	AX-200	T	C	.	.	.",
#         "chr21	1000100	AX-201	T	C	.	.	.",
#         "chr21	1000100	AX-300	T	G	.	.	.",
#         "chr21	1000100	AX-301	T	G	.	.	.",
#     ]
#     annotation_data = [
#         "##fileformat=VCFv4.1",
#         "##contig=<ID=chr21>",
#         "##INFO=<ID=AF,Number=A,Type=Float>",
#         "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO",
#         "chr21	1000100	.	T	C,G	.	.	AF=0.1,0.2",
#         "chr21	1000100	.	T	C	.	.	AF=0.1",
#         "chr21	1000100	.	T	G	.	.	AF=0.2",
#     ]
#
#     target_file = data_dir / 'target.vcf'
#     annotation_file = data_dir / 'annot.vcf'
#
#     _to_vcf(target_data, target_file)
#     _to_vcf(annotation_data, annotation_file)
#
#     target = Vcf(target_file, tmp_path).bgzip()
#     target.index()
#     annotation = Vcf(annotation_file, tmp_path).bgzip()
#     annotation.index()
#
#     result = target.annotate(annotation, 'INFO/AF')
#
#     df = result.to_df('%CHROM\t%POS\t%ID\t%REF\t%ALT\t%INFO\n')
#
#     assert 'AX-100' == df.iloc[0, :]['id']
#     assert 'C,G' == df.iloc[0, :]['alt']
#     assert 'AF=0.1,0.2' == df.iloc[0, :]['info']
#
#     assert 'AX-200' == df.iloc[1, :]['id']
#     assert 'C' == df.iloc[1, :]['alt']
#     assert 'AF=0.1' == df.iloc[1, :]['info']
#
#     assert 'AX-300' == df.iloc[2, :]['id']
#     assert 'G' == df.iloc[2, :]['alt']
#     assert 'AF=0.2' == df.iloc[2, :]['info']
#
#     assert 'AX-101' == df.iloc[3, :]['id']
#     assert 'C,G' == df.iloc[3, :]['alt']
#     assert 'AF=0.1,0.2' == df.iloc[3, :]['info']
#
#     assert 'AX-201' == df.iloc[4, :]['id']
#     assert 'C' == df.iloc[4, :]['alt']
#     assert 'AF=0.1' == df.iloc[4, :]['info']
#
#     assert 'AX-301' == df.iloc[5, :]['id']
#     assert 'G' == df.iloc[5, :]['alt']
#     assert 'AF=0.2' == df.iloc[5, :]['info']


def _to_vcf(data, filepath):
    with filepath.open('wt') as fd:
        for line in data:
            fd.write(line)
            fd.write('\n')


# def test__get_prefix_suffix():
#
#     assert ('01234', '89') == _get_prefix_suffix('0123456789', 100, 109,
#                                                  '567', 105)
#     assert ('', '3456789') == _get_prefix_suffix('0123456789', 100, 109,
#                                                  '012', 100)
#     assert ('0123456', '') == _get_prefix_suffix('0123456789', 100, 109,
#                                                  '789', 107)

# def test_sync_alleles(tmp_path: Path):
#     fixture_dir = Path(__file__).parents[0] / 'fixture_sync_alleles'
#
#     vcf_file_1 = Vcf(fixture_dir / 'one.vcf', tmp_path).bgzip()
#     vcf_file_2 = Vcf(fixture_dir / 'two.vcf', tmp_path).bgzip()
#
#     vcf_file_1_df = vcf_file_1.to_df('%CHROM\t%POS\t%ID\t%REF\t%ALT\n')
#     vcf_file_2_df = vcf_file_2.to_df('%CHROM\t%POS\t%ID\t%REF\t%ALT\n')
#
#     assert {'AGGAGTC'} == set(vcf_file_1_df['ref'])
#     assert {'AGGAGTC', 'A'} == set(vcf_file_2_df['ref'])
#
#     vcf_file_1_result, vcf_file_2_result = sync_alleles(
#         vcf_file_2.filepath, vcf_file_1.filepath, tmp_path / 'after')
#
#     vcf1_observed = Vcf(vcf_file_1_result, tmp_path
#                         / 'after').to_df('%CHROM\t%POS\t%ID\t%REF\t%ALT\n')
#     vcf2_observed = Vcf(vcf_file_2_result, tmp_path
#                         / 'after').to_df('%CHROM\t%POS\t%ID\t%REF\t%ALT\n')
#
#     assert {'AGGAGTC'} == set(vcf1_observed['ref'])
#     assert {'AGGAGTC'} == set(vcf2_observed['ref'])


def test_filter_variants():
    snvs = [
        Variant(chrom='chr1', pos=100, id_='rs100', ref='A', alt='C'),
        Variant(chrom='chr1', pos=200, id_='rs101', ref='T', alt='G'),
    ]

    ref_snv = Variant(chrom='chr1', pos=200, id_='AX-100', ref='T', alt='G')

    expected = Variant(chrom='chr1', pos=200, id_='rs101', ref='T', alt='G')

    actual = filter_variants(ref_snv, snvs)

    assert 1 == len(actual)
    assert expected == actual[0]

    snvs = [
        Variant(chrom='chr1', pos=100, id_='rs100', ref='A', alt='C'),
        Variant(chrom='chr1', pos=200, id_='rs101', ref='C', alt='G'),
    ]

    ref_snv = Variant(chrom='chr1', pos=200, id_='AX-100', ref='C', alt='A')

    actual = filter_variants(ref_snv, snvs)
    assert 0 == len(actual)

    expected = Variant(chrom='chr1', pos=200, id_='rs101', ref='C', alt='A,G')

    actual = filter_variants(ref_snv, snvs, fuzzy=True)
    assert 1 == len(actual)

    snvs = [
        Variant(chrom='chr1', pos=100, id_='rs101', ref='A', alt='C,G'),
    ]

    ref_snv = Variant(chrom='chr1', pos=100, id_='AX-100', ref='A', alt='T')

    expected = Variant(chrom='chr1', pos=100, id_='rs101', ref='A', alt='C,G')

    assert len(filter_variants(ref_snv, snvs, fuzzy=False)) == 0
    assert expected == filter_variants(ref_snv, snvs, fuzzy=True)[0]


def test_split_rtrim():
    snv = Variant(chrom='chr1', pos=100, id_='rs100', ref='AT', alt='CG')

    result = split_rtrim(snv)

    assert len(result)


def test_to_variants(tmp_path):
    vcf_file = Path(__file__).parents[0] / 'fixture/sample_id.vcf.bgz'

    variants = Vcf(vcf_file, tmp_path).to_variants(key='id')
    assert 'AX-100' in variants.keys()
    assert 'AX-101' in variants.keys()
    assert 'AX-200' in variants.keys()
    assert 'AX-201' in variants.keys()
    assert 'AX-300' in variants.keys()
    assert len(variants.keys()) == 5

    variants = Vcf(vcf_file, tmp_path).to_variants(key='coordinate')
    assert len(variants.keys()) == 3


def _write_bgz(path, text, tmp_path):
    plain = tmp_path / path
    plain.write_text(text)
    return Vcf(plain, tmp_path).bgzip()   # -> <name>.vcf.bgz, unindexed


def test_concat_allow_overlaps_false_needs_no_index(tmp_path):
    # Per-sample truth pieces are one-per-chromosome and disjoint, so concat can
    # run without --allow-overlaps and without indexing the inputs.
    a = _write_bgz(
        'a.vcf',
        '##fileformat=VCFv4.2\n'
        '##contig=<ID=chr1>\n'
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">\n'
        '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n'
        'chr1\t100\t.\tA\tG\t.\t.\t.\tGT\t0/1\n',
        tmp_path,
    )
    b = _write_bgz(
        'b.vcf',
        '##fileformat=VCFv4.2\n'
        '##contig=<ID=chr2>\n'
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">\n'
        '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n'
        'chr2\t200\t.\tC\tT\t.\t.\t.\tGT\t1/1\n',
        tmp_path,
    )
    # inputs are intentionally NOT indexed
    assert not (a.filepath.parent / f'{a.filepath.name}.csi').exists()
    assert not (b.filepath.parent / f'{b.filepath.name}.csi').exists()

    out = tmp_path / 'combined.vcf.bgz'
    result = concat([a.filepath, b.filepath], out, tmp_path / 'work',
                    preprocess=False, allow_overlaps=False)

    import subprocess
    rows = subprocess.run(
        ['bcftools', 'query', '-f', '%CHROM\t%POS\n', str(result.filepath)],
        capture_output=True, text=True, check=True).stdout.split()
    assert 'chr1' in rows and 'chr2' in rows
    assert '100' in rows and '200' in rows


def test_concat_writes_no_intermediate(tmp_path):
    a = _write_bgz(
        'a.vcf',
        '##fileformat=VCFv4.2\n'
        '##contig=<ID=chr1>\n'
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">\n'
        '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n'
        'chr1\t100\t.\tA\tG\t.\t.\t.\tGT\t0/1\n',
        tmp_path,
    )
    b = _write_bgz(
        'b.vcf',
        '##fileformat=VCFv4.2\n'
        '##contig=<ID=chr2>\n'
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">\n'
        '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n'
        'chr2\t200\t.\tC\tT\t.\t.\t.\tGT\t1/1\n',
        tmp_path,
    )

    work = tmp_path / 'work'
    out = tmp_path / 'out.vcf.bgz'
    result = concat([a.filepath, b.filepath], out, work,
                    preprocess=False, allow_overlaps=False)

    assert not list(work.rglob('*-concat.vcf.bgz'))
    assert not list(work.rglob('*-sort.vcf.bgz'))

    import subprocess
    rows = subprocess.run(
        ['bcftools', 'query', '-f', '%CHROM\t%POS\n', str(result.filepath)],
        capture_output=True, text=True, check=True).stdout.split()
    assert 'chr1' in rows and 'chr2' in rows
    assert '100' in rows and '200' in rows


def test_sort_plain_vcf_outputs_bgz(tmp_path):
    plain = tmp_path / 'unsorted.vcf'
    plain.write_text(
        '##fileformat=VCFv4.2\n'
        '##contig=<ID=chr1>\n'
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">\n'
        '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n'
        'chr1\t200\t.\tA\tG\t.\t.\t.\tGT\t0/1\n'
        'chr1\t100\t.\tC\tT\t.\t.\t.\tGT\t1/1\n'
    )

    res = Vcf(plain, tmp_path).sort()

    assert str(res.filepath).endswith('-sort.vcf.bgz')
    assert res.filepath.exists()

    import subprocess
    positions = subprocess.run(
        ['bcftools', 'query', '-f', '%POS\n', str(res.filepath)],
        capture_output=True, text=True, check=True).stdout.split()
    assert positions == ['100', '200']


def test_concat_preprocess_allow_overlaps(tmp_path):
    a = tmp_path / 'a.vcf'
    a.write_text(
        '##fileformat=VCFv4.2\n'
        '##contig=<ID=chr1>\n'
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">\n'
        '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n'
        'chr1\t100\t.\tA\tG\t.\t.\t.\tGT\t0/1\n'
        'chr1\t300\t.\tA\tG\t.\t.\t.\tGT\t0/1\n'
    )
    b = tmp_path / 'b.vcf'
    b.write_text(
        '##fileformat=VCFv4.2\n'
        '##contig=<ID=chr1>\n'
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">\n'
        '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n'
        'chr1\t200\t.\tC\tT\t.\t.\t.\tGT\t1/1\n'
        'chr1\t400\t.\tC\tT\t.\t.\t.\tGT\t1/1\n'
    )

    out = tmp_path / 'out.vcf.bgz'
    result = concat([a, b], out, tmp_path / 'work',
                    preprocess=True, allow_overlaps=True)

    assert result.filepath.exists()
    assert (result.filepath.parent / f'{result.filepath.name}.csi').exists()

    import subprocess
    positions = subprocess.run(
        ['bcftools', 'query', '-f', '%POS\n', str(result.filepath)],
        capture_output=True, text=True, check=True).stdout.split()
    assert positions == ['100', '200', '300', '400']


def test_concat_preprocess_tmp_is_readable(tmp_path):
    # concat(preprocess=True) must not create md5-named working dirs.
    import re
    a = tmp_path / 'a.vcf'
    a.write_text(
        '##fileformat=VCFv4.2\n##contig=<ID=chr1>\n'
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">\n'
        '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n'
        'chr1\t100\t.\tA\tG\t.\t.\t.\tGT\t0/1\n'
    )
    b = tmp_path / 'b.vcf'
    b.write_text(
        '##fileformat=VCFv4.2\n##contig=<ID=chr1>\n'
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">\n'
        '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n'
        'chr1\t200\t.\tC\tT\t.\t.\t.\tGT\t1/1\n'
    )
    work = tmp_path / 'work'
    concat([a, b], tmp_path / 'out.vcf.bgz', work, preprocess=True)

    md5 = [p.name for p in work.rglob('*')
           if p.is_dir() and re.fullmatch(r'[0-9a-f]{32}', p.name)]
    assert md5 == [], f'md5-named tmp dirs found: {md5}'


def test_split_by_samples_pieces_are_not_indexed(tmp_path):
    src = _write_bgz(
        'multi.vcf',
        '##fileformat=VCFv4.2\n'
        '##contig=<ID=chr1>\n'
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">\n'
        '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\n'
        'chr1\t100\t.\tA\tG\t.\t.\t.\tGT\t0/1\t1/1\n',
        tmp_path,
    ).index()

    pieces = Vcf(src.filepath, tmp_path).split_by_samples()

    assert set(pieces) == {'S1', 'S2'}
    for path in pieces.values():
        p = Path(path)
        assert p.exists()
        assert not (p.parent / f'{p.name}.csi').exists()
        assert not (p.parent / f'{p.name}.tbi').exists()
