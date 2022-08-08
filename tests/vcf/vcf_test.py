#!/usr/bin/env python
import pandas as pd
from vcf import Vcf
from pathlib import Path


def test_meta(tmp_path):
    vcf_file = Path(__file__).parents[0] / 'fixture/sample.vcf'
    vcf = Vcf(vcf_file, tmp_path)
    expected = (
        '##fileformat=VCFv4.1\n'
        '##FILTER=<ID=PASS,Description="All filters passed">\n'
        '##contig=<ID=chr21>\n'
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Phased Genotype">')
    assert expected == vcf.meta


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


def test_to_df(tmp_path):
    vcf_file = Path(__file__).parents[0] / 'fixture/sample.vcf.bgz'

    vcf = Vcf(vcf_file, tmp_path)
    df = vcf.to_df()
    assert 991 == len(df)


def test_to_tsv(tmp_path):
    vcf_file = Path(__file__).parents[0] / 'fixture/sample.vcf.bgz'
    vcf = Vcf(vcf_file, tmp_path)
    output_file = vcf.to_tsv(
        '[%CHROM\t%POS\t%ID\t%REF\t%ALT\t%SAMPLE\t%GT\n]')

    data = pd.read_csv(output_file, header=0, sep='\t', dtype='str')

    record1 = data.iloc[0, :]

    assert record1['CHROM'] == 'chr21'
    assert record1['POS'] == '5030578'
    assert record1['ID'] == '.'
    assert record1['REF'] == 'C'
    assert record1['ALT'] == 'T'
    assert record1['SAMPLE'] == 'HG00403'
    assert record1['GT'] == '0|0'


def test_annotate_tag(tmp_path):
    target_data = [
        "##fileformat=VCFv4.1",
        "##contig=<ID=chr21>",
        "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO",
        "chr21	1000000	AX-100	G	C	.	.	.",
    ]
    annotation_data = [
        "##fileformat=VCFv4.1",
        "##contig=<ID=chr21>",
        "##INFO=<ID=AF,Number=A,Type=Float>",
        "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO",
        "chr21	1000000	.	G	C	.	.	AF=0.5",
    ]

    data_dir = tmp_path / 'data'
    data_dir.mkdir()

    target_file = data_dir / 'target.vcf'
    annotation_file = data_dir / 'annot.vcf'

    _to_vcf(target_data, target_file)
    _to_vcf(annotation_data, annotation_file)

    target = Vcf(target_file, tmp_path).bgzip()
    target.index()
    annotation = Vcf(annotation_file, tmp_path).bgzip()
    annotation.index()

    result = target.annotate(annotation, 'INFO/AF')

    assert 'AX-100' == result.to_df().iloc[0, :]['ID']
    assert 'AF=0.5' == result.to_df().iloc[0, :]['INFO']


def test_annotate__match_ba(tmp_path):
    target_data = [
        "##fileformat=VCFv4.1",
        "##contig=<ID=chr21>",
        "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO",
        "chr21	1000000	.	A	C	.	.	.",
    ]
    annotation_data = [
        "##fileformat=VCFv4.1",
        "##contig=<ID=chr21>",
        "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO",
        "chr21	1000000	AX-100	A	C	.	.	.",
    ]

    data_dir = tmp_path / 'data'
    data_dir.mkdir()

    target_file = data_dir / 'target.vcf'
    annotation_file = data_dir / 'annot.vcf'

    _to_vcf(target_data, target_file)
    _to_vcf(annotation_data, annotation_file)

    target = Vcf(target_file, tmp_path).bgzip()
    target.index()
    annotation = Vcf(annotation_file, tmp_path).bgzip()
    annotation.index()

    result = target.annotate(annotation, 'ID')

    assert 'AX-100' == result.to_df().iloc[0, :]['ID']


def test_annotate__match_ma(tmp_path):
    target_data = [
        "##fileformat=VCFv4.1",
        "##contig=<ID=chr21>",
        "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO",
        "chr21	1000000	.	A	C,G	.	.	.",
    ]
    annotation_data = [
        "##fileformat=VCFv4.1",
        "##contig=<ID=chr21>",
        "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO",
        "chr21	1000000	AX-100	A	C,G	.	.	.",
    ]

    data_dir = tmp_path / 'data'
    data_dir.mkdir()

    target_file = data_dir / 'target.vcf'
    annotation_file = data_dir / 'annot.vcf'

    _to_vcf(target_data, target_file)
    _to_vcf(annotation_data, annotation_file)

    target = Vcf(target_file, tmp_path).bgzip()
    target.index()
    annotation = Vcf(annotation_file, tmp_path).bgzip()
    annotation.index()

    result = target.annotate(annotation, 'ID')

    assert 'AX-100' == result.to_df().iloc[0, :]['ID']


def test_annotate__match_ma_children(tmp_path):
    target_data = [
        "##fileformat=VCFv4.1",
        "##contig=<ID=chr21>",
        "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO",
        "chr21	1000000	.	A	C	.	.	.",
        "chr21	1000000	.	A	G	.	.	.",
    ]
    annotation_data = [
        "##fileformat=VCFv4.1",
        "##contig=<ID=chr21>",
        "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO",
        "chr21	1000000	AX-100	A	C	.	.	.",
        "chr21	1000000	AX-101	A	G	.	.	.",
    ]

    data_dir = tmp_path / 'data'
    data_dir.mkdir()

    target_file = data_dir / 'target.vcf'
    annotation_file = data_dir / 'annot.vcf'

    _to_vcf(target_data, target_file)
    _to_vcf(annotation_data, annotation_file)

    target = Vcf(target_file, tmp_path).bgzip()
    target.index()
    annotation = Vcf(annotation_file, tmp_path).bgzip()
    annotation.index()

    result = target.annotate(annotation, 'ID')

    assert 'AX-100' == result.to_df().iloc[0, :]['ID']
    assert 'AX-101' == result.to_df().iloc[1, :]['ID']


def test_annotate__subset(tmp_path):
    target_data = [
        "##fileformat=VCFv4.1",
        "##contig=<ID=chr21>",
        "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO",
        "chr21	1000000	.	A	C,G	.	.	.",
    ]
    annotation_data = [
        "##fileformat=VCFv4.1",
        "##contig=<ID=chr21>",
        "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO",
        "chr21	1000000	AX-100	A	C	.	.	.",
    ]

    data_dir = tmp_path / 'data'
    data_dir.mkdir()

    target_file = data_dir / 'target.vcf'
    annotation_file = data_dir / 'annot.vcf'

    _to_vcf(target_data, target_file)
    _to_vcf(annotation_data, annotation_file)

    target = Vcf(target_file, tmp_path).bgzip()
    target.index()
    annotation = Vcf(annotation_file, tmp_path).bgzip()
    annotation.index()

    result = target.annotate(annotation, 'ID')

    assert 'AX-100' == result.to_df().iloc[0, :]['ID']


def test_annotate__superset(tmp_path):
    target_data = [
        "##fileformat=VCFv4.1",
        "##contig=<ID=chr21>",
        "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO",
        "chr21	1000000	.	A	G	.	.	.",
    ]
    annotation_data = [
        "##fileformat=VCFv4.1",
        "##contig=<ID=chr21>",
        "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO",
        "chr21	1000000	AX-100	A	C,G	.	.	.",
    ]

    data_dir = tmp_path / 'data'
    data_dir.mkdir()

    target_file = data_dir / 'target.vcf'
    annotation_file = data_dir / 'annot.vcf'

    _to_vcf(target_data, target_file)
    _to_vcf(annotation_data, annotation_file)

    target = Vcf(target_file, tmp_path).bgzip()
    target.index()
    annotation = Vcf(annotation_file, tmp_path).bgzip()
    annotation.index()

    result = target.annotate(annotation, 'ID')

    assert 'AX-100' == result.to_df().iloc[0, :]['ID']


def test_annotate__span(tmp_path):
    target_data = [
        "##fileformat=VCFv4.1",
        "##contig=<ID=chr21>",
        "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO",
        "chr21	1000000	.	A	C,*	.	.	.",
        "chr21	1000010	.	A	*,C	.	.	.",
    ]
    annotation_data = [
        "##fileformat=VCFv4.1",
        "##contig=<ID=chr21>",
        "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO",
        "chr21	1000000	AX-100	A	C	.	.	.",
        "chr21	1000010	AX-200	A	C	.	.	.",
    ]

    data_dir = tmp_path / 'data'
    data_dir.mkdir()

    target_file = data_dir / 'target.vcf'
    annotation_file = data_dir / 'annot.vcf'

    _to_vcf(target_data, target_file)
    _to_vcf(annotation_data, annotation_file)

    target = Vcf(target_file, tmp_path).bgzip()
    target.index()
    annotation = Vcf(annotation_file, tmp_path).bgzip()
    annotation.index()

    result = target.annotate(annotation, 'ID')

    assert 'AX-100' == result.to_df().iloc[0, :]['ID']
    assert 'AX-200' == result.to_df().iloc[1, :]['ID']


def test_annotate__dup_target(tmp_path):
    target_data = [
        "##fileformat=VCFv4.1",
        "##contig=<ID=chr21>",
        "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO",
        "chr21	1000000	AX-100	G	C	.	.	.",
        "chr21	1000000	AX-100	G	C	.	.	.",
    ]
    annotation_data = [
        "##fileformat=VCFv4.1",
        "##contig=<ID=chr21>",
        "##INFO=<ID=AF,Number=A,Type=Float>",
        "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO",
        "chr21	1000000	.	G	C	.	.	AF=0.5",
    ]

    data_dir = tmp_path / 'data'
    data_dir.mkdir()

    target_file = data_dir / 'target.vcf'
    annotation_file = data_dir / 'annot.vcf'

    _to_vcf(target_data, target_file)
    _to_vcf(annotation_data, annotation_file)

    target = Vcf(target_file, tmp_path).bgzip()
    target.index()
    annotation = Vcf(annotation_file, tmp_path).bgzip()
    annotation.index()

    result = target.annotate(annotation, 'INFO/AF')

    assert 'AX-100' == result.to_df().iloc[0, :]['ID']
    assert 'AF=0.5' == result.to_df().iloc[0, :]['INFO']
    assert 'AX-100' == result.to_df().iloc[1, :]['ID']
    assert 'AF=0.5' == result.to_df().iloc[1, :]['INFO']

    target_data = [
        "##fileformat=VCFv4.1",
        "##contig=<ID=chr21>",
        "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO",
        "chr21	1000000	AX-100	G	C	.	.	.",
        "chr21	1000000	AX-200	G	C	.	.	.",
    ]
    annotation_data = [
        "##fileformat=VCFv4.1",
        "##contig=<ID=chr21>",
        "##INFO=<ID=AF,Number=A,Type=Float>",
        "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO",
        "chr21	1000000	.	G	C	.	.	AF=0.5",
    ]

    data_dir = tmp_path / 'data'
    data_dir.mkdir(exist_ok=True)

    target_file = data_dir / 'target.vcf'
    annotation_file = data_dir / 'annot.vcf'

    _to_vcf(target_data, target_file)
    _to_vcf(annotation_data, annotation_file)

    target = Vcf(target_file, tmp_path).bgzip()
    target.index()
    annotation = Vcf(annotation_file, tmp_path).bgzip()
    annotation.index()

    result = target.annotate(annotation, 'INFO')

    assert 'AX-100' == result.to_df().iloc[0, :]['ID']
    assert 'AF=0.5' == result.to_df().iloc[0, :]['INFO']
    assert 'AX-200' == result.to_df().iloc[1, :]['ID']
    assert 'AF=0.5' == result.to_df().iloc[1, :]['INFO']


def test_annotate__dup_annot(tmp_path):
    target_data = [
        "##fileformat=VCFv4.1",
        "##contig=<ID=chr21>",
        "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO",
        "chr21	1000000	.	T	G	.	.	.",
    ]
    annotation_data = [
        "##fileformat=VCFv4.1",
        "##contig=<ID=chr21>",
        "##INFO=<ID=AF,Number=A,Type=Float>",
        "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO",
        "chr21	1000000	AX-100	T	G	.	.	.",
        "chr21	1000000	AX-101	T	G	.	.	.",
    ]

    data_dir = tmp_path / 'data'
    data_dir.mkdir()

    target_file = data_dir / 'target.vcf'
    annotation_file = data_dir / 'annot.vcf'

    _to_vcf(target_data, target_file)
    _to_vcf(annotation_data, annotation_file)

    target = Vcf(target_file, tmp_path).bgzip()
    target.index()
    annotation = Vcf(annotation_file, tmp_path).bgzip()
    annotation.index()

    result = target.annotate(annotation, 'ID')

    assert 'AX-100' == result.to_df().iloc[0, :]['ID']


def test_annotate__dup_annot_ma(tmp_path):
    target_data = [
        "##fileformat=VCFv4.1",
        "##contig=<ID=chr21>",
        "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO",
        "chr21	1000000	.	A	C	.	.	.",
        "chr21	1000000	.	A	G	.	.	.",
    ]
    annotation_data = [
        "##fileformat=VCFv4.1",
        "##contig=<ID=chr21>",
        "##INFO=<ID=AF,Number=A,Type=Float>",
        "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO",
        "chr21	1000000	AX-100	A	C	.	.	.",
        "chr21	1000000	AX-101	A	C	.	.	.",
        "chr21	1000000	AX-200	A	G	.	.	.",
        "chr21	1000000	AX-201	A	G	.	.	.",
    ]

    data_dir = tmp_path / 'data'
    data_dir.mkdir()

    target_file = data_dir / 'target.vcf'
    annotation_file = data_dir / 'annot.vcf'

    _to_vcf(target_data, target_file)
    _to_vcf(annotation_data, annotation_file)

    target = Vcf(target_file, tmp_path).bgzip().index()
    annotation = Vcf(annotation_file, tmp_path).bgzip().index()

    result = target.annotate(annotation, 'ID')

    assert 'AX-100' == result.to_df().iloc[0, :]['ID']
    assert 'AX-200' == result.to_df().iloc[1, :]['ID']


def test_annotate__dup_target_ma(tmp_path):
    target_data = [
        "##fileformat=VCFv4.1",
        "##contig=<ID=chr21>",
        "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO",
        "chr21	1000100	AX-100	T	C,G	.	.	.",
        "chr21	1000100	AX-101	T	C	.	.	.",
    ]
    annotation_data = [
        "##fileformat=VCFv4.1",
        "##contig=<ID=chr21>",
        "##INFO=<ID=AF,Number=A,Type=Float>",
        "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO",
        "chr21	1000100	.	T	C,G	.	.	AF=0.1,0.2",
    ]

    data_dir = tmp_path / 'data'
    data_dir.mkdir()

    target_file = data_dir / 'target.vcf'
    annotation_file = data_dir / 'annot.vcf'

    _to_vcf(target_data, target_file)
    _to_vcf(annotation_data, annotation_file)

    target = Vcf(target_file, tmp_path).bgzip()
    target.index()
    annotation = Vcf(annotation_file, tmp_path).bgzip()
    annotation.index()

    result = target.annotate(annotation, 'ID', 'INFO/AF')

    result = result.to_df()

    assert 'AX-100' == result.iloc[0, :]['ID']
    assert 'C,G' == result.iloc[0, :]['ALT']
    assert 'AF=0.1,0.2' == result.iloc[0, :]['INFO']

    assert 'AX-101' == result.iloc[1, :]['ID']
    assert 'C' == result.iloc[1, :]['ALT']
    assert 'AF=0.1' == result.iloc[1, :]['INFO']

    target_data = [
        "##fileformat=VCFv4.1",
        "##contig=<ID=chr21>",
        "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO",
        "chr21	1000100	AX-100	T	C,G	.	.	.",
        "chr21	1000100	AX-101	T	C,G	.	.	.",
        "chr21	1000100	AX-200	T	C	.	.	.",
        "chr21	1000100	AX-201	T	C	.	.	.",
        "chr21	1000100	AX-300	T	G	.	.	.",
        "chr21	1000100	AX-301	T	G	.	.	.",
    ]
    annotation_data = [
        "##fileformat=VCFv4.1",
        "##contig=<ID=chr21>",
        "##INFO=<ID=AF,Number=A,Type=Float>",
        "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO",
        "chr21	1000100	.	T	C,G	.	.	AF=0.1,0.2",
        "chr21	1000100	.	T	C	.	.	AF=0.1",
        "chr21	1000100	.	T	G	.	.	AF=0.2",
    ]

    target_file = data_dir / 'target.vcf'
    annotation_file = data_dir / 'annot.vcf'

    _to_vcf(target_data, target_file)
    _to_vcf(annotation_data, annotation_file)

    target = Vcf(target_file, tmp_path).bgzip()
    target.index()
    annotation = Vcf(annotation_file, tmp_path).bgzip()
    annotation.index()

    result = target.annotate(annotation, 'INFO/AF')

    print(tmp_path)

    result = result.to_df()

    assert 'AX-100' == result.iloc[0, :]['ID']
    assert 'C,G' == result.iloc[0, :]['ALT']
    assert 'AF=0.1,0.2' == result.iloc[0, :]['INFO']

    assert 'AX-200' == result.iloc[1, :]['ID']
    assert 'C' == result.iloc[1, :]['ALT']
    assert 'AF=0.1' == result.iloc[1, :]['INFO']

    assert 'AX-300' == result.iloc[2, :]['ID']
    assert 'G' == result.iloc[2, :]['ALT']
    assert 'AF=0.2' == result.iloc[2, :]['INFO']

    assert 'AX-101' == result.iloc[3, :]['ID']
    assert 'C,G' == result.iloc[3, :]['ALT']
    assert 'AF=0.1,0.2' == result.iloc[3, :]['INFO']

    assert 'AX-201' == result.iloc[4, :]['ID']
    assert 'C' == result.iloc[4, :]['ALT']
    assert 'AF=0.1' == result.iloc[4, :]['INFO']

    assert 'AX-301' == result.iloc[5, :]['ID']
    assert 'G' == result.iloc[5, :]['ALT']
    assert 'AF=0.2' == result.iloc[5, :]['INFO']


def _to_vcf(data, filepath):
    with filepath.open('wt') as fd:
        for line in data:
            fd.write(line)
            fd.write('\n')
