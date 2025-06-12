from pathlib import Path

import polars as pl

GFF3_COLUMNS = [
    'seqid',
    'source',
    'type',
    'start',
    'end',
    'score',
    'strand',
    'phase',
    'attribute',
]

UCSC_GENE_COLUMNS = [
    'bin',
    'name',
    'chrom',
    'strand',
    'tx_start',
    'tx_end',
    'cds_start',
    'cds_end',
    'exon_count',
    'exon_starts',
    'exon_ends',
    'score',
    'name2',
    'cds_start_stat',
    'cds_end_stat',
    'exon_frames',
]

def export(input_file: Path, output_file, data_source: str, one_based):
    if data_source == 'UCSC':
        export_ucsc_gene(input_file, output_file, one_based)
    elif data_source == 'GENCODE':
        export_gencode_gene(input_file, output_file, one_based)

    else:
        raise Exception(data_source)

        

def export_gencode_gene(input_file: Path, output_file: Path, one_based):
    data = pl.read_csv(input_file, comment_prefix = '#', separator = '\t', new_columns = GFF3_COLUMNS)
    data = data.with_columns(pl.arange(0, len(data)).alias('idx'))

    genes = data.filter(pl.col('type') == 'gene')

    genes = genes.with_columns(pl.col('attribute').str.split(';').alias('attribute')).explode(pl.col('attribute'))
    genes = genes.with_columns(pl.col('attribute').str.split('=').alias('split'))\
            .with_columns(
                    pl.col('split').list.get(0).alias('key'), \
                    pl.col('split').list.get(1).alias('gene_name'), \
            )
    genes = genes.filter(pl.col('key') == 'gene_name')
    genes = genes.rename({'seqid': 'chrom', 'start': 'default_start', 'end': 'default_end'})
    if one_based:
        genes = genes.with_columns((pl.col('default_start')).alias('start'), pl.col('default_end').alias('end'))
    else:
        genes = genes.with_columns((pl.col('default_start') - 1).alias('start'), pl.col('default_end').alias('end'))
    genes = genes.select(['gene_name', 'chrom', 'start', 'end'])
    genes.write_csv(output_file, include_header = True, separator = '\t')


def export_ucsc_gene(
    input_file: Path,
    output_file:Path,
    one_based: bool,
    columns=['name2', 'name', 'chrom', 'tx_start', 'tx_end'],
):
    genes = pl.read_csv(
        input_file,
        has_header=False,
        new_columns=UCSC_GENE_COLUMNS,
        separator='\t',
    )

    genes = genes.select(pl.col(columns))

    # return genes.unique()

    bag = []
    for keys, data in genes.groupby('name2', 'chrom'):

        gene = keys[0]
        chrom = keys[1]

        if len(data) == 1:
            start = data['tx_start'][0]
            end = data['tx_end'][0]

        else:
            start = data['tx_start'].min()
            end = data['tx_end'].max()

        if one_based:
            start += 1

        bag.append({
            'gene': gene,
            'chrom': chrom,
            'start': start,
            'end': end,
        })

    data = pl.from_dicts(bag)
    data.write_csv(output_file, include_header = True, separator = '\t')

