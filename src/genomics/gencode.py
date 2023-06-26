from pathlib import Path

import polars as pl

GENE_COLUMNS = [
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


def load(
    gene_file: Path,
    zero_based: bool = True,
    columns=['name2', 'name', 'chrom', 'tx_start', 'tx_end'],
):
    genes = pl.read_csv(
        gene_file,
        has_header=False,
        new_columns=GENE_COLUMNS,
        separator='\t',
    )

    genes = genes.select(pl.col(columns))

    # return genes.unique()

    bag = []
    for gene, data in genes.groupby('name2'):

        if len(data) == 1:
            start = data['tx_start'][0]
            end = data['tx_end'][0]

        else:
            start = data['tx_start'].min()
            end = data['tx_end'].max()

        if not zero_based:
            start += 1

        bag.append({
            'gene': gene,
            'start': start,
            'end': end,
        })

    data = pl.from_dicts(bag)

    return data
