from pathlib import Path

import polars as pl


def load(acmg_file: Path, columns=['gene']):
    genes = pl.read_csv(acmg_file, has_header=True, separator='\t')

    genes.columns = [x.lower().replace(' ', '_') for x in genes.columns]

    genes = genes.select(pl.col(columns)).unique()

    return genes
