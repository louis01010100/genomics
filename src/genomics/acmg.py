from pathlib import Path

import polars as pl


def process(input_file: Path, output_file, columns=['gene']):
    genes = pl.read_csv(input_file, has_header=True, separator='\t')
    genes.columns = [x.lower().replace(' ', '_') for x in genes.columns]
    genes = genes.select(pl.col(columns)).unique()
    genes.write_csv(output_file, has_header=True, separator='\t')
