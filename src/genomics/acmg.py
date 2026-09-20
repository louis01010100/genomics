from collections import OrderedDict
from datetime import datetime
from importlib.metadata import version as _get_version
from pathlib import Path

import polars as pl

from .utils import init_logging, log_start, log_stop


def process(input_file: Path, output_file, columns=['gene']):
    start = datetime.now()
    log_dir = Path(output_file).parent
    log_dir.mkdir(parents=True, exist_ok=True)
    init_logging(log_dir / 'acmg.log')
    banner = f'genomics acmg {_get_version("genomics")}'
    log_start(banner=banner, info=OrderedDict([
        ('input-file', input_file),
        ('output-file', output_file),
    ]))
    genes = pl.read_csv(input_file, has_header=True, separator='\t')
    genes.columns = [x.lower().replace(' ', '_') for x in genes.columns]
    genes = genes.select(pl.col(columns)).unique()
    genes.write_csv(output_file, has_header=True, separator='\t')
    log_stop(banner, start, datetime.now())
