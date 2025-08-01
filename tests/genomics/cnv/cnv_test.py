from genomics.cnv import validate
from pathlib import Path

def test_validate(tmp_path):

    fixture_dir = Path(__file__).parent / 'fixture'
    predictions_file = fixture_dir / 'predictions.tsv'
    truths_file = fixture_dir / 'truths.tsv'
    sample_map_file = fixture_dir / 'sample_map.tsv'

    validate(
        predictions_file,
        truths_file,
        sample_map_file,
        tmp_path,
        reciprocal_overlap_cutoff=0.5,
        boundary_difference_cutoff=10000,
        window_size=1000000,
        step_size=10000,
        concordance_cutoff=0.5,
        n_threads=32,
    )

    with (tmp_path / 'regions_ppv-0.5.tsv').open('rt') as fh:
        i2c = {c: i for i, c in enumerate(next(fh).strip().split('\t'))}

        values = next(fh).strip().split('\t')

        assert values[i2c['chrom']] == 'chr1'
        assert values[i2c['start']] == '9010000'
        assert values[i2c['end']] == '11000000'
        assert values[i2c['length']] == '1990000'


