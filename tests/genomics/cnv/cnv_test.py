from genomics.cnv import validate, sort_matches, chop_region
from pathlib import Path

def test_chop_region():
    record = {
            'region_idx': 100,
            'start': 100,
            'end': 200,
    }
    matches = list()
    matches.append({ 'start': 50, 'end': 110, })
    matches.append({ 'start': 120, 'end': 140, })
    matches.append({ 'start': 150, 'end': 210, })

    fragments = chop_region(record, matches)

    assert len(fragments) == 2
    assert fragments[0]['start'] == 110
    assert fragments[0]['end'] == 120
    assert fragments[1]['start'] == 140
    assert fragments[1]['end'] == 150

def test_sort_matches():

    data = list()
    data.append({'chrom': 'chr1', 'start': 100, 'end': 200, 'category': 'TCR'})
    data.append({'chrom': 'chr1', 'start': 150, 'end': 250, 'category': 'TELEMERE'})
    data.append({'chrom': 'chr1', 'start': 300, 'end': 400, 'category': 'BCR'})

    result = sort_matches(data)

    assert len(result) == 2
    assert result[0]['start'] == 100
    assert result[0]['end'] == 250
    assert result[0]['category'] == 'TCR;TELEMERE'
    assert result[1]['start'] == 300
    assert result[1]['end'] == 400
    assert result[1]['category'] == 'BCR'

    data = list()
    data.append({'chrom': 'chr1', 'start': 100, 'end': 200, 'category': 'TCR'})
    data.append({'chrom': 'chr1', 'start': 50, 'end': 250, 'category': 'TELEMERE'})

    result = sort_matches(data)

    assert len(result) == 1
    assert result[0]['start'] == 50
    assert result[0]['end'] == 250
    assert result[0]['category'] == 'TELEMERE;TCR'



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
        breakpoint_tolerance_cutoff=100000,
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


