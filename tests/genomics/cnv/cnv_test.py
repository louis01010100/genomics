from genomics.cnv import (
    validate, sort_matches, exclude_regions, include_regions,
    _load_data, _filter_common_samples, _prepare_dataframe,
    _add_query_fields, _add_match_fields, _add_null_match_fields,
    _compute_and_add_tests, _create_no_match_result, _create_match_result,
    _count_matches_in_range, find_common_samples,
    _handle_no_region_overlap, _handle_region_overlap,
)
from genomics.gregion import GenomicRegion
from pathlib import Path
import polars as pl
import pytest

def test_include_regions():
    record = {
            'region_idx': 100,
            'chrom': 'chr1',
            'start': 100,
            'end': 200,
    }
    matches = list()
    matches.append({ 'chrom': 'chr1', 'start': 50, 'end': 110, })
    matches.append({ 'chrom': 'chr1', 'start': 120, 'end': 140, })
    matches.append({'chrom': 'chr1',  'start': 150, 'end': 210, })

    fragments = include_regions(record, matches)

    assert len(fragments) == 3
    assert fragments[0]['start'] == 100
    assert fragments[0]['end'] == 110
    assert fragments[1]['start'] == 120
    assert fragments[1]['end'] == 140
    assert fragments[2]['start'] == 150
    assert fragments[2]['end'] == 200

def test_exclude_regions():
    record = {
            'region_idx': 100,
            'chrom': 'chr1',
            'start': 100,
            'end': 200,
    }
    matches = list()
    matches.append({ 'chrom': 'chr1', 'start': 50, 'end': 110, })
    matches.append({ 'chrom': 'chr1', 'start': 120, 'end': 140, })
    matches.append({'chrom': 'chr1',  'start': 150, 'end': 210, })

    fragments = exclude_regions(record, matches)

    assert len(fragments) == 2
    assert fragments[0]['start'] == 110
    assert fragments[0]['end'] == 120
    assert fragments[1]['start'] == 140
    assert fragments[1]['end'] == 150

    record = {
            'region_idx': 100,
            'start': 100,
            'end': 200,
    }
    matches = list()
    matches.append({ 'start': 110, 'end': 150})

    fragments = exclude_regions(record, matches)

    assert len(fragments) == 2
    assert fragments[0]['start'] == 100
    assert fragments[0]['end'] == 110
    assert fragments[1]['start'] == 150
    assert fragments[1]['end'] == 200

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
        n_threads=32,
    )

    # with (tmp_path / 'regions_ppv-0.5.tsv').open('rt') as fh:
    #     i2c = {c: i for i, c in enumerate(next(fh).strip().split('\t'))}
    #
    #     values = next(fh).strip().split('\t')
    #
    #     assert values[i2c['chrom']] == 'chr1'
    #     assert values[i2c['start']] == '9010000'
    #     assert values[i2c['end']] == '11000000'
    #     assert values[i2c['length']] == '1990000'


# ============================================================================
# Tests for validate() helper functions
# ============================================================================

def test_load_data(tmp_path):
    """Test _load_data() function."""
    # Create test files
    predictions_file = tmp_path / "predictions.tsv"
    truths_file = tmp_path / "truths.tsv"

    predictions_file.write_text(
        "sample_name\tchrom\tstart\tend\tcn_state\n"
        "sample1\tchr1\t100\t200\t3\n"
    )
    truths_file.write_text(
        "sample_name\tchrom\tstart\tend\tcn_state\n"
        "sample1\tchr1\t100\t200\t3\n"
    )

    predictions, truths = _load_data(predictions_file, truths_file)

    assert len(predictions) == 1
    assert len(truths) == 1
    assert predictions['sample_name'][0] == 'sample1'
    assert truths['sample_name'][0] == 'sample1'


def test_prepare_dataframe():
    """Test _prepare_dataframe() function."""
    df = pl.DataFrame({
        'sample_name': ['sample1', 'sample2'],
        'chrom': ['chr1', 'chr1'],
        'start': ['100', '200'],
        'end': ['150', '250'],
    })

    result = _prepare_dataframe(df)

    assert result['start'].dtype == pl.Int64
    assert result['end'].dtype == pl.Int64
    assert 'region_idx' in result.columns
    assert list(result['region_idx']) == [0, 1]


def test_filter_common_samples():
    """Test _filter_common_samples() function."""
    predictions = pl.DataFrame({
        'sample_name': ['sample1', 'sample2', 'sample3'],
        'chrom': ['chr1', 'chr1', 'chr1'],
        'start': [100, 200, 300],
        'end': [150, 250, 350],
        'cn_state': ['3', '3', '3'],
    })

    truths = pl.DataFrame({
        'sample_name': ['sample1', 'sample4'],
        'chrom': ['chr1', 'chr1'],
        'start': [100, 400],
        'end': [150, 450],
        'cn_state': ['3', '3'],
    })

    sample_map = {'sample1': 'sample1', 'sample2': 'sample4'}

    pred_filtered, truth_filtered = _filter_common_samples(predictions, truths, sample_map)

    assert len(pred_filtered) == 2
    assert len(truth_filtered) == 2
    assert set(pred_filtered['sample_name']) == {'sample1', 'sample2'}
    assert set(truth_filtered['sample_name']) == {'sample1', 'sample4'}


def test_find_common_samples():
    """Test find_common_samples() function."""
    prediction_samples = ['sample1', 'sample2', 'sample3']
    truth_samples = ['sampleA', 'sampleB']
    sample_map = {'sample1': 'sampleA', 'sample2': 'sampleB', 'sample3': 'sampleC'}

    common_pred, common_truth = find_common_samples(
        prediction_samples, truth_samples, sample_map
    )

    assert common_pred == ['sample1', 'sample2']
    assert common_truth == ['sampleA', 'sampleB']


# ============================================================================
# Tests for _validate_cnv() helper functions
# ============================================================================

def test_add_query_fields():
    """Test _add_query_fields() function."""
    result = {}
    query = {
        'fragment_idx': 'frag_1',
        'chrom': 'chr1',
        'start': 100,
        'end': 200,
        'cn_state': '3',
        'n_markers': 50,
    }

    _add_query_fields(result, 'prediction', query)

    assert result['prediction_fragment_idx'] == 'frag_1'
    assert result['chrom'] == 'chr1'
    assert result['prediction_start'] == 100
    assert result['prediction_end'] == 200
    assert result['prediction_length'] == 100
    assert result['prediction_cn_state'] == '3'
    assert result['prediction_n_markers'] == 50


def test_add_match_fields():
    """Test _add_match_fields() function."""
    result = {}
    match = {
        'fragment_idx': 'frag_2',
        'start': 110,
        'end': 210,
        'cn_state': '3',
        'n_markers': 60,
    }

    _add_match_fields(result, 'truth', match)

    assert result['truth_fragment_idx'] == 'frag_2'
    assert result['truth_start'] == 110
    assert result['truth_end'] == 210
    assert result['truth_length'] == 100
    assert result['truth_cn_state'] == '3'
    assert result['truth_n_markers'] == 60


def test_add_null_match_fields():
    """Test _add_null_match_fields() function."""
    result = {}

    _add_null_match_fields(result, 'truth')

    assert result['truth_fragment_idx'] is None
    assert result['truth_start'] is None
    assert result['truth_end'] is None
    assert result['truth_length'] is None
    assert result['truth_cn_state'] is None
    assert result['truth_n_markers'] is None


def test_compute_and_add_tests():
    """Test _compute_and_add_tests() function."""
    result = {}
    query_region = GenomicRegion('chr1', 100, 200)
    db_region = GenomicRegion('chr1', 150, 250)
    query = {'cn_state': '3'}
    match = {'cn_state': '3'}

    _compute_and_add_tests(
        result, query_region, db_region, query, match,
        reciprocal_overlap_cutoff=0.3,
        breakpoint_tolerance_cutoff=100
    )

    assert 'reciprocal_overlap' in result
    assert 'breakpoint_difference' in result
    assert result['reciprocal_overlap_test'] in ['PASS', 'FAIL']
    assert result['breakpoint_tolerance_test'] in ['PASS', 'FAIL']
    assert result['cn_state_test'] == 'PASS'
    assert result['matched'] in [True, False]


def test_create_no_match_result():
    """Test _create_no_match_result() function."""
    query = {
        'fragment_idx': 'frag_1',
        'chrom': 'chr1',
        'start': 100,
        'end': 200,
        'cn_state': '3',
    }

    result = _create_no_match_result('prediction', 'truth', query)

    assert result['prediction_fragment_idx'] == 'frag_1'
    assert result['truth_fragment_idx'] is None
    assert result['matched'] is None
    assert result['reciprocal_overlap'] is None


def test_create_match_result():
    """Test _create_match_result() function."""
    query = {
        'fragment_idx': 'frag_1',
        'chrom': 'chr1',
        'start': 100,
        'end': 200,
        'cn_state': '3',
    }
    match = {
        'fragment_idx': 'frag_2',
        'chrom': 'chr1',
        'start': 150,
        'end': 250,
        'cn_state': '3',
    }
    query_region = GenomicRegion('chr1', 100, 200)

    result = _create_match_result(
        'prediction', 'truth', query, match, query_region,
        reciprocal_overlap_cutoff=0.3,
        breakpoint_tolerance_cutoff=100
    )

    assert result['prediction_fragment_idx'] == 'frag_1'
    assert result['truth_fragment_idx'] == 'frag_2'
    assert 'reciprocal_overlap' in result
    assert 'matched' in result


# ============================================================================
# Tests for report_by_size() helper functions
# ============================================================================

def test_count_matches_in_range():
    """Test _count_matches_in_range() function."""
    data = pl.DataFrame({
        'prediction_length': [50000, 75000, 150000, 500000],
        'prediction_n_markers': [100, 200, 300, 400],
        'matched': [True, False, True, True],
    })

    record, value = _count_matches_in_range(
        data, 'prediction', 50000, 100000, 0, 1000
    )

    assert record['min_length'] == 50000
    assert record['max_length'] == 100000
    assert record['n_detected'] == 2
    assert record['n_matched'] == 1
    assert value == 0.5


def test_count_matches_in_range_no_matches():
    """Test _count_matches_in_range() with no matches."""
    data = pl.DataFrame({
        'prediction_length': [500000, 600000],
        'prediction_n_markers': [100, 200],
        'matched': [True, False],
    })

    record, value = _count_matches_in_range(
        data, 'prediction', 50000, 100000, 0, 1000
    )

    assert record['n_detected'] == 0
    assert record['n_matched'] == 0
    assert value is None


# ============================================================================
# Tests for region annotation helper functions
# ============================================================================

def test_handle_no_region_overlap():
    """Test _handle_no_region_overlap() function."""
    record = {
        'region_idx': 100,
        'chrom': 'chr1',
        'start': 100,
        'end': 200,
        'cn_state': '3',
    }

    region, fragment = _handle_no_region_overlap(record, 'hotspot_region')

    assert region['hotspot_region'] is None
    assert fragment['fragment_idx'] == 100
    assert fragment['chrom'] == 'chr1'


def test_handle_region_overlap():
    """Test _handle_region_overlap() function."""
    record = {
        'region_idx': 100,
        'chrom': 'chr1',
        'start': 100,
        'end': 300,
    }
    matches = [
        {'chrom': 'chr1', 'start': 120, 'end': 180, 'region_name': 'hotspot1'},
        {'chrom': 'chr1', 'start': 200, 'end': 250, 'region_name': 'hotspot2'},
    ]

    region, fragments = _handle_region_overlap(
        record, matches, 'region_name', 'hotspot_region', use_include=True
    )

    assert 'hotspot1' in region['hotspot_region']
    assert 'hotspot2' in region['hotspot_region']
    assert len(fragments) >= 1


