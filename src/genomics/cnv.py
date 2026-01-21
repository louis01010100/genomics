from importlib import resources
import numpy as np
from copy import deepcopy

import polars as pl
from pathos.multiprocessing import ProcessPool

from genomics.gregion import GenomicRegion, create_database

from . import karyopype as kp
from .utils import chrom_sort_key

COMPLEX_REGIONS_FILE = (
    resources.files("genomics") / "resources" / "complex_regions-hg38.tsv"
)
CHROM_LENGTH_HG38_FILE = (
    resources.files("genomics") / "resources" / "chrom_lengths-hg38.tsv"
)


REQUIRED_COLUMNS = {
    "sample_name",
    "chrom",
    "start",
    "end",
    "cn_state",
}

HOTSPOT_COLUMNS = {
    'region_name',
    'chrom',
    'start',
    'end',
}

def validate(
    predictions_file,
    truths_file,
    sample_map_file,
    output_dir,
    hotspot_regions_file=None,
    reciprocal_overlap_cutoff=0.6,
    breakpoint_tolerance_cutoff=100000,
    n_threads=32,
):
    """Validate CNV predictions against truth data."""
    # Load and prepare data
    predictions, truths = _load_data(predictions_file, truths_file)
    sample_map = load_sample_map(sample_map_file)
    predictions, truths = _filter_common_samples(predictions, truths, sample_map)
    predictions = _prepare_dataframe(predictions)
    truths = _prepare_dataframe(truths)

    # Annotate regions and write to files
    pred_regions, pred_fragments, truth_regions, truth_fragments = \
        _annotate_and_get_regions(predictions, truths, hotspot_regions_file)
    _write_region_files(pred_regions, pred_fragments, truth_regions, truth_fragments, output_dir)

    # Validate samples
    samples = group_by_sample(pred_fragments, truth_fragments, sample_map)
    prediction_vs_truth, truth_vs_prediction = _validate_all_samples(
        samples, reciprocal_overlap_cutoff, breakpoint_tolerance_cutoff
    )

    # Write, reload, and cast results
    _write_validation_results(prediction_vs_truth, truth_vs_prediction, output_dir)

    # prediction_vs_truth, truth_vs_prediction = _reload_validation_results(output_dir)
    # prediction_vs_truth, truth_vs_prediction = _cast_validation_columns(
    #     prediction_vs_truth, truth_vs_prediction
    # )

    # Generate final reports
    _generate_reports(prediction_vs_truth, truth_vs_prediction, output_dir)


# ============================================================================
# Helper functions for validate()
# ============================================================================

def _load_data(predictions_file, truths_file):
    """Load prediction and truth files and validate required columns."""
    predictions = pl.read_csv(
        predictions_file, has_header=True, separator="\t", infer_schema=False
    )
    truths = pl.read_csv(
        truths_file, has_header=True, separator="\t", infer_schema=False
    )
    assert REQUIRED_COLUMNS.issubset(predictions.columns), predictions.columns
    assert REQUIRED_COLUMNS.issubset(truths.columns), truths.columns
    return predictions, truths


def _filter_common_samples(predictions, truths, sample_map):
    """Filter predictions and truths to only include common samples."""
    common_samples_prediction, common_samples_truth = find_common_samples(
        predictions['sample_name'], truths['sample_name'], sample_map
    )
    predictions = predictions.filter(
        pl.col('sample_name').is_in(common_samples_prediction)
    ).sort(REQUIRED_COLUMNS)
    truths = truths.filter(
        pl.col('sample_name').is_in(common_samples_truth)
    ).sort(REQUIRED_COLUMNS)
    return predictions, truths


def _prepare_dataframe(df):
    """Cast types, sort, and add region_idx to dataframe."""
    df = df.with_columns(
        pl.col("start").cast(pl.Int64),
        pl.col("end").cast(pl.Int64),
    )
    df = df.sort(['sample_name', 'chrom', 'start'])
    df = df.with_columns(
        pl.arange(0, len(df)).alias("region_idx"),
    )
    return df


def _annotate_and_get_regions(predictions, truths, hotspot_regions_file):
    """Annotate regions with hotspot/complex regions."""
    if hotspot_regions_file:
        hotspot_region_db = create_database(load_hotspot_regions(hotspot_regions_file))
        prediction_regions, prediction_fragments = annotate_hotspot_regions(
            predictions, hotspot_region_db
        )
        truth_regions, truth_fragments = annotate_hotspot_regions(truths, hotspot_region_db)
    else:
        complex_region_db = create_database(load_complex_regions(COMPLEX_REGIONS_FILE))
        prediction_regions, prediction_fragments = annotate_complex_regions(
            predictions, complex_region_db
        )
        truth_regions, truth_fragments = annotate_complex_regions(truths, complex_region_db)

    return prediction_regions, prediction_fragments, truth_regions, truth_fragments


def _write_region_files(prediction_regions, prediction_fragments,
                        truth_regions, truth_fragments, output_dir):
    """Write region and fragment files to disk."""
    prediction_regions.write_csv(
        output_dir / "prediction_regions.tsv", include_header=True, separator="\t"
    )
    prediction_fragments.write_csv(
        output_dir / "prediction_fragments.tsv", include_header=True, separator="\t"
    )
    truth_regions.write_csv(
        output_dir / "truth_regions.tsv", include_header=True, separator="\t"
    )
    truth_fragments.write_csv(
        output_dir / "truth_fragments.tsv", include_header=True, separator="\t"
    )


def _validate_all_samples(samples, reciprocal_overlap_cutoff, breakpoint_tolerance_cutoff):
    """Validate CNVs for all samples and return comparison results."""
    prediction_vs_truth = []
    truth_vs_prediction = []

    for sample in samples:
        sample_name_prediction = sample["sample_name_prediction"]
        sample_name_truth = sample["sample_name_truth"]

        prediction_fragments = sample["prediction_fragments"].filter(
            ~pl.col('chrom').is_null()
        ).to_dicts()
        truth_fragments = sample["truth_fragments"].filter(
            ~pl.col('chrom').is_null()
        ).to_dicts()

        data_ppv = _validate_cnv(
            prediction_fragments, create_database(truth_fragments),
            "prediction", "truth",
            reciprocal_overlap_cutoff, breakpoint_tolerance_cutoff
        )
        if data_ppv is not None:
            data_ppv = data_ppv.with_columns(
                pl.lit(sample_name_prediction).alias('sample_name')
            ).select(['sample_name'] + data_ppv.columns[:])
            prediction_vs_truth.append(data_ppv)

        data_sensitivity = _validate_cnv(
            truth_fragments, create_database(prediction_fragments),
            "truth", "prediction",
            reciprocal_overlap_cutoff, breakpoint_tolerance_cutoff
        )
        if data_sensitivity is not None:
            data_sensitivity = data_sensitivity.with_columns(
                pl.lit(sample_name_truth).alias('sample_name')
            ).select(['sample_name'] + data_sensitivity.columns[:])
            truth_vs_prediction.append(data_sensitivity)

    return (pl.concat(prediction_vs_truth, how="vertical_relaxed"),
            pl.concat(truth_vs_prediction, how="vertical_relaxed"))


def _write_validation_results(prediction_vs_truth, truth_vs_prediction, output_dir):
    """Write validation results to TSV files."""
    prediction_vs_truth.sort(['sample_name', 'chrom', 'prediction_start']).write_csv(
        output_dir / "prediction_vs_truth.tsv", include_header=True, separator="\t"
    )
    truth_vs_prediction.sort(['sample_name', 'chrom', 'truth_start']).write_csv(
        output_dir / "truth_vs_prediction.tsv", include_header=True, separator="\t"
    )


def _reload_validation_results(output_dir):
    """Reload validation results from TSV files."""
    prediction_vs_truth = pl.read_csv(
        output_dir / "prediction_vs_truth.tsv",
        has_header=True, separator="\t", infer_schema=False
    )
    truth_vs_prediction = pl.read_csv(
        output_dir / "truth_vs_prediction.tsv",
        has_header=True, separator="\t", infer_schema=False
    )
    return prediction_vs_truth, truth_vs_prediction


def _cast_validation_columns(prediction_vs_truth, truth_vs_prediction):
    """Cast validation result columns to appropriate types."""
    prediction_vs_truth = prediction_vs_truth.with_columns(
        pl.col("prediction_start").cast(pl.Int64).alias("start"),
        pl.col("prediction_end").cast(pl.Int64).alias("end"),
        pl.col("prediction_length").cast(pl.Int64),
        pl.col("reciprocal_overlap").cast(pl.Float64),
    )
    truth_vs_prediction = truth_vs_prediction.with_columns(
        pl.col("truth_start").cast(pl.Int64).alias("start"),
        pl.col("truth_end").cast(pl.Int64).alias("end"),
        pl.col("truth_length").cast(pl.Int64),
        pl.col("reciprocal_overlap").cast(pl.Float64),
    )
    return prediction_vs_truth, truth_vs_prediction


def _generate_reports(prediction_vs_truth, truth_vs_prediction, output_dir):
    """Generate all validation reports."""
    report(prediction_vs_truth, output_dir, 'ppv')
    report(truth_vs_prediction, output_dir, 'sensitivity')
    report_by_size(prediction_vs_truth, output_dir / 'ppv-markers_vs_len.tsv', 'prediction')
    report_by_size(truth_vs_prediction, output_dir / 'sensitivity-markers_vs_len.tsv', 'truth')


# ============================================================================
# Main validation function
# ============================================================================

# 0-based, half-open interval




def create_high_concordance_region(moving_average_data, concordance_cutoff):


    moving_average_data = moving_average_data.filter(
        pl.col("concordance") >= concordance_cutoff
    )
    moving_average_data = moving_average_data.sort(["chrom", "start"])

    regions = list()

    region = None

    for record in moving_average_data.to_dicts():
        fragment = GenomicRegion(
            chrom=record["chrom"],
            start=record["start"],
            end=record["end"],
        )

        if region is None:
            region = fragment
            continue

        if not region.is_close_to(fragment, 1):

            regions.append(
                {
                    "chrom": region.chrom,
                    "start": region.start,
                    "end": region.end,
                    "length": len(region),
                }
            )

            region = GenomicRegion(
                chrom=record["chrom"],
                start=record["start"],
                end=record["end"],
            )
        else:
            region = region.merge(fragment)

    if region:
        regions.append(
            {
                "chrom": region.chrom,
                "start": region.start,
                "end": region.end,
                "length": len(region),
            }
        )

    return pl.from_dicts(regions)


def _create_sliding_window_jobs(data, chrom_orders, chrom_lengths, complex_regions,
                                window_size, step_size, reciprocal_overlap_cutoff):
    """Generate sliding window jobs for each chromosome."""
    for chrom in chrom_orders:
        df = data.filter(pl.col("chrom") == chrom)
        if len(df) == 0:
            continue

        yield {
            "chrom": chrom,
            "chrom_length": chrom_lengths[chrom],
            "data": df.to_dicts(),
            "window_size": window_size,
            "step_size": step_size,
            "reciprocal_overlap_cutoff": reciprocal_overlap_cutoff,
            "complex_regions": complex_regions,
        }


def _process_sliding_window_job(job):
    """Process a single chromosome for sliding window analysis."""
    chrom = job["chrom"]
    chrom_length = job["chrom_length"]
    data = job["data"]
    window_size = job["window_size"]
    step_size = job["step_size"]
    complex_regions = job["complex_regions"]

    database = create_database(data)
    complex_region_db = create_database(complex_regions)

    bag = []

    for start in range(0, chrom_length, step_size):
        end = min(start + window_size, chrom_length)
        query = GenomicRegion(chrom, start, end)

        # Skip complex regions
        if complex_region_db.find_overlap(
            {"chrom": query.chrom, "start": query.start, "end": query.end}
        ):
            continue

        # Find overlapping CNVs
        matches = database.find_overlap(
            {"chrom": query.chrom, "start": query.start, "end": query.end}
        )
        if not matches:
            continue

        # Calculate concordance
        n_detected = 0
        n_matched = 0
        for match in matches:
            n_detected += 1
            if match["matched"] is True:
                n_matched += 1

        concordance = None if n_detected == 0 else n_matched / n_detected

        bag.append({
            "chrom": query.chrom,
            "start": query.start,
            "end": query.end,
            "concordance": concordance,
            "n_detected": n_detected,
            "n_matched": n_matched,
        })
    return bag


def summarize_sliding(
    data,
    chrom_length_file,
    complex_regions_file,
    window_size,
    step_size,
    reciprocal_overlap_cutoff,
    n_threads,
):
    """Summarize CNV concordance using sliding windows across chromosomes."""
    chrom_lengths = load_chrom_lengths(chrom_length_file)
    chrom_orders = sorted(chrom_lengths.keys(), key=chrom_sort_key)
    complex_regions = load_complex_regions(complex_regions_file)

    bag = []

    with ProcessPool(n_threads) as pool:
        for result in pool.uimap(
            _process_sliding_window_job,
            _create_sliding_window_jobs(
                data, chrom_orders, chrom_lengths, complex_regions,
                window_size, step_size, reciprocal_overlap_cutoff
            ),
        ):
            bag.extend(result)

    return pl.from_dicts(bag, infer_schema_length=None)


def load_sample_map(input_file):
    bag = dict()
    with input_file.open("rt") as fh:
        for line in fh:
            items = line.strip().split("\t")
            if items[0] == "sample_name_prediction":
                continue
            bag[items[0]] = items[1]

    return bag


def load_chrom_lengths(input_file):

    bag = dict()
    with input_file.open("rt") as fh:
        for line in fh:
            items = line.strip().split("\t")
            chrom = items[0]
            length = items[1]

            if chrom == "chrom":
                continue

            bag[chrom] = int(length)

    return bag


# ============================================================================
# Helper functions for _validate_cnv()
# ============================================================================

def _add_query_fields(result, qname, query):
    """Add query fields to result dictionary."""
    result[f"{qname}_fragment_idx"] = query["fragment_idx"]
    result["chrom"] = query["chrom"]
    result[f"{qname}_start"] = query["start"]
    result[f"{qname}_end"] = query["end"]
    result[f"{qname}_length"] = query["end"] - query["start"]
    result[f"{qname}_cn_state"] = query["cn_state"]
    result[f"{qname}_n_markers"] = query.get("n_markers")


def _add_match_fields(result, dbname, match):
    """Add match fields to result dictionary."""
    result[f"{dbname}_fragment_idx"] = match["fragment_idx"]
    result[f"{dbname}_start"] = match["start"]
    result[f"{dbname}_end"] = match["end"]
    result[f"{dbname}_length"] = match["end"] - match["start"]
    result[f"{dbname}_cn_state"] = match["cn_state"]
    result[f"{dbname}_n_markers"] = match.get("n_markers")


def _add_null_match_fields(result, dbname):
    """Add null match fields when no match found."""
    result[f"{dbname}_fragment_idx"] = None
    result[f"{dbname}_start"] = None
    result[f"{dbname}_end"] = None
    result[f"{dbname}_length"] = None
    result[f"{dbname}_cn_state"] = None
    result[f"{dbname}_n_markers"] = None


def _compute_and_add_tests(result, query_region, db_region, query, match,
                           reciprocal_overlap_cutoff, breakpoint_tolerance_cutoff):
    """Compute overlap metrics and add test results."""
    reciprocal_overlap = query_region.get_reciprocal_overlap(db_region)['ratio']
    breakpoint_difference = query_region.get_max_boundary_difference(db_region)

    result["reciprocal_overlap"] = reciprocal_overlap
    result["breakpoint_difference"] = breakpoint_difference
    result["reciprocal_overlap_test"] = (
        "PASS" if reciprocal_overlap >= reciprocal_overlap_cutoff else "FAIL"
    )
    result["breakpoint_tolerance_test"] = (
        "PASS" if breakpoint_difference <= breakpoint_tolerance_cutoff else "FAIL"
    )
    result["cn_state_test"] = (
        "PASS" if query["cn_state"] == match["cn_state"] else "FAIL"
    )

    # Determine if matched
    if result["reciprocal_overlap_test"] == "FAIL":
        result["matched"] = False
    elif result["cn_state_test"] == "FAIL":
        result["matched"] = False
    else:
        result["matched"] = True


def _create_no_match_result(qname, dbname, query):
    """Create result dictionary when no match is found."""
    result = {}
    _add_query_fields(result, qname, query)
    _add_null_match_fields(result, dbname)
    result["reciprocal_overlap"] = None
    result["breakpoint_difference"] = None
    result["reciprocal_overlap_test"] = None
    result["breakpoint_tolerance_test"] = None
    result["cn_state_test"] = None
    result["matched"] = None
    return result


def _create_match_result(qname, dbname, query, match, query_region,
                         reciprocal_overlap_cutoff, breakpoint_tolerance_cutoff):
    """Create result dictionary when match is found."""
    result = {}
    db_region = GenomicRegion(match["chrom"], match["start"], match["end"])

    _add_query_fields(result, qname, query)
    _add_match_fields(result, dbname, match)
    _compute_and_add_tests(result, query_region, db_region, query, match,
                          reciprocal_overlap_cutoff, breakpoint_tolerance_cutoff)

    return result


def _validate_cnv(
    queries,
    database,
    qname,
    dbname,
    reciprocal_overlap_cutoff,
    breakpoint_tolerance_cutoff,
):
    """Validate CNVs by comparing queries against a database."""
    bag = []
    for query in queries:
        matches = database.find_overlap(query)
        query_region = GenomicRegion(query["chrom"], query["start"], query["end"])

        if len(matches) == 0:
            result = _create_no_match_result(qname, dbname, query)
            bag.append(result)
        else:
            for match in matches:
                result = _create_match_result(
                    qname, dbname, query, match, query_region,
                    reciprocal_overlap_cutoff, breakpoint_tolerance_cutoff
                )
                bag.append(result)

    return pl.from_dicts(bag, infer_schema_length=None) if bag else None


def group_by_sample(prediction_fragments, truth_fragments, sample_map):

    prediction_bag = dict()
    for keys, df in prediction_fragments.group_by("sample_name"):
        sample_name = keys[0]
        prediction_bag[sample_name] = df

    truth_bag = dict()
    for keys, df in truth_fragments.group_by("sample_name"):
        sample_name = keys[0]
        truth_bag[sample_name] = df


    bag = list()

    for sample_name_prediction, df in prediction_bag.items():

        assert sample_name_prediction in sample_map

        sample_name_truth = sample_map[sample_name_prediction]

        if sample_name_truth not in truth_bag:
            print(f'{sample_name_truth} has no CNV after filtering')
            continue

        bag.append(
            {
                "sample_name_prediction": sample_name_prediction,
                "sample_name_truth": sample_name_truth,
                "prediction_fragments": df,
                "truth_fragments": truth_bag[sample_name_truth],
            }
        )

    return bag

def load_hotspot_regions(hotspot_regions_file):
    hotspot_regions = pl.read_csv(hotspot_regions_file, has_header = True, separator = '\t')
    assert HOTSPOT_COLUMNS.issubset(hotspot_regions.columns)

    return hotspot_regions.to_dicts()



def load_complex_regions(complex_regions_file):

    data = pl.read_csv(complex_regions_file, has_header=True, separator="\t")

    return data.to_dicts()

def _report(data, _type, label, output_dir):
    bag = list()

    if label == 'overall':
        bag = list()

        record = {
            'n_total' : len(data),
            'n_matched' : len(data.filter(pl.col('matched') == True)),
            _type: len(data.filter(pl.col('matched') == True)) / len(data),
        }

        if 'hotspot_overlap_ratio' in data.columns:
            on_hotspot = data.filter(pl.col('hotspot_overlap_ratio').cast(pl.Float64) > 0.5)

            record['n_total_hotspot'] = len(on_hotspot)
            record['n_matched_hotspot'] = len(on_hotspot.filter(pl.col('matched') == True))
            record[f'{_type}_hotspot'] = np.nan if len(on_hotspot) == 0 \
                    else record['n_matched_hotspot'] / record['n_total_hotspot']

        bag.append(record)

        data = pl.from_dicts(bag)

        print(data)

        data.write_csv(output_dir / f'{_type}.tsv', include_header = True, separator = '\t')

        return


    for keys, df in data.group_by([label]):

        record = {
            label: keys[0],
            'n_total': len(df),
            'n_matched': len(df.filter(pl.col('matched') == True)),
            _type: len(df.filter(pl.col('matched') == True)) / len(df),
        }

        if 'hotspot_overlap_ratio' in df.columns:
            on_hotspot = df.filter(pl.col('hotspot_overlap_ratio').cast(pl.Float64) > 0)

            record['n_total_hotspot'] = len(on_hotspot)
            record['n_matched_hotspot'] = len(on_hotspot.filter(pl.col('matched') == True))
            record[f'{_type}_hotspot'] = np.nan if len(on_hotspot) == 0 \
                    else record['n_matched_hotspot'] / record['n_total_hotspot']



        if label == 'barcode':
            record['n_samples'] = len(df['sample_name'].unique())

        bag.append(record)

    data = pl.from_dicts(bag).sort(_type)
    data.write_csv(
        output_dir / f'{_type}_by_{label}.tsv', include_header = True, separator = '\t')

def report(a_vs_b, output_dir, _type):


    if _type == 'ppv':
        frag_idx = 'prediction_fragment_idx'
    elif _type == 'sensitivity':
        frag_idx = 'truth_fragment_idx'


    data = a_vs_b.rename({frag_idx: 'fragment_idx'})

    bag = list()

    for keys, df in data.group_by(['fragment_idx']):
        if len(df) == 1:
            bag.append(df.row(0, named = True))
            continue

        if len(set(df['matched'])) == 1:
            bag.append(df.row(0, named = True))
            continue
        else:
            record = df.row(0, named = True)
            record['matched'] = True
            bag.append(record)

    data = pl.from_dicts(bag, infer_schema_length = None)

    if 'barcode' in data.columns:
        _report(data, _type, 'barcode', output_dir)
    _report(data, _type, 'sample_name', output_dir)
    _report(data, _type, 'overall', output_dir)






# def report_by_cn_segments(data, label):
#
#     lengths =  [0, 50000,100000,500000,1000000, np.inf]
#     n_markers = [0, 50, 100, 200, 300, 400, 500, 600, 700, 800, 900, 
#                     1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000, np.inf]
#
#     bag = list()
#
#     for i in range(0, len(lengths) - 1):
#         min_length = lengths[i]
#         max_length = lengths[i + 1]
#         target = data.filter(
#             (pl.col(f'prediction_length') >= min_length) & (pl.col(f'prediction_length') < max_length)
#         )
#
#         for j in range(0, len(n_markers) - 1):
#
#             min_n_markers = n_markers[j]
#             max_n_markers = n_markers[j + 1]
#
#             target2 = target.filter(
#                 (pl.col(f'prediction_n_markers') >= min_n_markers) & 
#                 (pl.col(f'prediction_n_markers') < max_n_markers)
#             )
#
#             n_detected = len(target2)
#             n_matched = len(target2.filter(pl.col('matched') == True))
#
#             if n_detected == 0:
#                 value = None
#             else:
#                 value = n_matched / n_detected
#
#             bag.append({
#                 'min_length': min_length,
#                 'max_length': max_length,
#                 'min_n_markers': min_n_markers,
#                 'max_n_markers': max_n_markers,
#                 'n_detected': n_detected,
#                 'n_matched': n_matched,
#                 'ppv': value,
#             })
#
#     return pl.from_dicts(bag)

def _count_matches_in_range(data, label, min_length, max_length, min_markers, max_markers):
    """Count matches in a specific size and marker range."""
    target = data.filter(
        (pl.col(f'{label}_length') >= min_length) &
        (pl.col(f'{label}_length') < max_length) &
        (pl.col(f'{label}_n_markers').cast(pl.Int64) >= min_markers) &
        (pl.col(f'{label}_n_markers').cast(pl.Int64) < max_markers)
    )

    n_detected = len(target)
    n_matched = len(target.filter(pl.col('matched') == True))
    value = None if n_detected == 0 else n_matched / n_detected

    return {
        'min_length': min_length,
        'max_length': max_length,
        'min_n_markers': min_markers,
        'max_n_markers': max_markers,
        'n_detected': n_detected,
        'n_matched': n_matched,
    }, value


def report_by_size(data, output_file, label):
    """Generate CNV validation report stratified by size and marker count."""
    metric = 'ppv' if label == 'prediction' else 'sensitivity'
    lengths = [0, 50000, 100000, 250000, 500000, 1000000, np.inf]
    n_markers = [0, np.inf]

    bag = []
    for i in range(len(lengths) - 1):
        for j in range(len(n_markers) - 1):
            record, value = _count_matches_in_range(
                data, label, lengths[i], lengths[i + 1],
                n_markers[j], n_markers[j + 1]
            )
            record[metric] = value
            bag.append(record)

    pl.from_dicts(bag).write_csv(output_file, include_header=True, separator='\t')

def sort_matches(data):
    data.sort(key = lambda x: x['start'])

    bag = list()
    for record in data:
        if not bag:
            bag.append(record)
            continue
        if bag[-1]['end'] > record['start']:

            bag[-1]['end'] = max(bag[-1]['end'], record['end'])
            bag[-1]['category'] = bag[-1]['category'] + ';' + record['category']

        else:
            bag.append(record)

    return bag

# ============================================================================
# Helper functions for region annotation
# ============================================================================

def _handle_no_region_overlap(record, region_label):
    """Handle case when CNV has no overlap with annotated regions."""
    region = record.copy()
    fragment = record.copy()
    region[region_label] = None
    fragment["fragment_idx"] = fragment["region_idx"]
    fragment['fragment_length'] = fragment['end'] - fragment['start'] if fragment['end'] else None
    return region, fragment


def _handle_region_overlap(record, matches, region_name_key, region_label, use_include):
    """Handle case when CNV overlaps with annotated regions."""
    matches = sort_matches(matches)
    region = deepcopy(record)
    region[region_label] = ';'.join([x[region_name_key] for x in matches])

    if use_include:
        fragments = include_regions(record, matches)
    else:
        fragments = exclude_regions(record, matches)

    return region, fragments


def _generic_region_annotation(cnvs, region_db, region_name_key, region_label, use_include):
    """
    Generic region annotation function.

    Args:
        cnvs: DataFrame of CNV calls
        region_db: Database of regions to annotate
        region_name_key: Key to extract region name ('region_name' for hotspot, 'category' for complex)
        region_label: Label for output column ('hotspot_region' or 'complex_region')
        use_include: If True, use include_regions; if False, use exclude_regions

    Returns:
        Tuple of (regions DataFrame, fragments DataFrame)
    """
    cnvs = cnvs.with_columns(
        pl.col("start").cast(pl.Int64),
        pl.col("end").cast(pl.Int64),
    )

    bag_region = []
    bag_fragment = []

    for record in cnvs.to_dicts():
        matches = region_db.find_overlap(record)

        if len(matches) == 0:
            region, fragment = _handle_no_region_overlap(record, region_label)
            bag_region.append(region)

            if not use_include:
                bag_fragment.append(fragment)
        else:
            region, fragments = _handle_region_overlap(
                record, matches, region_name_key, region_label, use_include
            )
            bag_region.append(region)
            bag_fragment.extend(fragments)

    regions = pl.from_dicts(bag_region, infer_schema_length=None)
    fragments = pl.from_dicts(bag_fragment, infer_schema_length=None)

    return regions, fragments


def annotate_hotspot_regions(cnvs, hotspot_region_db):
    """Annotate CNVs with hotspot region information."""
    return _generic_region_annotation(
        cnvs, hotspot_region_db, 'region_name', 'hotspot_region', use_include=True
    )


def annotate_complex_regions(cnvs, complex_region_db):
    """Annotate CNVs with complex region information."""
    return _generic_region_annotation(
        cnvs, complex_region_db, 'category', 'complex_region', use_include=False
    )

def include_regions(record, matches):

    bag = list()

    region_idx = record["region_idx"]
    frag_idx = 0

    region = GenomicRegion(record['chrom'], record['start'], record['end'])

    for match in matches:
        intersect = region.intersects(GenomicRegion(match['chrom'], match['start'], match['end']))
        fragment = deepcopy(record)
        fragment['start'] = intersect.start
        fragment['end'] = intersect.end
        fragment['fragment_length'] = fragment['end'] - fragment['start']
        fragment["fragment_idx"] = f"{region_idx}_{frag_idx}"
        bag.append(fragment)
        frag_idx += 1

    return bag


def exclude_regions(record, matches):

    bag = list()

    region_idx = record["region_idx"]
    start = record['start']
    end = record['end']
    frag_idx = 0


    for match in matches:
        if start < match['start']:
            fragment = deepcopy(record)
            fragment['start'] = start
            fragment['end'] = match['start']
            fragment['fragment_length'] = fragment['end'] - fragment['start']
            fragment["fragment_idx"] = f"{region_idx}_{frag_idx}"
            bag.append(fragment)
            frag_idx += 1
            start = match['end']
        elif match['start'] <= start <= match['end']:
            start = match['end']
        else:
            assert False, f'{start}-{end}-{record}-{match}'


    if start < end:
        fragment = deepcopy(record)
        fragment['start'] = start
        fragment['end'] = end
        fragment['fragment_length'] = fragment['end'] - fragment['start']
        fragment["fragment_idx"] = f"{region_idx}_{frag_idx}"
        bag.append(fragment)

    return bag

def find_common_samples(prediction_samples, truth_samples, sample_map):
    common_samples_p = list()
    common_samples_t = list()


    for p, t in sample_map.items():
        if p not in prediction_samples:
            continue

        if t not in truth_samples:
            continue

        common_samples_p.append(p)
        common_samples_t.append(t)

    return common_samples_p, common_samples_t

