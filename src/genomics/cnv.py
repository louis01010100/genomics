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

# 0-based, half-open interval
def validate(
    predictions_file,
    truths_file,
    sample_map_file,
    output_dir,
    hotspot_regions_file=None,
    # n_markers_cutoff = 50,
    # cnv_lenght_cutoff = 50000,
    reciprocal_overlap_cutoff=0.6,
    breakpoint_tolerance_cutoff=100000,
    # window_size=1000000,
    # step_size=1000,
    # concordance_cutoff=0.8,
    n_threads=32,
):

    predictions = pl.read_csv(
        predictions_file, has_header=True, separator="\t", infer_schema=False
    )
    truths = pl.read_csv(
        truths_file, has_header=True, separator="\t", infer_schema=False
    )

    assert REQUIRED_COLUMNS.issubset(predictions.columns), predictions.columns
    assert REQUIRED_COLUMNS.issubset(truths.columns), truths.columns

    sample_map = load_sample_map(sample_map_file)

    # assert '20191223_GT23_ARRAY_726_H02_TPMI.CEL' in sample_map
    # assert '00051461_HTCMA.CEL' == sample_map['20191223_GT23_ARRAY_726_H02_TPMI.CEL']
    # assert '20191223_GT23_ARRAY_726_H02_TPMI.CEL' in predictions['sample_name']
    # assert '00051461_HTCMA.CEL' in truths['sample_name']

    common_samples_prediction, common_samples_truth = find_common_samples(
            predictions['sample_name'], truths['sample_name'], sample_map)

    # assert '20191223_GT23_ARRAY_726_H02_TPMI.CEL' in common_samples_prediction
    # assert '00051461_HTCMA.CEL' in common_samples_truth
    # return


    predictions = predictions.filter(pl.col('sample_name').is_in(common_samples_prediction)).sort(REQUIRED_COLUMNS)
    truths = truths.filter(pl.col('sample_name').is_in(common_samples_truth)).sort(REQUIRED_COLUMNS)


    predictions = predictions.with_columns(
        pl.col("start").cast(pl.Int64),
        pl.col("end").cast(pl.Int64),
        # pl.col("n_markers").cast(pl.Int64),
    )
    predictions = predictions.sort(['sample_name', 'chrom', 'start'])
    predictions = predictions.with_columns(
        pl.arange(0, len(predictions)).alias("region_idx"),
    )

    truths = truths.with_columns(
        pl.col("start").cast(pl.Int64),
        pl.col("end").cast(pl.Int64),
        # pl.col("n_markers").cast(pl.Int64),
    )

    truths = truths.sort(['sample_name', 'chrom', 'start'])
    truths = truths.with_columns(
        pl.arange(0, len(truths)).alias("region_idx"),
    )


    if hotspot_regions_file:
        hotspot_region_db = create_database(load_hotspot_regions(hotspot_regions_file))
        predictions = annotate_hotspot_regions(predictions, hotspot_region_db)
        truths = annotate_hotspot_regions(truths, hotspot_region_db)

    complex_region_db = create_database(load_complex_regions(COMPLEX_REGIONS_FILE))

    prediction_regions, prediction_fragments = annotate_complex_regions(
        predictions,
        complex_region_db,
    )


    prediction_regions.write_csv(
        output_dir / "prediction_regions.tsv", include_header=True, separator="\t"
    )

    prediction_fragments.write_csv(
        output_dir / "prediction_fragments.tsv", include_header=True, separator="\t"
    )

    truth_regions, truth_fragments = annotate_complex_regions(truths, complex_region_db)

    truth_regions.write_csv(
        output_dir / "truth_regions.tsv", include_header=True, separator="\t"
    )
    truth_fragments.write_csv(
        output_dir / "truth_fragments.tsv", include_header=True, separator="\t"
    )

    samples = group_by_sample(prediction_fragments, truth_fragments, sample_map)


    prediction_vs_truth = list()
    truth_vs_prediction = list()

    for sample in samples:
        sample_name_prediction = sample["sample_name_prediction"]
        sample_name_truth = sample["sample_name_truth"]
        prediction_fragments = sample["prediction_fragments"].to_dicts()
        truth_fragments = sample["truth_fragments"].to_dicts()

        data_ppv = _validate_cnv(
            prediction_fragments,
            create_database(truth_fragments),
            "prediction",
            "truth",
            reciprocal_overlap_cutoff,
            breakpoint_tolerance_cutoff,
        ).with_columns(
                pl.lit(sample_name_prediction).alias('sample_name'),
        )

        data_ppv = data_ppv.select(['sample_name'] + data_ppv.columns[:-1])

        prediction_vs_truth.append(data_ppv)

        data_sensitivity = _validate_cnv(
            truth_fragments,
            create_database(prediction_fragments),
            "truth",
            "prediction",
            reciprocal_overlap_cutoff,
            breakpoint_tolerance_cutoff,
        ).with_columns(
                pl.lit(sample_name_truth).alias('sample_name'),
        )
        data_sensitivity = data_sensitivity.select(['sample_name'] + data_sensitivity.columns[:-1])

        truth_vs_prediction.append(data_sensitivity)

    prediction_vs_truth = pl.concat(prediction_vs_truth, how="vertical_relaxed")
    truth_vs_prediction = pl.concat(truth_vs_prediction, how="vertical_relaxed")

    prediction_vs_truth.sort(['sample_name', 'chrom', 'prediction_start']).write_csv(
        output_dir / "prediction_vs_truth.tsv", include_header=True, separator="\t"
    )

    truth_vs_prediction.sort(['sample_name', 'chrom', 'truth_start']).write_csv(
        output_dir / "truth_vs_prediction.tsv", include_header=True, separator="\t"
    )

    prediction_vs_truth = pl.read_csv(
        output_dir / "prediction_vs_truth.tsv", has_header=True, separator="\t",
        infer_schema = False,
        
    )
    truth_vs_prediction = pl.read_csv(
        output_dir / "truth_vs_prediction.tsv", has_header=True, separator="\t",
        infer_schema = False,
    )

    prediction_vs_truth = prediction_vs_truth.with_columns(
        pl.col("prediction_start").cast(pl.Int64).alias("start"),
        pl.col("prediction_end").cast(pl.Int64).alias("end"),
        pl.col("prediction_length").cast(pl.Int64),
        # pl.col("prediction_n_markers").cast(pl.Int64),
        pl.col("reciprocal_overlap").cast(pl.Float64),
    )

    truth_vs_prediction = truth_vs_prediction.with_columns(
        pl.col("truth_start").cast(pl.Int64).alias("start"),
        pl.col("truth_end").cast(pl.Int64).alias("end"),
        pl.col("truth_length").cast(pl.Int64),
        # pl.col("truth_n_markers").cast(pl.Int64),
        pl.col("reciprocal_overlap").cast(pl.Float64),
    )


    report(prediction_vs_truth, output_dir / 'prediction_fragments.tsv', output_dir , 'ppv')
    report(truth_vs_prediction, output_dir / 'truth_fragments.tsv', output_dir , 'sensitivity')

    report_by_size(prediction_vs_truth, output_dir / 'ppv-markers_vs_len.tsv', 'prediction')
    report_by_size(truth_vs_prediction, output_dir / 'sensitivity-markers_vs_len.tsv', 'truth')




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


def summarize_sliding(
    data,
    chrom_length_file,
    complex_regions_file,
    window_size,
    step_size,
    reciprocal_overlap_cutoff,
    n_threads,
):

    def jobs(
        data,
        chrom_orders,
        chrom_lengths,
        complex_regions,
        window_size,
        step_size,
        reciprocal_overlap_cutoff,
    ):

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

    def process(job):
        chrom = job["chrom"]
        chrom_length = job["chrom_length"]
        data = job["data"]
        window_size = job["window_size"]
        step_size = job["step_size"]
        reciprocal_overlap_cutoff = job["reciprocal_overlap_cutoff"]
        complex_regions = job["complex_regions"]

        database = create_database(data)
        complex_region_db = create_database(complex_regions)

        bag = list()

        for start in range(0, chrom_length, step_size):

            end = min(start + window_size, chrom_length)

            query = GenomicRegion(chrom, start, end)

            matches = complex_region_db.find_overlap(
                {"chrom": query.chrom, "start": query.start, "end": query.end}
            )

            if len(matches) > 0:
                continue

            matches = database.find_overlap(
                {"chrom": query.chrom, "start": query.start, "end": query.end}
            )

            if len(matches) == 0:
                continue

            n_detected = 0
            n_matched = 0

            for match in matches:

                n_detected += 1

                if match["matched"] is None:
                    continue
                if match["matched"] is False:
                    continue

                n_matched += 1
            if n_detected == 0:
                concordance = None
            else:
                concordance = n_matched / n_detected

            bag.append(
                {
                    "chrom": query.chrom,
                    "start": query.start,
                    "end": query.end,
                    "concordance": concordance,
                    "n_detected": n_detected,
                    "n_matched": n_matched,
                }
            )
        return bag

    bag = list()

    chrom_lengths = load_chrom_lengths(chrom_length_file)
    chrom_orders = sorted(chrom_lengths.keys(), key=chrom_sort_key)
    complex_regions = load_complex_regions(complex_regions_file)

    # for job in jobs(data, chrom_orders, chrom_lengths, window_size, step_size, reciprocal_overlap_cutoff):
    #     result = process(job)
    #     bag.extend(result)

    with ProcessPool(n_threads) as pool:
        for result in pool.uimap(
            process,
            jobs(
                data,
                chrom_orders=chrom_orders,
                chrom_lengths=chrom_lengths,
                complex_regions=complex_regions,
                window_size=window_size,
                step_size=step_size,
                reciprocal_overlap_cutoff=reciprocal_overlap_cutoff,
            ),
        ):
            bag.extend(result)

    result = pl.from_dicts(bag, infer_schema_length=None)

    return result


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


def _validate_cnv(
    queries,
    database,
    qname,
    dbname,
    reciprocal_overlap_cutoff,
    breakpoint_tolerance_cutoff,
):
    bag = list()
    for query in queries:

        matches = database.find_overlap(query)
        query_region = GenomicRegion(query["chrom"], query["start"], query["end"])

        if len(matches) == 0:
            result = dict()
            result[f"{qname}_fragment_idx"] = query["fragment_idx"]
            result[f"chrom"] = query["chrom"]
            result[f"{qname}_start"] = query["start"]
            result[f"{qname}_end"] = query["end"]
            result[f"{qname}_length"] = query["end"] - query["start"]
            result[f"{qname}_cn_state"] = query["cn_state"]
            # result[f"{qname}_n_markers"] = query["n_markers"] if 'n_markers' in query else None
            result[f"{dbname}_fragment_idx"] = None
            result[f"{dbname}_start"] = None
            result[f"{dbname}_end"] = None
            result[f"{dbname}_length"] = None
            result[f"{dbname}_cn_state"] = None
            # result[f"{dbname}_n_markers"] = None
            result["reciprocal_overlap"] = None
            result["breakpoint_difference"] = None
            result["reciprocal_overlap_test"] = None
            result["breakpoint_tolerance_test"] = None
            result["cn_state_test"] = None
            result["matched"] = None
            bag.append(result)
            continue

        for match in matches:
            result = dict()
            db_region = GenomicRegion(match["chrom"], match["start"], match["end"])
            reciprocal_overlap = query_region.get_reciprocal_overlap(db_region)['ratio']
            breakpoint_difference = query_region.get_max_boundary_difference(db_region)

            result[f"{qname}_fragment_idx"] = query["fragment_idx"]
            result[f"chrom"] = query["chrom"]
            result[f"{qname}_start"] = query["start"]
            result[f"{qname}_end"] = query["end"]
            result[f"{qname}_length"] = query["end"] - query["start"]
            result[f"{qname}_cn_state"] = query["cn_state"]
            # result[f"{qname}_n_markers"] = query["n_markers"] if 'n_markers' in query else None
            result[f"{dbname}_fragment_idx"] = match["fragment_idx"]
            result[f"{dbname}_start"] = match["start"]
            result[f"{dbname}_end"] = match["end"]
            result[f"{dbname}_length"] = match["end"] - match["start"]
            result[f"{dbname}_cn_state"] = match["cn_state"]
            # result[f"{dbname}_n_markers"] = match["n_markers"] if 'n_markers' in match else None
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

            if result["reciprocal_overlap_test"] == "FAIL":
                result["matched"] = False
            # elif result["breakpoint_tolerance_test"] == "FAIL":
            #     result["matched"] = False
            elif result["cn_state_test"] == "FAIL":
                result["matched"] = False
            else:
                result["matched"] = True

            bag.append(result)

    report = pl.from_dicts(bag, infer_schema_length=None)


    return report


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
            on_hotspot = data.filter(pl.col('hotspot_overlap_ratio').cast(pl.Float64) > 0)

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

def report(a_vs_b, fragments_file, output_dir, _type):

    fragments = pl.read_csv(
            fragments_file, has_header = True, separator = '\t', 
            infer_schema = False,)
            # schema_overrides = {'barcode': pl.Utf8, 'fragment_idx': pl.Utf8})
    if _type == 'ppv':
        frag_idx = 'prediction_fragment_idx'
    elif _type == 'sensitivity':
        frag_idx = 'truth_fragment_idx'


    data = a_vs_b.rename({frag_idx: 'fragment_idx'})

    # x = data.unique(subset = ['fragment_idx'])
    #
    # duplicates = x.filter(x.select(['sample_name', 'chrom', 'prediction_start', 'prediction_end']).is_duplicated())\
    #         .sort(['sample_name', 'chrom', 'prediction_start'])
    #
    # print('#' * 50)
    # print(duplicates)
    #
    # print('unique fragments' , len(data['fragment_idx'].unique()))
    # print('unique record', len(data.unique(subset=['sample_name', 'chrom', 'prediction_start', 'prediction_end'])))
    # print('unique record', len(data.unique(subset=['sample_name', 'chrom', 'truth_start', 'truth_end'])))


    # data = fragments.join(a_vs_b, left_on = 'fragment_idx', right_on = frag_idx, how = 'left')
    #
    # print(_type, len(data.filter(pl.col('matched') == True)))

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

    print(_type, len(data.filter(pl.col('matched') == True)))


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

def report_by_size(data, output_file, label):

    if label == 'prediction': 
        metric = 'ppv'
    elif label == 'truth':
        metric = 'sensitivity'

    lengths =  [0, 50000,100000,250000, 500000,1000000, np.inf]
    n_markers = [0, 50, 100, 500, 1000, 2000, 4000, 8000, 10000, np.inf]

    bag = list()

    for i in range(0, len(lengths) - 1):
        min_length = lengths[i]
        max_length = lengths[i + 1]
        target = data.filter(
            (pl.col(f'{label}_length') >= min_length) & 
            (pl.col(f'{label}_length') < max_length)
        )

        for j in range(0, len(n_markers) - 1):

            min_n_markers = n_markers[j]
            max_n_markers = n_markers[j + 1]

            target2 = target.filter(
                (pl.col(f'{label}_n_markers') >= min_n_markers) & 
                (pl.col(f'{label}_n_markers') < max_n_markers)
            )

            n_detected = len(target2)
            n_matched = len(target2.filter(pl.col('matched') == True))

            if n_detected == 0:
                value = None
            else:
                value = n_matched / n_detected

            bag.append({
                'min_length': min_length,
                'max_length': max_length,
                'min_n_markers': min_n_markers,
                'max_n_markers': max_n_markers,
                'n_detected': n_detected,
                'n_matched': n_matched,
                f'{metric}': value,
            })

    data = pl.from_dicts(bag)
    data.write_csv(output_file, include_header = True, separator = '\t')

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


def annotate_hotspot_regions(cnvs, hotspot_region_db):
    cnvs = cnvs.with_columns(
        pl.col("start").cast(pl.Int64).alias("start"),
        pl.col("end").cast(pl.Int64).alias("end"),
    )


    bag = list()

    for record in cnvs.to_dicts():


        # if record['sample_name'] == '20200915_GT3_ARRAY_435_H02_TPMI.CEL'and record['start_position'] == '18967500':

        matches = hotspot_region_db.find_overlap(record)

        tmp_record = deepcopy(record)


        if len(matches) == 0:
            overlap_ratio = 0
            hotspot_name = None
        else:

            x = GenomicRegion(record['chrom'], record['start'], record['end'])

            hotspot_names = list()

            overlap_length = 0

            for match in matches:
                intersect = x.intersects(GenomicRegion(match['chrom'], int(match['start']), int(match['end'])))


                overlap_length += len(intersect)
                hotspot_names.append(match['region_name'])

            overlap_ratio = overlap_length / len(x)
            hotspot_name = ','.join(hotspot_names)

        tmp_record['hotspot_overlap_ratio'] = overlap_ratio
        tmp_record['hotspot_name'] = hotspot_name
        bag.append(tmp_record)

    return pl.from_dicts(bag, infer_schema_length = None)


def annotate_complex_regions(cnvs, complex_region_db):

    cnvs = cnvs.with_columns(
        pl.col("start").cast(pl.Int64).alias("start"),
       pl.col("end").cast(pl.Int64).alias("end"),
    )

    bag_region = list()
    bag_fragment = list()


    for record in cnvs.to_dicts():


        matches = complex_region_db.find_overlap(record)

        if len(matches) == 0:
            region = record.copy()
            fragment = record.copy()

            region["complex_region"] = None
            bag_region.append(region)

            fragment["fragment_idx"] = fragment["region_idx"]
            bag_fragment.append(fragment)
        else:
            matches = sort_matches(matches)

            region = deepcopy(record)
            region['complex_region'] = ';'.join([x['category'] for x in matches])
            bag_region.append(region)

            fragments = chop_region(record, matches)

            bag_fragment.extend(fragments)


    regions = pl.from_dicts(bag_region, infer_schema_length=None)
    fragments = pl.from_dicts(bag_fragment, infer_schema_length=None)

    return regions, fragments


def chop_region(record, matches):

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

