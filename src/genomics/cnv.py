from importlib import resources

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


def validate(
    predictions_file,
    truths_file,
    sample_map_file,
    output_dir,
    cnvmix_regions_file=None,
    reciprocal_overlap_cutoff=0.5,
    boundary_difference_cutoff=10000,
    window_size=10000000,
    step_size=1000,
    concordance_cutoff=0.8,
    n_threads=32,
):

    predictions = pl.read_csv(
        predictions_file, has_header=True, separator="\t", infer_schema=False
    )
    truths = pl.read_csv(
        truths_file, has_header=True, separator="\t", infer_schema=False
    )

    sample_map = load_sample_map(sample_map_file)

    assert REQUIRED_COLUMNS.issubset(predictions.columns), predictions.columns
    assert REQUIRED_COLUMNS.issubset(truths.columns), truths.columns

    predictions = predictions.with_columns(
        pl.col("start").cast(pl.Int64),
        pl.col("end").cast(pl.Int64),
        pl.arange(0, len(predictions)).alias("region_idx"),
    )
    truths = truths.with_columns(
        pl.col("start").cast(pl.Int64),
        pl.col("end").cast(pl.Int64),
        pl.arange(0, len(truths)).alias("region_idx"),
    )

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
        sample_name = sample["sample_name_prediction"]
        prediction_fragments = sample["prediction_fragments"].to_dicts()
        truth_fragments = sample["truth_fragments"].to_dicts()

        data_ppv = _validate_cnv(
            prediction_fragments,
            create_database(truth_fragments),
            "prediction",
            "truth",
            reciprocal_overlap_cutoff,
            boundary_difference_cutoff,
        )

        prediction_vs_truth.append(data_ppv)

        data_sensitivity = _validate_cnv(
            truth_fragments,
            create_database(prediction_fragments),
            "truth",
            "prediction",
            reciprocal_overlap_cutoff,
            boundary_difference_cutoff,
        )

        truth_vs_prediction.append(data_sensitivity)

    prediction_vs_truth = pl.concat(prediction_vs_truth, how="vertical_relaxed")
    truth_vs_prediction = pl.concat(truth_vs_prediction, how="vertical_relaxed")

    prediction_vs_truth.write_csv(
        output_dir / "prediction_vs_truth.tsv", include_header=True, separator="\t"
    )
    truth_vs_prediction.write_csv(
        output_dir / "truth_vs_prediction.tsv", include_header=True, separator="\t"
    )

    prediction_vs_truth = pl.read_csv(
        output_dir / "prediction_vs_truth.tsv",
        has_header=True,
        separator="\t",
        infer_schema=False,
    )
    truth_vs_prediction = pl.read_csv(
        output_dir / "truth_vs_prediction.tsv",
        has_header=True,
        separator="\t",
        infer_schema=False,
    )
    prediction_vs_truth = prediction_vs_truth.with_columns(
        pl.col("prediction_start").cast(pl.Int64).alias("start"),
        pl.col("prediction_end").cast(pl.Int64).alias("end"),
        pl.col("reciprocal_overlap").cast(pl.Float64),
    )

    truth_vs_prediction = truth_vs_prediction.with_columns(
        pl.col("truth_start").cast(pl.Int64).alias("start"),
        pl.col("truth_end").cast(pl.Int64).alias("end"),
        pl.col("reciprocal_overlap").cast(pl.Float64),
    )

    ppv_summary = summarize_sliding(
        data=prediction_vs_truth,
        chrom_length_file=CHROM_LENGTH_HG38_FILE,
        complex_regions_file=COMPLEX_REGIONS_FILE,
        # bin_size=bin_size,
        window_size=window_size,
        step_size=step_size,
        reciprocal_overlap_cutoff=reciprocal_overlap_cutoff,
        n_threads=n_threads,
    )

    sensitivity_summary = summarize_sliding(
        data=truth_vs_prediction,
        chrom_length_file=CHROM_LENGTH_HG38_FILE,
        complex_regions_file=COMPLEX_REGIONS_FILE,
        # bin_size=bin_size,
        window_size=window_size,
        step_size=step_size,
        reciprocal_overlap_cutoff=reciprocal_overlap_cutoff,
        n_threads=n_threads,
    )

    ppv_summary.write_csv(
        output_dir / "ppv_summary.tsv", include_header=True, separator="\t"
    )
    sensitivity_summary.write_csv(
        output_dir / "sensitivity_summary.tsv",
        include_header=True,
        separator="\t",
    )

    ppv_summary = pl.read_csv(
        output_dir / "ppv_summary.tsv",
        has_header=True,
        separator="\t",
        schema_overrides={"concordance": pl.Float64},
    )

    sensitivity_summary = pl.read_csv(
        output_dir / "sensitivity_summary.tsv",
        has_header=True,
        separator="\t",
        schema_overrides={"concordance": pl.Float64},
    )

    ppv_high_concordance_regions = create_high_concordance_region(
        ppv_summary, concordance_cutoff
    )
    ppv_high_concordance_regions.write_csv(
        output_dir / f"regions_ppv-{concordance_cutoff}.tsv",
        include_header=True,
        separator="\t",
    )

    kp.plot_karyopype(
        [
            COMPLEX_REGIONS_FILE,
            output_dir / f"regions_ppv-{concordance_cutoff}.tsv",
        ],
        ["complex_regions", f"concordance_{concordance_cutoff}"],
        output_dir / f"regions_ppv-{concordance_cutoff}.png",
    )

    sensitivity_high_concordance_regions = create_high_concordance_region(
        sensitivity_summary, concordance_cutoff
    )

    sensitivity_high_concordance_regions.write_csv(
        output_dir / f"regions_sensitivity-{concordance_cutoff}.tsv",
        include_header=True,
        separator="\t",
    )

    kp.plot_karyopype(
        [
            COMPLEX_REGIONS_FILE,
            output_dir / f"regions_sensitivity-{concordance_cutoff}.tsv",
        ],
        ["complex_regions", f"concordance_{concordance_cutoff}"],
        output_dir / f"regions_sensitivity-{concordance_cutoff}.png",
    )


def create_high_concordance_region(moving_average_data, concordance_cutoff):

    # print(moving_average_data["concordance"].max())

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
            n_same = 0

            for match in matches:

                n_detected += 1

                if match["reciprocal_overlap"] is None:
                    continue
                if match["reciprocal_overlap_test"] == "FAIL":
                    continue
                # if match["breakpoint_tolerance_test"] == "FAIL":
                #     continue
                if match["cn_state_test"] == "FAIL":
                    continue

                n_same += 1
            if n_detected == 0:
                concordance = None
            else:
                concordance = n_same / n_detected

            bag.append(
                {
                    "chrom": query.chrom,
                    "start": query.start,
                    "end": query.end,
                    "concordance": concordance,
                    "n_detected": n_detected,
                    "n_same": n_same,
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
    boundary_difference_cutoff,
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
            result[f"{dbname}_fragment_idx"] = None
            result[f"{dbname}_start"] = None
            result[f"{dbname}_end"] = None
            result[f"{dbname}_length"] = None
            result[f"{dbname}_cn_state"] = None
            result["reciprocal_overlap"] = None
            result["boundary_difference"] = None
            result["reciprocal_overlap_test"] = None
            result["breakpoint_tolerance_test"] = None
            result["cn_state_test"] = None
            bag.append(result)
            continue

        for match in matches:
            result = dict()
            db_region = GenomicRegion(match["chrom"], match["start"], match["end"])
            reciprocal_overlap = query_region.get_reciprocal_overlap(db_region)
            boundary_difference = query_region.get_max_boundary_difference(db_region)

            result[f"{qname}_fragment_idx"] = query["fragment_idx"]
            result[f"chrom"] = query["chrom"]
            result[f"{qname}_start"] = query["start"]
            result[f"{qname}_end"] = query["end"]
            result[f"{qname}_length"] = query["end"] - query["start"]
            result[f"{qname}_cn_state"] = query["cn_state"]
            result[f"{dbname}_fragment_idx"] = match["fragment_idx"]
            result[f"{dbname}_start"] = match["start"]
            result[f"{dbname}_end"] = match["end"]
            result[f"{dbname}_length"] = match["end"] - match["start"]
            result[f"{dbname}_cn_state"] = match["cn_state"]
            result["reciprocal_overlap"] = reciprocal_overlap
            result["boundary_difference"] = boundary_difference

            result["reciprocal_overlap_test"] = (
                "PASS" if reciprocal_overlap >= reciprocal_overlap_cutoff else "FAIL"
            )
            result["breakpoint_tolerance_test"] = (
                "PASS" if boundary_difference <= boundary_difference_cutoff else "FAIL"
            )
            result["cn_state_test"] = (
                "PASS" if query["cn_state"] == match["cn_state"] else "FAIL"
            )

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
        if sample_name_prediction not in sample_map:
            continue

        sample_name_truth = sample_map[sample_name_prediction]

        bag.append(
            {
                "sample_name_prediction": sample_name_prediction,
                "sample_name_truth": sample_name_truth,
                "prediction_fragments": df,
                "truth_fragments": truth_bag[sample_name_truth],
            }
        )

    return bag


def load_complex_regions(complex_regions_file):

    data = pl.read_csv(complex_regions_file, has_header=True, separator="\t").to_dicts()

    return data


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
            record["complex_region"] = None
            bag_region.append(record)

            fragment = record.copy()
            fragment["fragment_idx"] = record["region_idx"]
            bag_fragment.append(fragment)
        else:
            region_idx = record["region_idx"]
            frag_idx = 0

            for match in matches:
                query = record.copy()
                match_region = GenomicRegion(
                    match["chrom"], match["start"], match["end"]
                )
                query_region = GenomicRegion(
                    query["chrom"], query["start"], query["end"]
                )

                query["complex_region"] = match["category"]
                query["reciprocal_overlap"] = query_region.get_reciprocal_overlap(
                    match_region
                )

                bag_region.append(query)

                fragment = record.copy()
                if record["start"] < match["start"]:
                    fragment["start"] = record["start"]
                    fragment["end"] = match["start"]
                    fragment["fragment_idx"] = f"{region_idx}_{frag_idx}"
                    bag_fragment.append(fragment)
                    frag_idx += 1

                if record["end"] > match["end"]:
                    fragment["start"] = match["end"]
                    fragment["end"] = record["end"]
                    fragment["fragment_idx"] = f"{region_idx}_{frag_idx}"
                    bag_fragment.append(fragment)
                    frag_idx += 1
    regions = pl.from_dicts(bag_region, infer_schema_length=None)
    fragments = pl.from_dicts(bag_fragment, infer_schema_length=None)

    return regions, fragments
