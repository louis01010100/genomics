import re
from importlib import resources
from pathlib import Path

import matplotlib.pyplot as plt
import polars as pl
from matplotlib.collections import PolyCollection
from matplotlib.lines import Line2D

CHROM_LENGTHS_FILE = (
    resources.files("genomics") / "resources" / "chrom_lengths-hg38.tsv"
)


def get_chromsizes(species, chromsizes=None, cannonical=True):
    bag = dict()
    with CHROM_LENGTHS_FILE.open("rt") as fh:
        next(fh)
        for line in fh:
            items = line.strip().split("\t")
            bag[items[0]] = int(items[1])
    return bag


def parse_regions(region_file):
    df = pl.read_csv(region_file, separator="\t", has_header=True)
    df = df.select(df.columns[:3])
    df.columns = ["chrom", "start", "end"]
    df = df.slice(1)
    df = df.with_columns(pl.col("start").cast(pl.Int64), pl.col("end").cast(pl.Int64))
    return df


def add_chromsize_colors(sizes_dict, color):
    return {k: (v, color) for k, v in sizes_dict.items()}


def chrom_sort_key(chrom):
    match = re.match(r"chr(\d+|X|Y|M)$", chrom)
    if not match:
        return (1000,)
    val = match.group(1)
    if val.isdigit():
        return (int(val),)
    elif val == "X":
        return (23,)
    elif val == "Y":
        return (24,)
    elif val == "M":
        return (25,)
    else:
        return (999,)


def chromosome_collections(df: pl.DataFrame, y_positions, height, **kwargs):
    if "width" not in df.columns:
        df = df.with_columns((pl.col("end") - pl.col("start")).alias("width"))

    for chrom in df["chrom"].unique().to_list():
        group = df.filter(pl.col("chrom") == chrom)
        ybase = y_positions[chrom]
        segments = []
        colors = []

        for row in group.iter_rows(named=True):
            segments.append(
                [
                    [row["start"], ybase],
                    [row["start"] + row["width"], ybase],
                    [row["start"] + row["width"], ybase + height],
                    [row["start"], ybase + height],
                ]
            )
            colors.append(row["colors"])

        poly = PolyCollection(segments, facecolors=colors, **kwargs)
        yield poly


def add_regions(ax, chromsizes, region_files):
    chromsize_colors = add_chromsize_colors(chromsizes, color="whitesmoke")
    chrom_height = 1
    chrom_spacing = 1
    gene_height = 0.4
    gene_padding = 0.1
    chromosome_list = sorted(chromsizes.keys(), key=chrom_sort_key)

    ybase = 0
    chrom_ybase = {}
    gene_ybase = {}
    chrom_centers = {}

    for chrom in chromosome_list[::-1]:
        chrom_ybase[chrom] = ybase
        chrom_centers[chrom] = ybase + chrom_height / 2.0
        gene_ybase[chrom] = ybase - gene_height - gene_padding
        ybase += chrom_height + chrom_spacing

    ideo = pl.DataFrame(
        {
            "chrom": list(chromsize_colors.keys()),
            "end": [v[0] for v in chromsize_colors.values()],
            "colors": [v[1] for v in chromsize_colors.values()],
        }
    ).with_columns(pl.lit(0).alias("start"))
    ideo = ideo.with_columns((pl.col("end") - pl.col("start")).alias("width"))

    for collection in chromosome_collections(
        ideo, chrom_ybase, chrom_height, edgecolor="k"
    ):
        ax.add_collection(collection)

    leg = []
    leg_lab = []

    for i, r in enumerate(region_files):
        print(i, r)
        color = f"C{i}"
        leg_lab.append(f"regionss_file{i+1}")
        leg.append(Line2D([0], [0], color=color, lw=4))
        rdf = parse_regions(r)

        rdf = rdf.with_columns(pl.lit(color).alias("colors"))
        for collection in chromosome_collections(
            rdf, chrom_ybase, chrom_height, edgecolor=color
        ):
            ax.add_collection(collection)

    ax.set_yticks([chrom_centers[i] for i in chromosome_list])
    ax.set_yticklabels(chromosome_list)
    ax.legend(leg, leg_lab, loc=4)
    ax.axis("tight")
    return ax


def plot_karyopype(region_files, output_file, figsize=(10, 7)):
    chromsizes = get_chromsizes(CHROM_LENGTHS_FILE)
    fig, ax = plt.subplots(1, 1, figsize=figsize)

    add_regions(ax=ax, chromsizes=chromsizes, region_files=region_files)

    ax.xaxis.set_major_formatter(plt.FuncFormatter(lambda x, n: int(x / 1e6)))
    ax.set_title(f"Human GRCh38 Karyopype", fontsize=14)
    plt.xlabel("Chromosome Position (Mbp)", fontsize=14)
    plt.ylabel(f"Chromosome", fontsize=14)
    plt.savefig(output_file)
