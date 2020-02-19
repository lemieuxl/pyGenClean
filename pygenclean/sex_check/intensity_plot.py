"""Creates an intensity plot for sexual chromosomes."""


import logging
import argparse
from os import path

import pandas as pd

import matplotlib.pyplot as plt

from ..utils import decode_chrom
from ..utils.plot import generate_html_scatter
from ..utils.plink import get_markers_on_chrom, get_sample_sexes

from ..error import ProgramError

from ..version import pygenclean_version as __version__


SCRIPT_NAME = "intensity-plot"
DESCRIPTION = "Plots summarized X intensities versus summarized Y intensities."


logger = logging.getLogger(__name__)


def main(args=None, argv=None):
    """Creates an intensity plot for sexual chromosomes.

    Args:
        args (argparse.Namespace): the arguments and options.
        argv (list): the arguments as list.

    These are the steps:

    1. If there are ``summarized_intensities`` provided, reads the file and
       skips to step 6.
    2. Gets markers on the sexual chromosomes
       (:py:func:`get_markers_on_chrom`).
    3. Gets sex for each sample (:py:func:`get_sample_sexes`).
    4. Reads the file containing samples with sex mismatches
       (:py:func:`read_sex_mismatches`).
    5. Reads the intensities and summarizes them (:py:func:`read_intensities`).
    6. Plots the summarized intensities
       (:py:func:`plot_summarized_intensities`).

    Note:
        Since the intensity file only contain a single sample ID (instead of
        the traditional ``FID`` and ``IID``, only the ``IID`` is kept.

    """
    if args is None:
        args = parse_args(argv)
    check_args(args)

    data = None
    if args.summarized_intensities is None:
        # Retrieving the markers on the sexual chromosomes
        logger.info("Finding markers on chromosomes 23 and 24")
        xy_markers = get_markers_on_chrom(args.bfile + ".bim", {23, 24})

        # Reading the FAM file
        logger.info("Fetching sample's sex")
        sample_sex = get_sample_sexes(args.bfile + ".fam", only_iid=True)

        # Reading the sex problem file
        logger.info("Fetching samples with sex mismatches")
        sex_mismatches = read_sex_mismatches(args.sex_mismatches)

        # Reading the intensity file
        logger.info("Fetching intensities")
        data = read_intensities(
            filename=args.intensities,
            required_markers=xy_markers,
            sample_sex=sample_sex,
            mismatches=sex_mismatches,
            args=args,
        )

        # Saving the summarized intensities
        data.to_csv(f"{args.out}.summarized_intensities.tsv", sep="\t",
                    index=False)

    else:
        logger.info("Retrieving summarized intensities from %s",
                    args.summarized_intensities)
        data = pd.read_csv(args.summarized_intensities, sep="\t")

    # Plotting
    logger.info("Generating the intensity plot")
    plot_summarized_intensities(data, args)


def read_intensities(filename, required_markers, sample_sex, mismatches, args):
    """Reads the intensities from a file.

    Args:
        filename (str): the name of the input file.
        required_markers (set): the set of required markers.
        sample_sex (dict): the sex of each sample.
        mismatches (set): the samples with issues.

    Returns:
        pd.DataFrame: the summarized data containing the following information:
        sample_id, chr23, chr24, sex and status.

    Reads the normalized intensities from a CSV file, keeping only the required
    markers (*i.e.* those on chromosomes 23 and 24).

    The final dara set contains the following information:

    - ``sample_id``: the sample ID.
    - ``chr23``: the summarized intensities for chromosome 23.
    - ``chr24``: the summarized intensities for chromosome 24.
    - ``sex``: the sex of the sample (``Male``, ``Female`` or ``Unknown``).
    - ``status``: the status of the sample (``OK`` or ``Mismatch``).

    The summarized intensities for a chromosome (:math:`S_{chr}`) is computed
    using this formula (where :math:`I_{chr}` is the set of all marker
    intensities on chromosome :math:`chr`):

    .. math::
        S_{chr} = \\frac{\\sum{I_{chr}}}{||I_{chr}||}

    """
    # The columns
    cols = {
        "snp": args.intensity_snp_col,
        "sample": args.intensity_sample_col,
        "chrom": args.intensity_chrom_col,
        "x": args.intensity_x_col,
        "y": args.intensity_y_col,
    }

    df = None
    if args.intensity_nb_lines == "all":
        df = pd.read_csv(
            filename,
            sep=args.intensity_delimiter,
            dtype={cols["chrom"]: str, cols["sample"]: str},
            usecols=[cols["snp"], cols["sample"], cols["x"], cols["y"],
                     cols["chrom"]],
        )
        df = process_df(
            df=df, snp_col=cols["snp"], sample_col=cols["sample"],
            chrom_col=cols["chrom"], sample_sex=sample_sex,
            required_markers=required_markers,
        )

    else:
        # Will read chunk by chunk using pandas in case we have a huge file
        reader = pd.read_csv(
            filename,
            sep=args.intensity_delimiter,
            dtype={cols["chrom"]: str, cols["sample"]: str},
            usecols=[cols["snp"], cols["sample"], cols["x"], cols["y"],
                     cols["chrom"]],
            chunksize=args.intensity_nb_lines,
        )

        # Reading the chunks and concatenating the values
        df = pd.concat(
            [
                process_df(
                    df=sub_df, snp_col=cols["snp"], sample_col=cols["sample"],
                    chrom_col=cols["chrom"], sample_sex=sample_sex,
                    required_markers=required_markers,
                ) for sub_df in reader
            ],
            ignore_index=True,
        )

    logger.info("Computing summarized intensities")
    # Sum X and Y, group by sample and chromosomes and average the sum
    df["sum_x_y"] = df.loc[:, [cols["x"], cols["y"]]].sum(axis=1)
    df = df.groupby([cols["sample"], "__chrom"]).mean().reset_index().pivot(
        index=cols["sample"], columns="__chrom", values="sum_x_y",
    )

    # Renaming the columns
    df = df.reset_index().rename(
        columns={cols["sample"]: "sample_id", 23: "chr23", 24: "chr24"},
    )

    # Adding the required information (sex and status)
    df["sex"] = [sample_sex[sample] for sample in df.sample_id]
    df["status"] = [
        "Mismatch" if sample in mismatches else "OK" for sample in df.sample_id
    ]

    return df.loc[:, ["sample_id", "chr23", "chr24", "sex", "status"]]


def process_df(df, snp_col, sample_col, chrom_col, sample_sex,
               required_markers):
    """Pre-process a data frame to keep only required information.

    Args:
        df (pandas.DataFrame): the data.
        snp_col (str): the name of the marker column.
        sample_col (str): the name of the sample column.
        chrom_col (str): the name of the chromosome column.
        sample_sex (dict): the sex of each sample.
        required_markers (set): the set of required markers.

    Returns:
        pandas.DataFrame: the processed data.

    """
    # Keeping only the required markers and samples
    df = df.loc[df.loc[:, snp_col].isin(required_markers), :]
    df = df.loc[df.loc[:, sample_col].isin(sample_sex), :]

    # Decoding the chromosome
    df["__chrom"] = df.loc[:, chrom_col].map(decode_chrom)

    return df.dropna()


def read_sex_mismatches(filename):
    """Reads the samples with sex mismatches.

    Args:
        filename (str): the name of the file.

    Returns:
        set: samples with sex issues.

    Note:
        If ``filename`` is ``None``, an empty set is returned.

    """
    mismatches = set()

    if filename is None:
        return mismatches

    header = None
    with open(filename) as f:
        for line in f:
            row = line.rstrip("\r\n").split()

            if header is None:
                header = {name: i for i, name in enumerate(row)}
                continue

            mismatches.add(row[header["IID"]])

    return mismatches


def plot_summarized_intensities(df, args):
    """Plots the summarized intensities.

    Args:
        df (pandas.DataFrame): the summarized_intensities.
        args (argparse.Namespace): the options

    Plots the summarized intensities of the markers on the sexual chromosomes.
    Samples with a sex mismatch will be highlighted.

    """
    figure, axe = plt.subplots(1, 1)

    # Changing the spines
    axe.xaxis.set_ticks_position("bottom")
    axe.yaxis.set_ticks_position("left")
    axe.spines["top"].set_visible(False)
    axe.spines["right"].set_visible(False)

    # Changing the titles
    axe.set_xlabel(args.xlabel)
    axe.set_ylabel(args.ylabel)

    # The colors
    plot_config = {
        ("OK", "Male"): {
            "color": "#0099CC", "edgecolor": "#0099CC", "marker": "o",
            "size": 20, "plotly_size": 12, "plotly_symbol": "circle",
        },
        ("OK", "Female"): {
            "color": "#CC0000", "edgecolor": "#CC0000", "marker": "o",
            "size": 20, "plotly_size": 12, "plotly_symbol": "circle",
        },
        ("OK", "Unknown"): {
            "color": "#555555", "edgecolor": "#555555", "marker": "o",
            "size": 20, "plotly_size": 12, "plotly_symbol": "circle",
        },
        ("Mismatch", "Male"): {
            "color": "#669900", "edgecolor": "#000000", "marker": "^",
            "size": 30, "plotly_size": 12, "plotly_symbol": "triangle-up",
            "plotly_line": {"color": "black", "width": 1},
        },
        ("Mismatch", "Female"): {
            "color": "#9933CC", "edgecolor": "#000000", "marker": "v",
            "size": 30, "plotly_size": 12, "plotly_symbol": "triangle-down",
            "plotly_line": {"color": "black", "width": 1},
        },
        ("Mismatch", "Unknown"): {
            "color": "#555555", "edgecolor": "#000000", "marker": ">",
            "size": 30, "plotly_size": 12, "plotly_symbol": "triangle-right",
            "plotly_line": {"color": "black", "width": 1},
        },
    }

    html_data = {
        "page_title": "pyGenClean - sex-check - Summarized Intensities",
        "scatter_title": "Summarized Intensities",
        "xlabel": args.xlabel,
        "ylabel": args.ylabel,
        "labels": [],
        "data": {},
        "colors": {},
        "sizes": {},
        "symbols": {},
        "lines": {},
    }
    for status in ("OK", "Mismatch"):
        for sex in ("Male", "Female", "Unknown"):
            subset = (df.sex == sex) & (df.status == status)
            nb_subset = subset.sum()
            sub_df = df.loc[subset, :]

            config = plot_config[(status, sex)]
            label = f"{status} {sex} (n={nb_subset:,d})"
            axe.scatter(
                sub_df.chr23, sub_df.chr24, s=config["size"],
                color=config["color"], edgecolors=config["edgecolor"],
                marker=config["marker"], label=label,
            )

            # Saving the HTML data
            html_data["labels"].append(label)
            html_data["data"][label] = {
                "x": sub_df.chr23.tolist(),
                "y": sub_df.chr24.tolist(),
                "ids": sub_df.sample_id.tolist(),
            }
            html_data["symbols"][label] = config["plotly_symbol"]
            html_data["colors"][label] = config["color"]
            html_data["sizes"][label] = config["plotly_size"]
            if "plotly_line" in config:
                html_data["lines"][label] = config["plotly_line"]

    axe.legend(loc=8, fancybox=True, ncol=2, borderaxespad=0,
               bbox_to_anchor=(0, 1.02, 1, 0.102))

    plt.savefig(f"{args.out}.{args.format}", bbox_inches="tight", dpi=args.dpi)
    plt.close(figure)

    # Creating the HTML file
    generate_html_scatter(f"{args.out}.html", **html_data)


def check_args(args):
    """Checks the arguments and options.

    Args:
        args (argparse.Namespace): the arguments and options to check.

    """
    # Checking if there are input files
    if args.summarized_intensities is None:
        if args.bfile is None or args.intensities is None:
            raise ProgramError(
                "need to specify either '--bfile' and '--intensities', or "
                "'--summarized-intensities'"
            )

        # Checking the fam file bim
        for suffix in (".bim", ".fam"):
            filename = args.bfile + suffix
            if not path.isfile(filename):
                raise ProgramError(f"{filename}: no such file")

        # Checking the intensity file
        if not path.isfile(args.intensities):
            raise ProgramError(f"{args.intensities}: no such file")

    else:
        if args.bfile is not None or args.intensities is not None:
            raise ProgramError(
                "need to specify either '--bfile' and '--intensities', or "
                "'--summarized-intensities'"
            )

        if not path.isfile(args.summarized_intensities):
            raise ProgramError(f"{args.summarized_intensities}: no such file")

    # Checking the sex problem file
    if args.sex_mismatches is not None:
        if not path.isfile(args.sex_mismatches):
            raise ProgramError(f"{args.sex_mismatches}: no such file")

    if args.intensity_nb_lines != "all":
        try:
            args.intensity_nb_lines = int(args.intensity_nb_lines)
        except ValueError:
            raise ProgramError(
                f"{args.intensity_nb_lines}: invalid number of lines",
            )


def parse_args(argv=None):
    """Parses the arguments and function."""
    parser = argparse.ArgumentParser(description=DESCRIPTION)

    parser.add_argument(
        "-v", "--version", action="version",
        version=f"pyGenClean {SCRIPT_NAME} {__version__}",
    )

    # Adding the arguments
    add_args(parser)

    if argv is None:
        return parser.parse_args()
    return parser.parse_args(argv)


def add_args(parser):
    """Adds argument to the parser."""
    # Input files
    group = parser.add_argument_group("Input files")
    group.add_argument(
        "--bfile", type=str, metavar="FILE",
        help="The Plink binary prefix containing information about markers "
             "and individuals. Must be specified if "
             "'--summarized-intensities' is not.",
    )
    group.add_argument(
        "--intensities", type=str, metavar="FILE",
        help="A file containing alleles intensities for each of the markers "
             "located on the X and Y chromosome. Must be specified if "
             "'--summarized-intensities' is not.",
    )
    group.add_argument(
        "--summarized-intensities", type=str, metavar="FILE",
        help="The name of the file containing the summarized intensities.",
    )
    group.add_argument(
        "--sex-mismatches", type=str, metavar="FILE",
        help="The file containing individuals with sex problems. This file is "
             "not read if the option 'summarized-intensities' is used.",
    )

    # The plot options
    group = parser.add_argument_group("Plot options")
    group.add_argument(
        "--format", type=str, metavar="FORMAT", default="png",
        choices=["png", "ps", "pdf"],
        help="The output file format (png, ps or pdf formats are available). "
             "[default: %(default)s]",
    )
    group.add_argument(
        "--dpi", type=int, metavar="DPI", default=300,
        help="The pixel density of the figure(s) (DPI). [%(default)d]",
    )
    group.add_argument(
        "--xlabel", type=str, metavar="STRING", default="X intensity",
        help="The label of the X axis. [default: %(default)s]",
    )
    group.add_argument(
        "--ylabel", type=str, metavar="STRING", default="Y intensity",
        help="The label of the Y axis. [default: %(default)s]",
    )

    # The intensity file options
    group = parser.add_argument_group("Intensity file options")
    group.add_argument(
        "--intensity-snp-col", type=str, metavar="COL", default="SNP Name",
        help="The name of the column containing marker names in the intensity "
             "file. [%(default)s]",
    )
    group.add_argument(
        "--intensity-sample-col", type=str, metavar="COL",
        default="Sample Name",
        help="The name of the column containing sample IDs in the intensity "
             "file. [%(default)s]",
    )
    group.add_argument(
        "--intensity-x-col", type=str, metavar="COL", default="X",
        help="The name of the column containing X values in the intensity "
             "file. [%(default)s]",
    )
    group.add_argument(
        "--intensity-y-col", type=str, metavar="COL", default="Y",
        help="The name of the column containing Y values in the intensity "
             "file. [%(default)s]",
    )
    group.add_argument(
        "--intensity-chrom-col", type=str, metavar="COL", default="Chr",
        help="The name of the column containing chromosomes in the intensity "
             "file. [%(default)s]",
    )
    group.add_argument(
        "--intensity-delimiter", type=str, metavar="SEP", default="\t",
        help="The field separator. [tabulation]",
    )
    group.add_argument(
        "--intensity-nb-lines", type=str, metavar="INT", default=10e6,
        help="The number of line to read from the intensity file at a time. "
             "Use 'all' to read the file in one go (use at your own risk)."
             "[%(default)d]",
    )

    # The OUTPUT files
    group = parser.add_argument_group("Output files")
    group.add_argument(
        "--out", type=str, metavar="FILE", default="sexcheck",
        help="The prefix of the output files. [default: %(default)s]",
    )
