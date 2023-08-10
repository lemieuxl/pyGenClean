"""Create an MDS plot."""


import argparse
import logging
from os import path
from typing import List, Optional

import matplotlib.pyplot as plt
import pandas as pd

from ...error import ProgramError
from ...utils import timer
from ...version import pygenclean_version as __version__


SCRIPT_NAME = "plot-mds"
DESCRIPTION = "Create an MDS plot."


logger = logging.getLogger(__name__)


@timer(logger)
def main(args: Optional[argparse.Namespace] = None,
         argv: Optional[List[str]] = None) -> None:
    """The main function of the module.

    Args:
        args (argparse.Namespace): the arguments and options.
        argv (list): the argument as a list.

    These are the steps:

    1. Reads the population file (:py:func:`readPopulations`).
    2. Extract the MDS data (:py:func:`extractData`).
    3. Plots the MDS data (:py:func:`plotMDS`).

    """
    if args is None:
        args = parse_args(argv)
    check_args(args)

    # Reads the population file
    populations = read_populations(args.population_file)

    # Reading the MDS data
    mds = read_mds(args.file, populations)

    # Plot the data
    plot_mds(mds, args)


def read_populations(filename: str) -> pd.DataFrame:
    """Reads a population file.

    Args:
        filename (str): the name of the population file.

    Returns:
        dict: Dictionary of population for each of the samples.

    """
    return pd.read_csv(
        filename,
        sep="\t",
        names=["fid", "iid", "population"],
        dtype={"fid": str, "iid": str},
    ).set_index(["fid", "iid"], verify_integrity=True).sort_index()


def plot_mds(df: pd.DataFrame, args: argparse.Namespace) -> None:
    """Plot the MDS data.

    Args:
        df (pd.DataFrame): the data to plot (MDS values and population).
        args (argparse.Namespace): the options.

    """
    figure, axe = plt.subplots()

    for i, pop_name in enumerate(args.population_order):
        # The data for the population
        mds = df.loc[df.population == pop_name, :]

        # The scatter
        axe.scatter(
            mds[args.xaxis],
            mds[args.yaxis],
            c=args.population_color[i],
            s=args.population_size[i],
            marker=args.population_marker[i],
            label=pop_name,
        )

    # The legend
    axe.legend(loc="best", fancybox=True, fontsize=10)

    # The title and XY labels
    axe.set_title(args.title)
    axe.set_xlabel(args.xaxis)
    axe.set_ylabel(args.yaxis)

    if args.format == "X11":
        plt.show()

    else:
        plt.savefig(args.out + "." + args.format, dpi=300)

    plt.close(figure)


def read_mds(filename: str, populations: pd.DataFrame) -> pd.DataFrame:
    """Reads the MDS data and adds the population.

    Args:
        filename (str): the name of the MDS file.
        populations (pd.DataFrame): the population of each sample.

    Returns:
        pd.DataFrame: the MDS data with information about the population.

    """
    # Reading the MDS file
    mds = pd.read_csv(
        filename,
        delim_whitespace=True,
        dtype={"FID": str, "IID": str},
    ).set_index(["FID", "IID"], verify_integrity=True).sort_index()

    # Adding the population
    mds = mds.assign(population=populations.population)

    if mds.isnull().any().any():
        raise ProgramError(
            f"{filename}: some samples dont't have a population",
        )

    return mds


def check_args(args: argparse.Namespace) -> None:
    """Checks the arguments and options.

    Args:
        args (argparse.Namespace): the arguments and options.

    If there is a problem with an option, an exception is raised using the
    :py:class:`ProgramError` class, a message is printed to the
    :py:class:`sys.stderr` and the program exists with code 1.

    """
    # Check in input file
    if not path.isfile(args.file):
        raise ProgramError(f"{args.file}: no such file")

    # Check the population file
    if not path.isfile(args.population_file):
        raise ProgramError(f"{args.population_file}: no such file")

    # The number of population
    nb_pop = len(args.population_order)

    # Checking the colors
    if len(args.population_color) < nb_pop:
        raise ProgramError("--population-color: no enough marker colors")

    # Checking the markers
    if len(args.population_marker) < nb_pop:
        raise ProgramError("--population-marker: no enough markers")

    if len(args.population_size) < nb_pop:
        raise ProgramError("--population-size: no enough marker sizes")


def parse_args(argv: Optional[List[str]] = None) -> argparse.Namespace:
    """Parses the command line options and arguments."""
    parser = argparse.ArgumentParser(description=DESCRIPTION)

    parser.add_argument(
        "-v", "--version", action="version",
        version=f"pyGenClean {SCRIPT_NAME} {__version__}",
    )

    # Adding the arguments and options
    add_args(parser)

    return parser.parse_args(argv)


def add_args(parser: argparse.ArgumentParser) -> None:
    """Add arguments and options to the parser."""
    # The INPUT files
    group = parser.add_argument_group("Input File")
    group.add_argument(
        "--file", type=str, metavar="FILE", required=True,
        help="The MBS file.",
    )
    group.add_argument(
        "--population-file", type=str, metavar="FORMAT", required=True,
        help="A file containing population information. There must be three "
             "columns: famID, indID and population information.",
    )

    # The graphical options
    group = parser.add_argument_group("Graphical Options")
    add_graphical_options(group)

    # The population options
    group = parser.add_argument_group("Population Options")
    group.add_argument(
        "--population-order", type=str, nargs="+", metavar="POP",
        default=["CEU", "YRI", "JPT-CHB", "SOURCE"],
        help="The population order. [%(default)s]",
    )
    group.add_argument(
        "--population-color", type=str, nargs="+", metavar="COLOR",
        default=["#3700B8", "#4DAF4A", "#984EA3", "#E41A1C"],
        help="The population marker color. [%(default)s]",
    )
    group.add_argument(
        "--population-marker", type=str, nargs="+", metavar="MARKER",
        default=[".", ".", ".", "+"],
        help="The population marker. [%(default)s]",
    )
    group.add_argument(
        "--population-size", type=int, nargs="+", metavar="SIZE",
        default=[60, 60, 60, 60],
        help="The population marker size. [%(default)s]",
    )

    # The OUTPUT files
    group = parser.add_argument_group("Output File")
    group.add_argument(
        "--out", type=str, metavar="FILE", default="mds",
        help="The prefix of the output files. [default: %(default)s]",
    )


def add_graphical_options(parser: argparse._ArgumentGroup,
                          prefix: str = "") -> None:
    """Adds the graphical options."""
    parser.add_argument(
        f"--{prefix}format", type=str, metavar="FORMAT", default="png",
        choices=["png", "ps", "pdf", "X11"],
        help="The output file format (png, ps, pdf, or X11 formats are "
             "available). [default: %(default)s]",
    )
    parser.add_argument(
        f"--{prefix}title", type=str, metavar="STRING",
        default="C2 in function of C1 - MDS",
        help="The title of the MDS plot. [default: %(default)s]",
    )
    parser.add_argument(
        f"--{prefix}xaxis", type=str, metavar="STRING", default="C1",
        help="The component to use for the X axis. [default: %(default)s]",
    )
    parser.add_argument(
        f"--{prefix}yaxis", type=str, metavar="STRING", default="C2",
        help="The component to use for the Y axis. [default: %(default)s]",
    )
