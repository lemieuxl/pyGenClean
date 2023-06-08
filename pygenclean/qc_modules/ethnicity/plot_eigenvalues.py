"""Plot eigenvalues (scree plot)."""


import argparse
import logging
import re
from os import path
from typing import List, Optional

import matplotlib.pyplot as plt
import numpy as np

from ...error import ProgramError
from ...utils import timer
from ...version import pygenclean_version as __version__


SCRIPT_NAME = "plot-eigenvalues"
DESCRIPTION = "Plot eigenvalues (scree plot)."


logger = logging.getLogger(__name__)


@timer(logger)
def main(args: Optional[argparse.Namespace] = None,
         argv: Optional[List[str]] = None) -> None:
    """The main function.

    Args:
        args (argparse.Namespace): the arguments and options.
        argv (list): the argument as a list.

    The purpose of this module is to plot Eigenvectors provided by the
    Eigensoft software.

    Here are the steps of this module:

    1. Reads the Eigenvector (:py:func:`read_eigenvalues`).
    2. Plots the Scree Plot (:py:func:`create_scree_plot`).

    """
    if args is None:
        args = parse_args(argv)
    check_args(args)

    # Reads the eigenvalues
    eigenvalues = read_eigenvalues(args.evec)

    # Creates the plot
    create_scree_plot(eigenvalues, args.out, args)


def create_scree_plot(data: np.ndarray, filename: str,
                      args: argparse.Namespace) -> None:
    """Creates the scree plot.

    Args:
        data (np.ndarray): the eigenvalues.
        filename (str): the name of the output files.
        args (argparse.Namespace): the options.

    """
    # Computing the cumulative sum
    cumul_data = np.cumsum(data)

    # Creating the figure and axes
    figure, axes = plt.subplots(2, 1, figsize=(8, 16))

    # The title of the figure
    figure.suptitle(args.title, fontsize=16, weight="bold")

    figure.subplots_adjust(hspace=0.27, top=0.93, bottom=0.06)

    # First, plotting the eigenvalues
    axes[0].set_title("Scree Plot", weight="bold")
    axes[0].set_xlabel("Component")
    axes[0].set_ylabel("Eigenvalue")
    axes[0].plot(np.arange(len(data)) + 1, data, marker="o", c="#0099CC",
                 mec="#0099CC", ls="-", lw=2, clip_on=False)

    # Then, plotting the annotation
    for i, value in enumerate(data):
        axes[0].annotate(
            np.round(value, 3),
            xy=(i + 1, value),
            xytext=(1, 10),
            textcoords="offset points",
            ha="left",
            va="bottom",
            bbox=dict(boxstyle="round,pad=0.5", fc="#FFFFFF"),
        )

    # Next plot the cumulative values
    axes[1].set_title(
        f"Cumulative explained variance (max={np.sum(data):.3f})",
        weight="bold",
    )
    axes[1].set_xlabel("Component")
    axes[1].set_ylabel("Cumulative explained variance")
    axes[1].axhline(np.sum(data) * 0.8, ls="--", lw="2", c="#999999")
    axes[1].plot(
        np.arange(len(data)) + 1,
        cumul_data,
        marker="o",
        c="#CC0000",
        mec="#CC0000",
        mfc="#CC0000",
        ls="-",
        lw=2,
    )

    # Then, plotting the annotation
    for i, value in enumerate(cumul_data):
        axes[1].annotate(
            np.round(value, 3),
            xy=(i + 1, value),
            xytext=(1, -10),
            textcoords="offset points",
            ha="left",
            va="top",
            bbox=dict(boxstyle="round,pad=0.5", fc="#FFFFFF"),
        )

    # Saving the file
    plt.savefig(filename, dpi=300)
    plt.close(figure)


def read_eigenvalues(filename: str) -> np.ndarray:
    """Reads the eigenvalues from EIGENSOFT results.

    Args:
        filename (str): the name of the input file.

    Returns:
        pd.DataFrame: the eigenvalues.

    """
    re_splitter = re.compile(r"\s+")

    # The data is the first line of the result file (should begin with #
    # "#eigvals"
    data = None
    with open(filename, "r") as f:
        data = re_splitter.split(re.sub(r"(^\s+)|(\s+$)", "", f.readline()))
        if not data[0].startswith("#eigvals"):
            raise ProgramError(f"{filename}: not a evec file")

    return np.array(data[1:], dtype=float)


def check_args(args:  argparse.Namespace) -> None:
    """Checks the arguments and options.

    Args:
        args (argparse.Namespace): the arguments and options.

    If there is a problem with an option, an exception is raised using the
    :py:class:`ProgramError` class, a message is printed to the
    :class:`sys.stderr` and the program exits with error code 1.

    """
    # Checking that the input file exists
    if not path.isfile(args.evec):
        raise ProgramError(f"{args.evec}: no such file")


def parse_args(argv: Optional[List[str]] = None) -> argparse.Namespace:
    """Parses the command line options and arguments."""
    parser = argparse.ArgumentParser(description=DESCRIPTION)

    parser.add_argument(
        "-v", "--version", action="version",
        version=f"pyGenClean {SCRIPT_NAME} {__version__}",
    )

    add_args(parser)

    return parser.parse_args(argv)


def add_args(parser: argparse.ArgumentParser) -> None:
    """Add arguments and options to the parser."""
    # The input files
    group = parser.add_argument_group("Input File")
    group.add_argument(
        "--evec", type=str, metavar="FILE", required=True,
        help="The EVEC file from EIGENSOFT",
    )

    # The plot options
    group = parser.add_argument_group("Plot Options")
    add_graphical_options(group)

    # The output
    group = parser.add_argument_group("Output Options")
    group.add_argument(
        "--out", metavar="FILE", default="scree_plot.png", type=str,
        help="The name of the output file [%(default)s]",
    )


def add_graphical_options(parser: argparse._ArgumentGroup,
                          prefix: str = "") -> None:
    """Adds the graphical options."""
    parser.add_argument(
        f"--{prefix}title", type=str, metavar="TITLE",
        default="EIGENSOFT results",
        help="The main title of the scree plot [%(default)s]",
    )
