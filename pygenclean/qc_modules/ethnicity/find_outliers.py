"""Finds outliers in SOURCE from CEU samples."""


import argparse
import logging
import re
from os import path
from typing import Dict, List, Optional, Set, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.cluster import KMeans
from sklearn.metrics.pairwise import euclidean_distances

from ...error import ProgramError
from ...utils import timer
from ...version import pygenclean_version as __version__
from .plot_mds import read_mds, read_populations


SCRIPT_NAME = "find-outliers"
DESCRIPTION = "Finds outliers in SOURCE from CEU samples."


logger = logging.getLogger(__name__)


@timer(logger)
def main(args: Optional[argparse.Namespace] = None,
         argv: Optional[List[str]] = None) -> None:
    """The main function of the module.

    Args:
        args (argparse.Namespace): the arguments and options.
        argv (list): the argument as a list.

    These are the steps of the modules:

    1. Prints the options.
    2. Reads the population file (:py:func:`read_population_file`).
    3. Reads the ``mds`` file (:py:func:`read_mds_file`).
    4. Computes the three reference population clusters' centers
       (:py:func:`find_ref_centers`).
    5. Computes three clusters according to the reference population clusters'
       centers, and finds the outliers of a given reference population
       (:py:func:`find_outliers`). This steps also produce three different
       plots.
    6. Writes outliers in a file (``prefix.outliers``).

    """
    # Getting and checking the options
    if args is None:
        args = parse_args(argv)
    check_args(args)

    # Reads the population file
    logger.info("Reading population file")
    populations = read_populations(args.population_file)

    # Reads the MDS file
    logger.info("Reading MDS file")
    mds = read_mds(args.mds, populations)

    # Finds the population centers
    logger.info("Finding reference population centers")
    centers, center_info = find_ref_centers(mds, args.xaxis, args.yaxis)

    # Computes three clusters using KMeans and the reference cluster centers
    logger.info("Finding outliers")
    outliers = find_outliers(mds, centers, center_info, args)
    logger.info(
        "  - There are %d outliers for the %s population",
        len(outliers), args.outliers_of,
    )

    # Printing the outlier file
    with open(args.out + ".outliers", 'w') as f:
        for sample_id in outliers:
            print(*sample_id, sep="\t", file=f)

    # Printing the outlier population file
    with open(args.out + ".population_file_outliers", "w") as f:
        for sample_id, value in populations.iterrows():
            population = value.population
            if sample_id in outliers:
                population = "OUTLIER"
            print(*sample_id, population, sep="\t", file=f)


def find_outliers(
    mds: pd.DataFrame,
    centers: np.ndarray,
    center_info: Dict[str, int],
    args: argparse.Namespace,
) -> Set[Tuple[str, str]]:
    """Finds the outliers for a given population.

    Args:
        mds (pd.DataFrame): the ``mds`` information about each samples.
        centers (np.ndarray): the centers of the three reference population
                              clusters.
        center_info (dict): the label of the three reference population
                            clusters.
        args (argparse.Namespace): the arguments and options.

    Returns:
        set: the outliers from the ``ref_pop`` population.

    Perform a ``KMeans`` classification using the three centers from the three
    reference population cluster.

    Samples are outliers of the required reference population (``ref_pop``) if:

    * the sample is part of another reference population cluster;
    * the sample is an outlier of the desired reference population
      (``ref_pop``).

    A sample is an outlier of a given cluster :math:`C_j` if the distance
    between this sample and the center of the cluster :math:`C_j` (:math:`O_j`)
    is bigger than a constant times the cluster's standard deviation
    :math:`\\sigma_j`.

    .. math::
        \\sigma_j = \\sqrt{\\frac{\\sum{d(s_i,O_j)^2}}{||C_j|| - 1}}

    where :math:`||C_j||` is the number of samples in the cluster :math:`C_j`,
    and :math:`d(s_i,O_j)` is the distance between the sample :math:`s_i` and
    the center :math:`O_j` of the cluster :math:`C_j`.

    .. math::
        d(s_i, O_j) = \\sqrt{(x_{O_j} - x_{s_i})^2 + (y_{O_j} - y_{s_i})^2}

    Using a constant equals to one ensure we remove 100% of the outliers from
    the cluster. Using a constant of 1.6 or 1.9 ensures we remove 99% and 95%
    of outliers, respectively (an error rate of 1% and 5%, respectively).

    """
    # Formatting the data
    data = mds.loc[:, [args.xaxis, args.yaxis]].to_numpy()

    # Configuring and running the KMeans
    k_means = KMeans(init=centers, n_clusters=3, n_init=1)
    k_means.fit_predict(data)

    # Creating the figure and axes
    fig_before, axe_before = plt.subplots()
    fig_after, axe_after = plt.subplots()
    fig_outliers, axe_outliers = plt.subplots()

    # Setting the title and labels for "before"
    axe_before.set_title("Before finding outliers", weight="bold")
    axe_before.set_xlabel(args.xaxis)
    axe_before.set_ylabel(args.yaxis)

    # Setting the title and labels for "after"
    axe_after.set_title(f"After finding outliers\n($> "
                        f"{args.multiplier} \\sigma$)", weight="bold")
    axe_after.set_xlabel(args.xaxis)
    axe_after.set_ylabel(args.yaxis)

    # Setting the title and labels for "outliers"
    axe_outliers.set_title("Outliers", weight="bold")
    axe_outliers.set_xlabel(args.xaxis)
    axe_outliers.set_ylabel(args.yaxis)

    # The population name
    ref_pop_name = {v: k for k, v in center_info.items()}

    # The colors
    colors = ["#CC0000", "#669900", "#0099CC"]
    outlier_colors = ["#FFCACA", "#E2F0B6", "#C5EAF8"]

    # Plotting each of the clusters with the initial center
    outliers_set = set()
    for label in range(3):
        # Subsetting the data
        subset_mds = mds.loc[k_means.labels_ == label, :]
        subset_data = data[k_means.labels_ == label]

        # Plotting the cluster
        axe_before.scatter(
            subset_mds[args.xaxis],
            subset_mds[args.yaxis],
            c=colors[label],
            s=1,
            label=ref_pop_name[label],
        )

        # Plotting the cluster center (the real one)
        axe_before.plot(
            centers[label][0],
            centers[label][1],
            "o",
            mec="#000000",
            mfc="#FFBB33",
            ms=5,
        )

        # Computing the distances
        distances = euclidean_distances(
            subset_data, centers[label].reshape(1, 2),
        )

        # Finding the outliers (that are not in the reference populations
        sigma = np.sqrt(
            np.true_divide(np.sum(distances ** 2), len(distances) - 1),
        )
        outliers = np.logical_and(
            (distances > args.multiplier * sigma).flatten(),
            subset_mds["population"] != ref_pop_name[label],
        )
        logger.info(
            "  - %d outliers for the %s cluster",
            np.count_nonzero(outliers), ref_pop_name[label],
        )

        not_ref_pop = subset_mds["population"] != ref_pop_name[label]

        if ref_pop_name[label] != args.outliers_of:
            # This is not the population we want, hence everybody is an outlier
            # (we don't include the reference population).
            outliers_set |= set(not_ref_pop.index[not_ref_pop].to_list())

            # Plotting all samples that are not part of the reference
            # populations
            axe_outliers.scatter(
                subset_mds.loc[not_ref_pop, args.xaxis],
                subset_mds.loc[not_ref_pop, args.yaxis],
                c="#555555",
                s=1,
            )

        else:
            # Plotting the outliers
            axe_outliers.scatter(
                subset_mds.loc[not_ref_pop & outliers, args.xaxis],
                subset_mds.loc[not_ref_pop & outliers, args.yaxis],
                c="#555555",
                s=1,
                label="OUTLIERS",
            )

            # Plotting the not outliers
            axe_outliers.scatter(
                subset_mds.loc[not_ref_pop & (~outliers), args.xaxis],
                subset_mds.loc[not_ref_pop & (~outliers), args.yaxis],
                c="#FFBB33",
                s=1,
                label="SOURCE",
            )
            outliers_set |= set(
                subset_mds.index[not_ref_pop & outliers].to_list()
            )

        # Plotting the cluster (without outliers)
        axe_after.scatter(
            subset_mds.loc[~outliers, args.xaxis],
            subset_mds.loc[~outliers, args.yaxis],
            c=colors[label],
            s=1,
            label=ref_pop_name[label],
        )

        # Plotting the cluster (only outliers)
        axe_after.scatter(
            subset_mds.loc[outliers, args.xaxis],
            subset_mds.loc[outliers, args.yaxis],
            c=outlier_colors[label],
            s=1,
        )

        # Plotting only the reference populations
        axe_outliers.scatter(
            subset_mds.loc[~not_ref_pop, args.xaxis],
            subset_mds.loc[~not_ref_pop, args.yaxis],
            c=colors[label],
            s=1,
            label=ref_pop_name[label],
        )

        # Plotting the cluster center (the real one)
        axe_after.plot(
            centers[label][0],
            centers[label][1],
            "o",
            mec="#000000",
            mfc="#FFBB33",
            ms=5,
        )

    # The legends
    axe_before.legend(loc="best", fancybox=True, fontsize=6)
    axe_after.legend(loc="best", fancybox=True, fontsize=6)
    axe_outliers.legend(loc="best", fancybox=True, fontsize=6)

    # Saving the figure
    fig_before.savefig(f"{args.out}.before.{args.format}", dpi=300)
    fig_after.savefig(f"{args.out}.after.{args.format}", dpi=300)
    fig_outliers.savefig(f"{args.out}.outliers.{args.format}", dpi=300)

    return outliers_set


def find_ref_centers(
    mds: pd.DataFrame,
    xaxis: str,
    yaxis: str,
) -> Tuple[np.ndarray, Dict[str, int]]:
    """Finds the center of the three reference clusters.

    Args:
        mds (pd.DataFrame): the ``mds`` information about each samples.

    Returns:
        tuple: a tuple with a :py:class:`numpy.array` containing the centers of
               the three reference population cluster as first element, and a
               :py:class:`dict` containing the label of each of the three
               reference population clusters.

    First, we extract the ``mds`` values of each of the three reference
    populations. The, we compute the center of each of those clusters by
    computing the means.

    .. math::
        \\textrm{Cluster}_\\textrm{pop} = \\left(
            \\frac{\\sum_{i=1}^n x_i}{n}, \\frac{\\sum_{i=1}^n y_i}{n}
        \\right)

    """
    # Computing the centers of each of the reference clusters
    ceu_mds = mds.loc[mds.population == "CEU", [xaxis, yaxis]]
    yri_mds = mds.loc[mds.population == "YRI", [xaxis, yaxis]]
    asn_mds = mds.loc[mds.population == "JPT-CHB", [xaxis, yaxis]]

    # Computing the centers
    centers = [[ceu_mds[xaxis].mean(), ceu_mds[yaxis].mean()],
               [yri_mds[xaxis].mean(), yri_mds[yaxis].mean()],
               [asn_mds[xaxis].mean(), asn_mds[yaxis].mean()]]

    return np.array(centers), {"CEU": 0, "YRI": 1, "JPT-CHB": 2}


def check_args(args: argparse.Namespace) -> None:
    """Checks the arguments and options.

    Args:
        args (argparse.Namespace): the arguments and options

    If there is a problem with an option, an exception is raised using the
    :py:class:`ProgramError` class, a message is printed to the
    :class:`sys.stderr` and the program exists with code 1.

    """
    # Checking the input files
    if not path.isfile(args.mds):
        raise ProgramError(f"{args.mds}: no such file")

    if not path.isfile(args.population_file):
        raise ProgramError(f"{args.population_file}: no such file")

    # Checking the chosen components
    component_re = re.compile(r"C\d+$")
    for axis in (args.xaxis, args.yaxis):
        if not component_re.match(axis):
            raise ProgramError(f"{axis}: invalid component")
    if args.xaxis == args.yaxis:
        raise ProgramError("xaxis must be different than yaxis")


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
        "--mds", type=str, metavar="FILE", required=True,
        help="The MDS file from Plink",
    )
    group.add_argument(
        "--population-file", type=str, metavar="FILE", required=True,
        help="A population file containing the following columns (without a "
             "header): FID, IID and POP. POP should be one of 'CEU', "
             "'JPT-CHB', 'YRI' and SOURCE.",
    )

    # The options
    group = parser.add_argument_group("Graphical Options")
    add_graphical_options(group)

    # The OUTPUT files
    group = parser.add_argument_group("Output File")
    group.add_argument(
        "--out", type=str, metavar="FILE", default="ethnicity",
        help="The prefix of the output files. [default: %(default)s]",
    )


def add_graphical_options(parser: argparse._ArgumentGroup,
                          prefix: str = "") -> None:
    """Adds the graphical options."""
    parser.add_argument(
        f"--{prefix}format", type=str, metavar="FORMAT", default="png",
        choices=["png", "ps", "pdf"],
        help="The output file format (png, ps, or pdf formats are available). "
             "[default: %(default)s]",
    )
    parser.add_argument(
        "--outliers-of", type=str, metavar="POP", default="CEU",
        choices=["CEU", "YRI", "JPT-CHB"],
        help="Finds the outliers of this population. [default: %(default)s]",
    )
    parser.add_argument(
        "--multiplier", type=float, metavar="FLOAT", default=1.9,
        help="To find the outliers, we look for more than x times the cluster "
             "standard deviation. [default: %(default).1f]",
    )
    parser.add_argument(
        f"--{prefix}xaxis", type=str, metavar="COMPONENT", default="C1",
        help="The component to use for the X axis. [default: %(default)s]",
    )
    parser.add_argument(
        f"--{prefix}yaxis", type=str, metavar="COMPONENT", default="C2",
        help="The component to use for the Y axis. [default: %(default)s]",
    )
