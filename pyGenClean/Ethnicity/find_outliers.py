#!/usr/bin/env python2.7

# This file is part of pyGenClean.
#
# pyGenClean is free software: you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# pyGenClean is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
# A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with
# pyGenClean.  If not, see <http://www.gnu.org/licenses/>.


import os
import re
import sys
import logging
import argparse
from glob import glob

import numpy as np

from .. import __version__
from ..PlinkUtils import createRowFromPlinkSpacedOutput as create_row


logger = logging.getLogger("find_outliers")


def main(argString=None):
    """The main function.

    :param argString: the options.

    :type argString: list of strings

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
    args = parseArgs(argString)
    checkArgs(args)

    logger.info("Options used:")
    for key, value in vars(args).iteritems():
        logger.info("  --{} {}".format(key.replace("_", "-"), value))

    # Reads the population file
    logger.info("Reading population file")
    populations = read_population_file(args.population_file)

    # Reads the MDS file
    logger.info("Reading MDS file")
    mds = read_mds_file(args.mds, args.xaxis, args.yaxis, populations)

    # Finds the population centers
    logger.info("Finding reference population centers")
    centers, center_info = find_ref_centers(mds)

    # Computes three clusters using KMeans and the reference cluster centers
    logger.info("Finding outliers")
    outliers = find_outliers(mds, centers, center_info, args.outliers_of, args)
    logger.info("  - There are {} outliers for the {} population".format(
        len(outliers),
        args.outliers_of,
    ))

    # Printing the outlier file
    try:
        with open(args.out + ".outliers", 'w') as output_file:
            for sample_id in outliers:
                print >>output_file, "\t".join(sample_id)
    except IOError:
        msg = "{}: can't write file".format(args.out + ".outliers")
        raise ProgramError(msg)

    # Printing the outlier population file
    try:
        with open(args.out + ".population_file_outliers", "w") as output_file:
            for sample_id, population in populations.iteritems():
                if sample_id in outliers:
                    population = "OUTLIER"
                print >>output_file, "\t".join(list(sample_id) + [population])
    except IOError:
        msg = "{}: can't write file".format(
            args.out + ".population_file_outliers",
        )
        raise ProgramError(msg)

    # If there is a summary file in the working directory (for LaTeX), we want
    # to modify it, because it means that this script is run after the pipeline
    # (to modify the multiplier, for example).
    if args.overwrite_tex:
        summary_file = glob(os.path.join(os.getcwd(), "*.summary.tex"))
        if len(summary_file) == 0:
            logger.warning("No TeX summary file found")
        if len(summary_file) > 1:
            raise ProgramError("More than one TeX summary file found")
        summary_file = summary_file[0]

        # Overwriting
        overwrite_tex(summary_file, len(outliers), args)


def overwrite_tex(tex_fn, nb_outliers, script_options):
    """Overwrites the TeX summary file with new values.

    :param tex_fn: the name of the TeX summary file to overwrite.
    :param nb_outliers: the number of outliers.
    :param script_options: the script options.

    :type tex_fn: str
    :type nb_outliers: int
    :type options: argparse.Namespace

    """
    # Getting the content of the TeX file
    content = None
    with open(tex_fn, "r") as i_file:
        content = i_file.read()

    # Is there a figure
    has_figure = re.search("includegraphics", content) is not None

    # Changing the first sentence
    content, n = re.subn(
        r"(Using\s[0-9,]+\smarkers?\sand\sa\smultiplier\sof\s)[0-9.]+"
        r"(,\sthere\swas\sa\stotal\sof\s)[0-9,]+(\soutliers?\sof\sthe\s)"
        r"\w+(\spopulation.)",
        r"\g<1>{}\g<2>{:,d}\g<3>{}\g<4>".format(
            script_options.multiplier,
            nb_outliers,
            script_options.outliers_of,
        ),
        content
    )
    if n != 1:
        raise ProgramError("{}: invalid TeX summary file".format(tex_fn))

    # Do we need to change the principal components in the text?
    c1_c2 = {"C1", "C2"}
    if has_figure and ({script_options.xaxis, script_options.yaxis} != c1_c2):
        content, n = re.subn(
            r"(shows\s)the\sfirst\stwo\sprincipal\scomponents(\sof\sthe\sMDS"
            r"\sanalysis)",
            r"\g<1>components {} versus {}\g<2>".format(
                script_options.yaxis.replace("C", ""),
                script_options.xaxis.replace("C", ""),
            ),
            content,
        )
        if n != 1:
            raise ProgramError("{}: invalid TeX summary file".format(tex_fn))

        content, n = re.subn(
            r"(MDS\splots\sshowing\s)the\sfirst\stwo\sprincipal\scomponents"
            r"(\sof\sthe\ssource\sdataset\swith\sthe\sreference\spanels.)",
            r"\g<1>components {} versus {}\g<2>".format(
                script_options.yaxis.replace("C", ""),
                script_options.xaxis.replace("C", ""),
            ),
            content,
        )
        if n != 1:
            raise ProgramError("{}: invalid TeX summary file".format(tex_fn))

    if has_figure:
        # Changing the population in the figure description
        content, n = re.subn(
            r"(where\soutliers\sof\sthe\s)\w+(\spopulation\sare\sshown\s"
            r"in\sgrey.)",
            r"\g<1>{}\g<2>".format(script_options.outliers_of),
            content,
        )
        if n != 1:
            raise ProgramError("{}: invalid TeX summary file".format(tex_fn))

        content, n = re.subn(
            r"(The\soutliers\sof\sthe\s)\w+(\spopulation\sare\sshown\sin\s"
            r"grey,\swhile\ssamples\sof\sthe\ssource\sdataset\sthat\s"
            r"resemble\sthe\s)\w+(\spopulation\sare\sshown\sin\sorange.\sA\s"
            r"multiplier\sof\s)[0-9.]+(\swas\sused\sto\sfind\sthe\s)[0-9,]+(\s"
            r"outliers?.)",
            r"\g<1>{pop}\g<2>{pop}\g<3>{mult}\g<4>{nb:,d}\g<5>".format(
                pop=script_options.outliers_of,
                mult=script_options.multiplier,
                nb=nb_outliers,
            ),
            content,
        )
        if n != 1:
            raise ProgramError("{}: invalid TeX summary file".format(tex_fn))

        # Changing the figure (path)
        content, n = re.subn(
            r"(\s\\includegraphics\[.+\]\{)\{ethnicity.outliers\}.png(\}\s)",
            r"\g<1>{}\g<2>".format(
                "{" + script_options.out + ".outliers}." +
                script_options.format
            ),
            content,
        )
        if n != 1:
            raise ProgramError("{}: invalid TeX summary file".format(tex_fn))

    # Saving the new content
    with open(tex_fn, "w") as o_file:
        o_file.write(content)


def find_outliers(mds, centers, center_info, ref_pop, options):
    """Finds the outliers for a given population.

    :param mds: the ``mds`` information about each samples.
    :param centers: the centers of the three reference population clusters.
    :param center_info: the label of the three reference population clusters.
    :param ref_pop: the reference population for which we need the outliers
                    from.
    :param options: the options

    :type mds: numpy.recarray
    :type centers: numpy.array
    :type center_info: dict
    :type ref_pop: str
    :type options: argparse.Namespace

    :returns: a :py:class:`set` of outliers from the ``ref_pop`` population.

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
    # Importing matplotlib for plotting purposes
    import matplotlib as mpl
    if options.format != "X11" and mpl.get_backend() != "agg":
        mpl.use("Agg")
    import matplotlib.pyplot as plt
    if options.format != "X11":
        plt.ioff()
    import matplotlib.mlab as mlab

    from sklearn.cluster import KMeans
    from sklearn.metrics.pairwise import euclidean_distances

    # Formatting the data
    data = np.array(zip(mds["c1"], mds["c2"]))

    # Configuring and running the KMeans
    k_means = KMeans(init=centers, n_clusters=3, n_init=1)
    k_means.fit_predict(data)

    # Creating the figure and axes
    fig_before, axe_before = plt.subplots(1, 1)
    fig_after, axe_after = plt.subplots(1, 1)
    fig_outliers, axe_outliers = plt.subplots(1, 1)

    # Setting the axe
    axe_before.xaxis.set_ticks_position("bottom")
    axe_before.yaxis.set_ticks_position("left")
    axe_before.spines["top"].set_visible(False)
    axe_before.spines["right"].set_visible(False)
    axe_before.spines["bottom"].set_position(("outward", 9))
    axe_before.spines["left"].set_position(("outward", 9))
    axe_after.xaxis.set_ticks_position("bottom")
    axe_after.yaxis.set_ticks_position("left")
    axe_after.spines["top"].set_visible(False)
    axe_after.spines["right"].set_visible(False)
    axe_after.spines["bottom"].set_position(("outward", 9))
    axe_after.spines["left"].set_position(("outward", 9))
    axe_outliers.xaxis.set_ticks_position("bottom")
    axe_outliers.yaxis.set_ticks_position("left")
    axe_outliers.spines["top"].set_visible(False)
    axe_outliers.spines["right"].set_visible(False)
    axe_outliers.spines["bottom"].set_position(("outward", 9))
    axe_outliers.spines["left"].set_position(("outward", 9))

    # Setting the title and labels
    axe_before.set_title("Before finding outliers", weight="bold")
    axe_before.set_xlabel(options.xaxis)
    axe_before.set_ylabel(options.yaxis)
    axe_after.set_title("After finding outliers\n($> {} "
                        "\\sigma$)".format(options.multiplier), weight="bold")
    axe_after.set_xlabel(options.xaxis)
    axe_after.set_ylabel(options.yaxis)
    axe_outliers.set_title("Outliers", weight="bold")
    axe_outliers.set_xlabel(options.xaxis)
    axe_outliers.set_ylabel(options.yaxis)

    # The population name
    ref_pop_name = ["CEU", "YRI", "JPT-CHB"]

    # The colors
    colors = ["#CC0000", "#669900", "#0099CC"]
    outlier_colors = ["#FFCACA", "#E2F0B6", "#C5EAF8"]

    # Plotting each of the clusters with the initial center
    plots_before = []
    plots_after = []
    plots_outliers = []
    plot_outliers = None
    plot_not_outliers = None
    outliers_set = set()
    for label in xrange(3):
        # Subsetting the data
        subset_mds = mds[k_means.labels_ == label]
        subset_data = data[k_means.labels_ == label]

        # Plotting the cluster
        p, = axe_before.plot(subset_mds["c1"], subset_mds["c2"], "o",
                             mec=colors[label], mfc=colors[label], ms=2,
                             clip_on=False)
        plots_before.append(p)

        # Plotting the cluster center (the real one)
        axe_before.plot(centers[label][0], centers[label][1], "o",
                        mec="#000000", mfc="#FFBB33", ms=6, clip_on=False)

        # Computing the distances
        distances = euclidean_distances(subset_data, centers[label])

        # Finding the outliers (that are not in the reference populations
        sigma = np.sqrt(np.true_divide(np.sum(distances ** 2),
                                       len(distances) - 1))
        outliers = np.logical_and(
            (distances > options.multiplier * sigma).flatten(),
            subset_mds["pop"] != ref_pop_name[label],
        )
        logger.info("  - {} outliers for the {} cluster".format(
            np.sum(outliers),
            ref_pop_name[label],
        ))

        # Saving the outliers
        if ref_pop_name[label] != options.outliers_of:
            # This is not the population we want, hence everybody is an outlier
            # (we don't include the reference population).
            outlier_mds = subset_mds[subset_mds["pop"] != ref_pop_name[label]]
            outliers_set |= set([(i["fid"], i["iid"]) for i in outlier_mds])

            # Plotting all samples that are not part of the reference
            # populations
            axe_outliers.plot(
                subset_mds["c1"][subset_mds["pop"] != ref_pop_name[label]],
                subset_mds["c2"][subset_mds["pop"] != ref_pop_name[label]],
                "o",
                mec="#555555",
                mfc="#555555",
                ms=2,
                clip_on=False,
            )
        else:
            # This is the population we want, hence only the real outliers are
            # outliers (we don't include the reference population)
            outlier_mds = subset_mds[
                np.logical_and(subset_mds["pop"] != ref_pop_name[label],
                               outliers)
            ]

            # Plotting the outliers
            plot_outliers, = axe_outliers.plot(outlier_mds["c1"],
                                               outlier_mds["c2"], "o",
                                               mec="#555555", mfc="#555555",
                                               ms=2, clip_on=False)

            # Plotting the not outliers
            plot_not_outliers, = axe_outliers.plot(
                subset_mds["c1"][np.logical_and(
                        ~outliers,
                        subset_mds["pop"] != ref_pop_name[label]
                )],
                subset_mds["c2"][np.logical_and(
                    ~outliers,
                    subset_mds["pop"] != ref_pop_name[label]
                )],
                "o",
                mec="#FFBB33",
                mfc="#FFBB33",
                ms=2,
                clip_on=False,
            )
            outliers_set |= set([(i["fid"], i["iid"]) for i in outlier_mds])

        # Plotting the cluster (without outliers)
        p, = axe_after.plot(
            subset_mds[~outliers]["c1"],
            subset_mds[~outliers]["c2"],
            "o",
            mec=colors[label],
            mfc=colors[label],
            ms=2,
            clip_on=False,
        )
        plots_after.append(p)

        # Plotting the cluster (only outliers)
        axe_after.plot(subset_mds[outliers]["c1"], subset_mds[outliers]["c2"],
                       "o", mec=outlier_colors[label],
                       mfc=outlier_colors[label], ms=2, clip_on=False)

        # Plotting only the reference populations
        p, = axe_outliers.plot(
            subset_mds["c1"][subset_mds["pop"] == ref_pop_name[label]],
            subset_mds["c2"][subset_mds["pop"] == ref_pop_name[label]],
            "o",
            mec=colors[label],
            mfc=colors[label],
            ms=2,
            clip_on=False,
        )
        plots_outliers.append(p)

        # Plotting the cluster center (the real one)
        axe_after.plot(centers[label][0], centers[label][1], "o",
                       mec="#000000", mfc="#FFBB33", ms=6)

    # The legends
    axe_before.legend(plots_before, ref_pop_name, loc="best", numpoints=1,
                      fancybox=True, fontsize=8).get_frame().set_alpha(0.5)
    axe_after.legend(plots_after, ref_pop_name, loc="best", numpoints=1,
                     fancybox=True, fontsize=8).get_frame().set_alpha(0.5)
    axe_outliers.legend(plots_outliers + [plot_not_outliers, plot_outliers],
                        ref_pop_name + ["SOURCE", "OUTLIERS"], loc="best",
                        numpoints=1, fancybox=True,
                        fontsize=8).get_frame().set_alpha(0.5)

    # Saving the figure
    fig_before.savefig("{}.before.{}".format(options.out, options.format),
                       dpi=300)
    fig_after.savefig("{}.after.{}".format(options.out, options.format),
                      dpi=300)
    fig_outliers.savefig("{}.outliers.{}".format(options.out, options.format),
                         dpi=300)

    return outliers_set


def find_ref_centers(mds):
    """Finds the center of the three reference clusters.

    :param mds: the ``mds`` information about each samples.

    :type mds: numpy.recarray

    :returns: a tuple with a :py:class:`numpy.array` containing the centers of
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
    ceu_mds = mds[mds["pop"] == "CEU"]
    yri_mds = mds[mds["pop"] == "YRI"]
    asn_mds = mds[mds["pop"] == "JPT-CHB"]

    # Computing the centers
    centers = [[np.mean(ceu_mds["c1"]), np.mean(ceu_mds["c2"])],
               [np.mean(yri_mds["c1"]), np.mean(yri_mds["c2"])],
               [np.mean(asn_mds["c1"]), np.mean(asn_mds["c2"])]]

    return np.array(centers), {"CEU": 0, "YRI": 1, "JPT-CHB": 2}


def read_mds_file(file_name, c1, c2, pops):
    """Reads a MDS file.

    :param file_name: the name of the ``mds`` file.
    :param c1: the first component to read (x axis).
    :param c2: the second component to read (y axis).
    :param pops: the population of each sample.

    :type file_name: str
    :type c1: str
    :type c2: str
    :type pops: dict

    :returns: a :py:class:`numpy.recarray` (one sample per line) with the
              information about the family ID, the individual ID, the first
              component to extract, the second component to extract and the
              population.

    The ``mds`` file is the result of Plink (as produced by the
    :py:mod:`pyGenClean.Ethnicity.check_ethnicity` module).

    """
    mds = []
    max_fid = 0
    max_iid = 0
    with open(file_name, 'rb') as input_file:
        # Getting and checking the header
        header_index = dict([
            (col_name, i) for i, col_name in
            enumerate(create_row(input_file.readline()))
        ])
        for col_name in {"FID", "IID", c1, c2}:
            if col_name not in header_index:
                msg = "{}: no column named {}".format(file_name, col_name)
                raise ProgramError(msg)

        for row in map(create_row, input_file):
            # Getting the sample ID
            sample_id = (row[header_index["FID"]], row[header_index["IID"]])

            # Checking we have a population for this sample
            if sample_id not in pops:
                msg = "{} {}: not in population file".format(sample_id[0],
                                                             sample_id[1])
                raise ProgramError(msg)

            # Saving the data
            mds.append((sample_id[0], sample_id[1], row[header_index[c1]],
                        row[header_index[c2]], pops[sample_id]))

            # Cheking to find the max FID
            if len(sample_id[0]) > max_fid:
                max_fid = len(sample_id[0])
            if len(sample_id[1]) > max_iid:
                max_iid = len(sample_id[1])

    # Creating the numpy array
    mds = np.array(mds, dtype=[("fid", "a{}".format(max_fid)),
                               ("iid", "a{}".format(max_iid)),
                               ("c1", float), ("c2", float),
                               ("pop", "a7")])

    return mds


def read_population_file(file_name):
    """Reads the population file.

    :param file_name: the name of the population file.

    :type file_name: str

    :returns: a :py:class:`dict` containing the population for each of the
              samples.

    The population file should contain three columns:

    1. The family ID.
    2. The individual ID.
    3. The population of the file (one of ``CEU``, ``YRI``, ``JPT-CHB`` or
       ``SOURCE``).

    The outliers are from the ``SOURCE`` population, when compared to one of
    the three reference population (``CEU``, ``YRI`` or ``JPT-CHB``).

    """
    pops = {}
    required_pops = {"CEU", "YRI", "JPT-CHB", "SOURCE"}
    with open(file_name, 'rb') as input_file:
        for line in input_file:
            row = line.rstrip("\r\n").split("\t")

            # The data
            sample_id = tuple(row[:2])
            pop = row[-1]

            # Checking the pop
            if pop not in required_pops:
                msg = ("{}: sample {}: unknown population "
                       "{}".format(file_name, " ".join(sample_id), pop))
                raise ProgramError(msg)

            # Saving the population file
            pops[tuple(row[:2])] = row[-1]

    return pops


def checkArgs(args):
    """Checks the arguments and options.

    :param args: a :py:class:`argparse.Namespace` object containing the options
                 of the program.

    :type args: argparse.Namespace

    :returns: ``True`` if everything was OK.

    If there is a problem with an option, an exception is raised using the
    :py:class:`ProgramError` class, a message is printed to the
    :class:`sys.stderr` and the program exists with code 1.

    """
    # Checking the input files
    if not os.path.isfile(args.mds):
        msg = "{}: no such file".format(args.mds)
        raise ProgramError(msg)

    if not os.path.isfile(args.population_file):
        msg = "{}: no such file".format(args.population_file)
        raise ProgramError(msg)

    # Checking the chosen components
    if args.xaxis == args.yaxis:
        msg = "xaxis must be different than yaxis"
        raise ProgramError(msg)

    return True


def parseArgs(argString=None):  # pragma: no cover
    """Parses the command line options and arguments.

    :param argString: the options.

    :type argString: list

    :returns: A :py:class:`argparse.Namespace` object created by the
              :py:mod:`argparse` module. It contains the values of the
              different options.

    ===================== ====== ==============================================
            Options        Type                   Description
    ===================== ====== ==============================================
    ``--mds``             string The MDS file from Plink.
    ``--population-file`` string A population file from
                                 :py:mod:`pyGenClean.Ethnicity.check_ethnicity`
                                 module.
    ``--format``          string The output file format (png, ps, or pdf.
    ``--outliers-of``     string Finds the outliers of this population.
    ``--multiplier``      float  To find the outliers, we look for more than
                                 :math:`x` times the cluster standard
                                 deviation.
    ``--xaxis``           string The component to use for the X axis.
    ``--yaxis``           string The component to use for the Y axis.
    ``--format``          string The output file format (png, ps, or pdf
                                 formats are available).
    ``--out``             string The prefix of the output files.
    ===================== ====== ==============================================

    .. note::
        No option check is done here (except for the one automatically done by
        argparse). Those need to be done elsewhere (see :py:func:`checkArgs`).

    """
    args = None
    if argString is None:
        args = parser.parse_args()
    else:
        args = parser.parse_args(argString)

    return args


def add_custom_options(parser):
    """Adds custom options to a parser.

    :param parser: the parser to which to add options.

    :type parser: argparse.ArgumentParser

    """
    parser.add_argument("--outliers-of", type=str, metavar="POP",
                        default="CEU", choices=["CEU", "YRI", "JPT-CHB"],
                        help=("Finds the outliers of this population. "
                              "[default: %(default)s]"))
    parser.add_argument("--multiplier", type=float, metavar="FLOAT",
                        default=1.9,
                        help=("To find the outliers, we look for more than "
                              "x times the cluster standard deviation. "
                              "[default: %(default).1f]"))
    parser.add_argument("--xaxis", type=str, metavar="COMPONENT", default="C1",
                        choices=["C{}".format(i) for i in xrange(1, 11)],
                        help=("The component to use for the X axis. "
                              "[default: %(default)s]"))
    parser.add_argument("--yaxis", type=str, metavar="COMPONENT", default="C2",
                        choices=["C{}".format(i) for i in xrange(1, 11)],
                        help=("The component to use for the Y axis. "
                              "[default: %(default)s]"))


class ProgramError(Exception):
    """An :py:class:`Exception` raised in case of a problem.

    :param msg: the message to print to the user before exiting.

    :type msg: str

    """
    def __init__(self, msg):
        """Construction of the :py:class:`ProgramError` class.

        :param msg: the message to print to the user

        :type msg: str

        """
        self.message = str(msg)

    def __str__(self):
        return self.message

# The parser object
desc = "Finds outliers in SOURCE from CEU samples."
long_desc = None
parser = argparse.ArgumentParser(description=desc)
parser.add_argument("-v", "--version", action="version",
                    version="pyGenClean version {}".format(__version__))


# The INPUT files
group = parser.add_argument_group("Input File")
group.add_argument("--mds", type=str, metavar="FILE", required=True,
                   help=("The MDS file from Plink"))
group.add_argument("--population-file", type=str, metavar="FILE",
                   required=True,
                   help=("A population file containing the following columns "
                         "(without a header): FID, IID and POP. POP should be "
                         "one of 'CEU', 'JPT-CHB', 'YRI' and SOURCE."))
# The options
group = parser.add_argument_group("Options")
add_custom_options(group)
group.add_argument("--format", type=str, metavar="FORMAT", default="png",
                   choices=["png", "ps", "pdf"],
                   help=("The output file format (png, ps, or pdf "
                         "formats are available). [default: %(default)s]"))
group.add_argument("--overwrite-tex", action="store_true",
                   help=("Using this option will overwrite any file that "
                         "match '*.summary.tex' (if any). This file is "
                         "automatically generated by the pyGenClean main "
                         "pipeline."))
# The OUTPUT files
group = parser.add_argument_group("Output File")
group.add_argument("--out", type=str, metavar="FILE",
                   default="ethnicity",
                   help=("The prefix of the output files. [default: "
                         "%(default)s]"))


def safe_main():
    """A safe version of the main function (that catches ProgramError)."""
    try:
        main()
    except KeyboardInterrupt:
        logger.info("Cancelled by user")
        sys.exit(0)
    except ProgramError as e:
        logger.error(e.message)
        parser.error(e.message)


if __name__ == "__main__":
    safe_main()
