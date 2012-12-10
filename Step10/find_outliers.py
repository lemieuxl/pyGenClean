#!/usr/bin/env python2.7

import os
import sys
import argparse

import numpy as npy

from PlinkUtils import createRowFromPlinkSpacedOutput as create_row

def main(argString=None):
    # Getting and checking the options
    args = parseArgs(argString)
    checkArgs(args)

    print "   - Options used:"
    for key, value in vars(args).iteritems():
        print "      --{} {}".format(key, value)

    # Reads the population file
    print "   - Reading population file"
    populations = read_population_file(args.population_file)

    # Reads the MDS file
    print "   - Reading MDS file"
    mds = read_mds_file(args.mds, args.xaxis, args.yaxis, populations)

    # Finds the population centers
    print "   - Finding reference population centers"
    centers, center_info = find_ref_centers(mds)

    # Computes three clusters using KMeans and the reference cluster centers
    print "   - Finding outliers"
    outliers = find_outliers(mds, centers, center_info, args.outliers_of, args)
    print ("   - There are {} outliers for the {} "
           "population".format(len(outliers), args.outliers_of))
    try:
        with open(args.out + ".outliers", 'wb') as output_file:
            for sample_id in outliers:
                print >>output_file, "\t".join(sample_id)
    except IOError:
        msg = "{}: can't write file".format(args.out + ".outliers")
        raise ProgramError(msg)


def find_outliers(mds, centers, center_info, ref_pop, options):
    """Finds the outliers for a given population."""
    # Importing matplotlib for plotting purpuses
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
    data = npy.array(zip(mds["c1"], mds["c2"]))

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
        sigma = npy.sqrt(npy.true_divide(npy.sum(distances ** 2),
                                         len(distances) - 1))
        outliers = npy.logical_and((distances > options.multiplier * sigma).flatten(),
                                   subset_mds != ref_pop_name[label])
        print "      - {} outliers for the {} cluster".format(npy.sum(outliers),
                                                              ref_pop_name[label])

        # Saving the outliers
        if ref_pop_name[label] != options.outliers_of:
            # This is not the population we want, hence everybody is an outlier
            # (we don't include the reference population).
            outlier_mds = subset_mds[subset_mds["pop"] != ref_pop_name[label]]
            outliers_set |= set([(i["fid"], i["iid"]) for i in outlier_mds])

            # Plotting all samples that are not part of the reference
            # populations
            axe_outliers.plot(subset_mds["c1"][subset_mds["pop"] != ref_pop_name[label]],
                              subset_mds["c2"][subset_mds["pop"] != ref_pop_name[label]],
                              "o", mec="#555555", mfc="#555555", ms=2,
                              clip_on=False)
        else:
            # This is the population we want, hence only the real outliers are
            # outliers (we don't include the reference population)
            outlier_mds = subset_mds[npy.logical_and(subset_mds["pop"] != ref_pop_name[label],
                                                     outliers)]

            # Plotting the outliers
            plot_outliers, = axe_outliers.plot(outlier_mds["c1"],
                                               outlier_mds["c2"], "o",
                                               mec="#555555", mfc="#555555",
                                               ms=2, clip_on=False)

            # Plotting the not outliers
            plot_not_outliers, = axe_outliers.plot(subset_mds["c1"][npy.logical_and(~outliers, subset_mds["pop"] != ref_pop_name[label])],
                                                   subset_mds["c2"][npy.logical_and(~outliers, subset_mds["pop"] != ref_pop_name[label])],
                                                   "o", mec="#FFBB33",
                                                   mfc="#FFBB33", ms=2,
                                                   clip_on=False)
            outliers_set |= set([(i["fid"], i["iid"]) for i in outlier_mds])

        # Plotting the cluster (without outliers)
        p, = axe_after.plot(subset_mds[~outliers]["c1"],
                            subset_mds[~outliers]["c2"], "o", mec=colors[label],
                            mfc=colors[label], ms=2, clip_on=False)
        plots_after.append(p)

        # Plotting the cluster (only outliers)
        axe_after.plot(subset_mds[outliers]["c1"], subset_mds[outliers]["c2"],
                       "o", mec=outlier_colors[label],
                       mfc=outlier_colors[label], ms=2, clip_on=False)

        # Plotting only the reference populations
        p, = axe_outliers.plot(subset_mds["c1"][subset_mds["pop"] == ref_pop_name[label]],
                               subset_mds["c2"][subset_mds["pop"] == ref_pop_name[label]],
                               "o", mec=colors[label], mfc=colors[label], ms=2,
                               clip_on=False)
        plots_outliers.append(p)

        # Plotting the cluster center (the real one)
        axe_after.plot(centers[label][0], centers[label][1], "o", mec="#000000",
                       mfc="#FFBB33", ms=6)

    # The legends
    axe_before.legend(plots_before, ref_pop_name, "best", numpoints=1,
                      fancybox=True, fontsize=8).get_frame().set_alpha(0.5)
    axe_after.legend(plots_after, ref_pop_name, "best", numpoints=1,
                     fancybox=True, fontsize=8).get_frame().set_alpha(0.5)
    axe_outliers.legend(plots_outliers + [plot_not_outliers, plot_outliers],
                        ref_pop_name + ["SOURCE", "OUTLIERS"], "best",
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
    """Finds the center of the three reference clusters."""
    # Computing the centers of each of the reference clusters
    ceu_mds = mds[mds["pop"] == "CEU"]
    yri_mds = mds[mds["pop"] == "YRI"]
    asn_mds = mds[mds["pop"] == "JPT-CHB"]

    # Computing the centers
    centers = [[npy.mean(ceu_mds["c1"]), npy.mean(ceu_mds["c2"])],
               [npy.mean(yri_mds["c1"]), npy.mean(yri_mds["c2"])],
               [npy.mean(asn_mds["c1"]), npy.mean(asn_mds["c2"])]]

    return npy.array(centers), {"CEU": 0, "YRI": 1, "JPT-CHB": 2}


def read_mds_file(file_name, c1, c2, pops):
    """Reads a MDS file."""
    mds = []
    max_fid = 0
    max_iid = 0
    with open(file_name, 'rb') as input_file:
        # Getting and checking the header
        header_index = dict([(col_name, i) for i, col_name in
                                enumerate(create_row(input_file.readline()))])
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
    mds = npy.array(mds, dtype=[("fid", "a{}".format(max_fid)),
                                ("iid", "a{}".format(max_iid)),
                                ("c1", float), ("c2", float),
                                ("pop", "a7")])

    return mds


def read_population_file(file_name):
    """Reads the population file."""
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
                msg = ("{}: sample {}: unknow population "
                       "{}".format(file_name, " ".join(sample_id), pop))
                raise ProgramError(msg)

            # Saving the population file
            pops[tuple(row[:2])] = row[-1]

    return pops


def checkArgs(args):
    """Checks the arguments and options.

    :param args: a :py:class:`Namespace` object containing the options of the
                 program.
    :type args: :py:class:`argparse.Namespace`

    :returns: ``True`` if everything was OK.

    If there is a problem with an option, an exception is raised using the
    :py:class:`ProgramError` class, a message is printed
    to the :class:`sys.stderr` and the program exists with code 1.

    """
    # Checking the input files
    if not os.path.isfile(args.mds):
        msg = "{}: no such file".format(args.mds)
        raise ProgramError(msg)

    if not os.path.isfile(args.population_file):
        msg = "{}: no such file".format(args.population_file)
        raise ProgramError(msg)

    # Checking the choosen components
    if args.xaxis == args.yaxis:
        msg = "xaxis must be different than yaxis"
        raise ProgramError(msg)

    return True


def parseArgs(argString=None): # pragma: no cover
    """Parses the command line options and arguments.

    :returns: A :py:class:`numpy.Namespace` object created by
              the :py:mod:`argparse` module. It contains the values of the
              different options.

    ===============  ======  ===================================================
        Options       Type                     Description
    ===============  ======  ===================================================
    ===============  ======  ===================================================

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
    """Adds custom options to a parser."""
    parser.add_argument("--outliers-of", type=str, metavar="POP", default="CEU",
                        choices=["CEU", "YRI", "JPT-CHB"],
                        help=("Finds the ouliers of this population. "
                              "[default: %(default)s]"))
    parser.add_argument("--multiplier", type=float, metavar="FLOAT", default=1.9,
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
    :type msg: string

    """
    def __init__(self, msg):
        """Construction of the :py:class:`ProgramError` class.

        :param msg: the message to print to the user
        :type msg: string

        """
        self.message = str(msg)

    def __str__(self):
        return self.message

# The parser object
prog = "find_outliers"
desc = """Finds outliers in SOURCE from CEU samples."""
parser = argparse.ArgumentParser(description=desc, prog=prog)


# The INPUT files
group = parser.add_argument_group("Input File")
group.add_argument("--mds", type=str, metavar="FILE", required=True,
                   help=("The MDS file from Plink"))
group.add_argument("--population-file", type=str, metavar="FILE", required=True,
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
# The OUTPUT files
group = parser.add_argument_group("Output File")
group.add_argument("--out", type=str, metavar="FILE",
                    default="ethnic",
                    help=("The prefix of the output files. [default: "
                          "%(default)s]"))

if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print >>sys.stderr, "Cancelled by user"
        sys.exit(0)
    except ProgramError as e:
        parser.error(e.message)
