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

##     # Computing 20 clusters to find biggest radius
##     print "   - Finding biggest radius of 20 clusters"
##     radius = find_biggest_radius(mds)
##     print "      - Max={}".format(radius)

    # Computes three clusters using KMeans and the reference cluster centers
    print "   - Finding outliers for the {} population".format(args.outliers_of)
    find_outliers(mds, centers, center_info, args.outliers_of, args)


def find_outliers(mds, centers, center_info, ref_pop, options):
    """Finds the outliers for a given population."""
    # Importing matplotlib for plotting purpuses
    import matplotlib as mpl
    if options.format != "X11" and mpl.get_backend() != "agg":
        mpl.use("Agg")
    import matplotlib.pyplot as plt
    if options.format != "X11":
        plt.ioff()

    from sklearn.cluster import KMeans
    from sklearn.metrics.pairwise import euclidean_distances

    # Formatting the data
    data = npy.array(zip(mds["c1"], mds["c2"]))

    # Configuring and running the KMeans
    k_means = KMeans(init=centers, n_clusters=3, n_init=1)
    k_means.fit_predict(data)

    # Creating the figure and axes
    fig, axe = plt.subplots(1, 1)

    # Setting the axe
    axe.xaxis.set_ticks_position("bottom")
    axe.yaxis.set_ticks_position("left")
    axe.spines["top"].set_visible(False)
    axe.spines["right"].set_visible(False)
    axe.spines["bottom"].set_position(("outward", 9))
    axe.spines["left"].set_position(("outward", 9))

    # Setting the title and labels
    axe.set_title("Before finding outliers", weight="bold")
    axe.set_xlabel(options.xaxis)
    axe.set_ylabel(options.yaxis)

    # The population name
    ref_pop_name = ["CEU", "YRI", "JPT-CHB"]

    # The colors
    colors = ["#CC0000", "#669900", "#0099CC"]
    outlier_colors = ["#FFCACA", "#E2F0B6", "#C5EAF8"]

    # Plotting each of the clusters with the initial center
    cluster_outlyingness = []
    cluster_max_distances = []
    for label in xrange(3):
        # Subsetting the data
        subset_mds = mds[k_means.labels_ == label]
        subset_data = data[k_means.labels_ == label]

        # Plotting the cluster
        axe.plot(subset_mds["c1"], subset_mds["c2"], "o", mec=colors[label],
                 mfc=colors[label], ms=1)

        # Plotting the cluster center (the real one)
        axe.plot(centers[label][0], centers[label][1], "o", mec="#000000",
                 mfc="#FFBB33", ms=6)

        # Computing the distances and outlyingness
        distances = euclidean_distances(subset_data, centers[label])
        max_distance = npy.max(distances)
        cluster_outlyingness.append(npy.true_divide(distances, max_distance))
        cluster_max_distances.append(max_distance)

    # Saving the figure
    if options.format == "X11":
        plt.show()
    else:
        plt.savefig("{}.before.{}".format(options.out, options.format), dpi=300)

    # Creating the figure and axes
    fig, axe = plt.subplots(1, 1)

    # Setting the axe
    axe.xaxis.set_ticks_position("bottom")
    axe.yaxis.set_ticks_position("left")
    axe.spines["top"].set_visible(False)
    axe.spines["right"].set_visible(False)
    axe.spines["bottom"].set_position(("outward", 9))
    axe.spines["left"].set_position(("outward", 9))

    # Setting the title and labels
    axe.set_title("After finding outliers", weight="bold")
    axe.set_xlabel(options.xaxis)
    axe.set_ylabel(options.yaxis)

    # Plotting each of the clusters with the initial center
    for label in xrange(3):
        # Subsetting the data
        subset_mds = mds[k_means.labels_ == label]
        subset_data = data[k_means.labels_ == label]
        outlyingness = cluster_outlyingness[label]

        # Finding the outliers
        outliers = (outlyingness > 0.216).flatten()

        # Plotting the cluster (without outliers)
        axe.plot(subset_mds[~outliers]["c1"], subset_mds[~outliers]["c2"], "o",
                 mec=colors[label], mfc=colors[label], ms=1)

        # Plotting the cluster (only outliers)
        axe.plot(subset_mds[outliers]["c1"], subset_mds[outliers]["c2"], "o",
                 mec=outlier_colors[label], mfc=outlier_colors[label], ms=1)

        # Plotting the cluster center (the real one)
        axe.plot(centers[label][0], centers[label][1], "o", mec="#000000",
                 mfc="#FFBB33", ms=6)

        # Computing the distances and outlyingness
        distances = euclidean_distances(subset_data, centers[label])
        max_distance = npy.max(distances)
        cluster_outlyingness.append(npy.true_divide(distances, max_distance))

    # Saving the figure
    if options.format == "X11":
        plt.show()
    else:
        plt.savefig("{}.after.{}".format(options.out, options.format), dpi=300)

    # Finding the cluster number for CEU
    ceu_mds = mds["pop"] == "CEU"
    ceu_label = npy.unique(k_means.labels_[ceu_mds])

##     # Changing the labels
##     new_labels = k_means.labels_.copy()
##     for label in xrange(3):
##         outliers = npy.logical_or(cluster_distances[label] >
##                                     cluster_means[label] +
##                                     options.multiply_std_by *
##                                     cluster_stds[label],
##                                   cluster_distances[label] <
##                                     cluster_means[label] -
##                                     options.multiply_std_by *
##                                     cluster_stds[label])
##         print new_labels[k_means.labels_ == label][outliers.flatten()] = label - 100
## 
##     print npy.unique(new_labels)

def find_biggest_radius(mds):
    """Finds the biggest radius after KMeans with 20 clusters."""
    # Formatting the data
    data = npy.array(zip(mds["c1"], mds["c2"]))

    # Configuring and runnning the KMeans
    k_means = KMeans(init="k-means++", n_clusters=20, n_init=300)
    k_means.fit(data)

    # Finding the bigguest radius (skipping clusters with less than 5
    # invididuals)
    biggest_radius = 0
    for label in npy.unique(k_means.labels_):
        # Getting the subset of MDS for this cluster
        current_cluster_mds = mds[k_means.labels_ == label]

        # Skip if less than 5 individuals
        if len(current_cluster_mds) <= 5:
            continue

        # This is the center for this cluster
        current_cluster_center = k_means.cluster_centers_[label]

        # Computing the distance between each point of this cluster and the
        # center
        current_distances = euclidean_distances(data[k_means.labels_ == label],
                                                current_cluster_center)

        # The maximal distance
        current_max_distance = npy.max(current_distances)

        # Is this the maximal distance
        if current_max_distance > biggest_radius:
            biggest_radius = current_max_distance

    return biggest_radius


def find_ref_centers(mds):
    """Finds the center of the three reference clusters."""
    # Computing the centers of each of the reference clusters
    ceu_mds = mds[mds["pop"] == "CEU"]
    yri_mds = mds[mds["pop"] == "YRI"]
    asn_mds = mds[mds["pop"] == "JPT-CHB"]

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


## def plot_baf_lrr(file_names, options):
##     """Plot BAF and LRR for a list of files."""
##     # importing important stuff
##     import matplotlib as mpl
##     if options.format != "X11" and mpl.get_backend() != "agg":
##         mpl.use("Agg")
##     import matplotlib.pyplot as plt
##     if options.format != "X11":
##         plt.ioff()
## 
##     # For each of the sample/files
##     for sample, file_name in file_names.iteritems():
##         data = []
## 
##         # Reading the file
##         open_func = open
##         if file_name.endswith(".gz"):
##             open_func = gzip.open
##         with open_func(file_name, 'rb') as input_file:
##             header_index = dict([(col_name, i) for i, col_name in
##                                     enumerate(input_file.readline().rstrip("\r\n").split("\t"))])
##             for col_name in {"Chr", "Position", "B Allele Freq", "Log R Ratio"}:
##                 if col_name not in header_index:
##                     msg = "{}: no column named {}".format(file_name, col_name)
##                     raise ProgramError(msg)
## 
##             # Reading the dat
##             for line in input_file:
##                 row = line.rstrip("\r\n").split("\t")
##                 
##                 # We only need X and Y chromosomes
##                 chromosome = encode_chromosome(row[header_index["Chr"]])
##                 if chromosome not in {"X", "Y"}:
##                     continue
## 
##                 # The position
##                 position = row[header_index["Position"]]
##                 try:
##                     position = int(position)
##                 except ValueError:
##                     msg = "{}: impossible position {}".format(file_name,
##                                                               position)
##                     raise ProgramError(msg)
## 
##                 # The BAF
##                 baf = row[header_index["B Allele Freq"]]
##                 try:
##                     baf = float(baf)
##                 except ValueError:
##                     msg = "{}: impossible baf {}".format(file_name, baf)
##                     raise ProgramError(msg)
## 
##                 # The LRR
##                 lrr = row[header_index["Log R Ratio"]]
##                 try:
##                     lrr = float(lrr)
##                 except ValueError:
##                     msg = "{}: impossible lrr {}".format(file_name, lrr)
##                     raise ProgramError(msg)
## 
##                 # Saving the data
##                 data.append((chromosome, position, lrr, baf))
## 
##         # Creating the numpy array
##         data = npy.array(data, dtype=[("chr", "a1"), ("pos", int),
##                                       ("lrr", float), ("baf", float)])
## 
##         # Creating the figure and axes
##         fig, axes = plt.subplots(2, 2, figsize=(20, 8))
##         plt.subplots_adjust(left=0.05, right=0.97, wspace=0.15, hspace=0.3)
##         fig.suptitle(sample, fontsize=16, weight="bold")
## 
##         # Setting subplot properties
##         for ax in axes.flatten():
##             ax.xaxis.set_ticks_position("bottom")
##             ax.yaxis.set_ticks_position("left")
##             ax.spines["top"].set_visible(False)
##             ax.spines["right"].set_visible(False)
##             ax.spines["bottom"].set_position(("outward", 9))
##             ax.spines["left"].set_position(("outward", 9))
## 
##         # Separating the axes
##         x_lrr_ax, x_baf_ax, y_lrr_ax, y_baf_ax = axes.flatten(order='F')
## 
##         # Printing the X chromosome
##         curr_chr = data["chr"] == "X"
##         x_lrr_ax.plot(data["pos"][curr_chr]/1000000.0, data["lrr"][curr_chr],
##                       "o", ms=1, mec="#0099CC",
##                       mfc="#0099CC")[0].set_clip_on(False)
##         x_baf_ax.plot(data["pos"][curr_chr]/1000000.0, data["baf"][curr_chr],
##                       "o", ms=1, mec="#669900",
##                       mfc="#669900")[0].set_clip_on(False)
##         x_lrr_ax.axhline(y=0, color="#000000", ls="--", lw=1.2)
##         x_baf_ax.axhline(y=0.5, color="#000000", ls="--", lw=1.2)
##         x_lrr_ax.set_ylabel("LRR", weight="bold")
##         x_baf_ax.set_ylabel("BAF", weight="bold")
##         x_baf_ax.set_xlabel("Position (Mb)", weight="bold")
##         x_lrr_ax.set_title("Chromosome X", weight="bold")
## 
##         # Printing the X chromosome
##         curr_chr = data["chr"] == "Y"
##         y_lrr_ax.plot(data["pos"][curr_chr]/1000000.0, data["lrr"][curr_chr],
##                       "o", ms=1, mec="#0099CC",
##                       mfc="#0099CC")[0].set_clip_on(False)
##         y_baf_ax.plot(data["pos"][curr_chr]/1000000.0, data["baf"][curr_chr],
##                       "o", ms=1, mec="#669900",
##                       mfc="#669900")[0].set_clip_on(False)
##         y_lrr_ax.axhline(y=0, color="#000000", ls="--", lw=1.2)
##         y_baf_ax.axhline(y=0.5, color="#000000", ls="--", lw=1.2)
##         y_lrr_ax.set_ylabel("LRR", weight="bold")
##         y_baf_ax.set_ylabel("BAF", weight="bold")
##         y_baf_ax.set_xlabel("Position (Mb)", weight="bold")
##         y_lrr_ax.set_title("Chromosome Y", weight="bold")
## 
##         # Saving the figure
##         if options.format == "X11":
##             plt.show()
##         else:
##             plt.savefig("{}_{}_lrr_baf.{}".format(options.out, sample,
##                                                   options.format),
##                         dpi=300)

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
group.add_argument("--outliers-of", type=str, metavar="POP", default="CEU",
                   choices=["CEU", "YRI", "JPT-CHB"],
                   help=("Finds the ouliers of this population. "
                         "[default: %(default)s]"))
group.add_argument("--multiply-std-by", type=int, metavar="INT", default=4,
                   help=("To find the outliers, we look for more than "
                         "mean +- (x * std). [default: %(default)d]"))
# The options
group = parser.add_argument_group("Options")
group.add_argument("--format", type=str, metavar="FORMAT", default="png",
                    choices=["png", "ps", "pdf", "X11"],
                    help=("The output file format (png, ps, pdf, or X11 "
                          "formats are available). [default: %(default)s]"))
group.add_argument("--xaxis", type=str, metavar="COMPONENT", default="C1",
                   choices=["C{}".format(i) for i in xrange(1, 11)],
                   help=("The component to use for the X axis. "
                         "[default: %(default)s]"))
group.add_argument("--yaxis", type=str, metavar="COMPONENT", default="C2",
                   choices=["C{}".format(i) for i in xrange(1, 11)],
                   help=("The component to use for the Y axis. "
                         "[default: %(default)s]"))
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
