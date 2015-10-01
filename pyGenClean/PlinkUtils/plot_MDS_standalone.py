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
import sys
import logging
import argparse

import numpy as np

from .. import __version__
from . import createRowFromPlinkSpacedOutput


logger = logging.getLogger("plot_MDS")


def main():
    """The main function of the module.

    These are the steps:

    1. Reads the population file (:py:func:`readPopulations`).
    2. Extracts the MDS values (:py:func:`extractData`).
    3. Plots the MDS values (:py:func:`plotMDS`).

    """
    # Getting and checking the options
    args = parseArgs()
    checkArgs(args)

    # Reads the population file
    populations = readPopulations(args.population_file, args.population_order)

    # Acquire the data
    theData, theLabels = extractData(args.file, populations,
                                     args.population_order, args.xaxis,
                                     args.yaxis)

    # Plot the data
    plotMDS(theData, args.population_order, theLabels, args.population_colors,
            args.population_alpha, args.population_sizes,
            args.population_markers, args)


def readPopulations(inputFileName, requiredPopulation):
    """Reads a population file.

    :param inputFileName: the name of the population file.
    :param requiredPopulation: the required population.

    :type inputFileName: str
    :type requiredPopulation: list

    :returns: a :py:class:`dict` containing the population of each samples.

    """
    populations = {}
    requiredPopulation = set(requiredPopulation)
    with open(inputFileName, "r") as inputFile:
        for line in inputFile:
            row = line.rstrip("\r\n").split("\t")

            # Getting the informations
            famID = row[0]
            indID = row[1]
            pop = row[2]

            # Check if we already saw this sample
            if (famID, indID) in populations:
                if pop != populations[(famID, indID)]:
                    msg = ("{} {}: sample has multiple population ({} and "
                           "{})".format(famID, indID, pop,
                                        populations[(famID, indID)]))
                    raise ProgramError(msg)

            # Save the population if we need it
            if pop in requiredPopulation:
                # We need this population
                populations[(famID, indID)] = pop

    popMissing = requiredPopulation - set(populations.values())
    if len(popMissing) != 0:
        msg = "Population that were asked for doesn't exists in " \
              "population file: %s" % str(popMissing)
        raise ProgramError(msg)

    return populations


def plotMDS(data, theOrders, theLabels, theColors, theAlphas, theSizes,
            theMarkers, options):
    """Plot the MDS data.

    :param data: the data to plot (MDS values).
    :param theOrders: the order of the populations to plot.
    :param theLabels: the names of the populations to plot.
    :param theColors: the colors of the populations to plot.
    :param theAlphas: the alpha value for the populations to plot.
    :param theSizes: the sizes of the markers for each population to plot.
    :param theMarkers: the type of marker for each population to plot.
    :param options: the options.

    :type data: list of numpy.array
    :type theOrders: list
    :type theLabels: list
    :type theColors: list
    :type theAlphas: list
    :type theSizes: list
    :type theMarkers: list
    :type options: argparse.Namespace

    """
    # Do the import
    import matplotlib as mpl
    if options.format != "X11" and mpl.get_backend() != "agg":
        mpl.use("Agg")
    import matplotlib.pyplot as plt
    if options.format != "X11":
        plt.ioff()

    fig = plt.figure()
    try:
        fig.subplots_adjust(right=options.adjust_right,
                            left=options.adjust_left,
                            bottom=options.adjust_bottom,
                            top=options.adjust_top)
    except ValueError as e:
        raise ProgramError(e)
    ax = fig.add_subplot(111)

    # Setting the axis
    ax.xaxis.set_ticks_position("bottom")
    ax.yaxis.set_ticks_position("left")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["bottom"].set_position(("outward", 9))
    ax.spines["left"].set_position(("outward", 9))

    # The plot
    plotObject = []
    labels = []
    for i, index in enumerate(theOrders):
        try:
            tmp, = ax.plot(data[0][i], data[1][i], theMarkers[i],
                           color=theColors[i], mec=theColors[i],
                           markersize=theSizes[i], alpha=theAlphas[i])
        except ValueError as e:
            msg = "Problem with markers: %(e)s" % locals()
            raise ProgramError(msg)
        plotObject.append(tmp)
        labels.append(index)

    # The legend
    prop = mpl.font_manager.FontProperties(size=options.legend_size)
    leg = ax.legend(plotObject, labels, loc=options.legend_position,
                    numpoints=1, fancybox=True, prop=prop,
                    ncol=options.legend_ncol)
    leg.get_frame().set_alpha(0.5)

    # The title and XY labels
    ax.set_title(options.title, fontsize=options.title_fontsize, weight="bold")
    ax.set_xlabel(options.xlabel, fontsize=options.label_fontsize)
    ax.set_ylabel(options.ylabel, fontsize=options.label_fontsize)

    # Changing the size of the tick labels
    for tick in ax.yaxis.get_major_ticks() + ax.xaxis.get_major_ticks():
        tick.label.set_fontsize(options.axis_fontsize)

    if options.format == "X11":
        # Show the plot
        plt.show()
    else:
        fileName = options.out + "." + options.format
        try:
            plt.savefig(fileName, dpi=300)
        except IOError:
            msg = "%(fileName)s: can't write file" % locals()
            raise ProgramError(msg)
        except ValueError as e:
            colorError = False
            for errorMsg in str(e).split("\n"):
                if errorMsg.startswith("to_rgb"):
                    colorError = True
            if colorError:
                msg = "problem with the population colors"
                raise ProgramError(msg)
            else:
                print str(e)


def extractData(fileName, populations, population_order, xaxis, yaxis):
    """Extract the C1 and C2 columns for plotting.

    :param fileName: the name of the MDS file.
    :param populations: the population of each sample in the MDS file.
    :param population_order: the required population order.
    :param xaxis: the component to print as the X axis.
    :param yaxis: the component to print as the Y axis.

    :type fileName: str
    :type populations: dict
    :type population_order: list
    :type xaxis: str
    :type yaxis: str

    :returns: the MDS data with information about the population of each
              sample. The first element of the returned tuple is a tuple. The
              last element of the returned tuple is the list of the populations
              (the order is the same as in the first element). The first
              element of the first tuple is the C1 data, and the last element
              is the C2 data.

    .. note::
        If a sample in the MDS file is not in the population file, it is skip.

    """
    # The different population labels
    nbPossibleLabels = len(population_order)

    c1 = [[] for i in xrange(nbPossibleLabels)]
    c2 = [[] for i in xrange(nbPossibleLabels)]
    with open(fileName, 'r') as inputFile:
        headerIndex = None
        for i, line in enumerate(inputFile):
            row = createRowFromPlinkSpacedOutput(line)

            if i == 0:
                # This is the header
                headerIndex = dict([(row[j], j) for j in xrange(len(row))])

                for columnName in ["FID", "IID", xaxis, yaxis]:
                    if columnName not in headerIndex:
                        msg = "%(fileName)s: no column named " \
                              "%(columnName)s" % locals()
                        raise ProgramError(msg)

            else:
                # Getting the component 1 and 2
                currC1 = row[headerIndex[xaxis]]
                currC2 = row[headerIndex[yaxis]]

                # Getting the individual informations
                famID = row[headerIndex["FID"]]
                indID = row[headerIndex["IID"]]

                curLabel = ""
                if (famID, indID) in populations:
                    curLabel = populations[(famID, indID)]
                else:
                    continue

                c1[population_order.index(curLabel)].append(currC1)
                c2[population_order.index(curLabel)].append(currC2)

    return (np.array(c1), np.array(c2)), population_order


def checkArgs(args):
    """Checks the arguments and options.

    :param args: an object containing the options of the program.
    :type args: argparse.Namespace

    :returns: ``True`` if everything was OK.

    If there is a problem with an option, an exception is raised using the
    :py:class:`ProgramError` class, a message is printed to the
    :class:`sys.stderr` and the program exists with code 1.

    """
    # Check in input file
    if not os.path.isfile(args.file):
        msg = "%s: no such file" % args.file
        raise ProgramError(msg)

    # Check the population file
    if args.population_file is None:
        msg = "population-file: no population file"
        raise ProgramError(msg)
    elif not os.path.isfile(args.population_file):
        msg = "%s: no such file" % args.population_file
        raise ProgramError(msg)

    # Split the order
    args.population_order = args.population_order.split(",")
    nbPop = len(args.population_order)

    # Format the legend position
    args.legend_position = args.legend_position.replace("-", " ")

    # Format the color
    colors = []
    for color in args.population_colors.split(","):
        colors.append("#" + color)
    args.population_colors = colors
    if len(args.population_colors) < nbPop:
        msg = "not enough colors for the required populations"
        raise ProgramError(msg)

    # Format the point size
    sizes = []
    for size in args.population_sizes.split(","):
        try:
            sizes.append(int(size))
        except ValueError:
            msg = "%(size)s: not a valid size" % locals()
            raise ProgramError(msg)
    args.population_sizes = sizes
    if len(args.population_sizes) < nbPop:
        msg = "not enough sizes for the required populations"
        raise ProgramError(msg)

    # Format the point marker
    args.population_markers = args.population_markers.split(",")

    # Format the alpha values
    alphas = []
    for alpha in args.population_alpha.split(","):
        try:
            alphas.append(float(alpha))
        except ValueError:
            msg = "%(alpha)s: not a valid alpha value" % locals()
            raise ProgramError(msg)
    args.population_alpha = alphas

    # Check the legend alpha value
    if (args.legend_alpha < 0) or (args.legend_alpha > 1):
        msg = "%s: alpha for legend must be between " \
              "0 and 1" % str(args.legend_alpha)
        raise ProgramError(msg)

    return True


def parseArgs():  # pragma: no cover
    """Parses the command line options and arguments.

    :returns: A :py:class:`argparse.Namespace` object created by the
              :py:mod:`argparse` module. It contains the values of the
              different options.

    ======================== ====== =========================================
             Options          Type                Description
    ======================== ====== =========================================
    ``--file``               string The MBS file.
    ``--population-file``    string A file containing population information.
    ``--population-order``   string The order to print the different
                                    populations.
    ``--population-colors``  string The population point color in the plot.
    ``--population-sizes``   string The population point size in the plot.
    ``--population-markers`` string The population point marker in the plot.
    ``--population-alpha``   string The population alpha value in the plot.
    ``--format``             string The output file format.
    ``--title``              string The title of the MDS plot.
    ``--xaxis``              string The component to print on the X axis.
    ``--xlabel``             string The label of the X axis.
    ``--yaxis``              string The component to print on the Y axis.
    ``--ylabel``             string The label of the Y axis.
    ``--legend-position``    string The position of the legend.
    ``--legend-size``        int    The size of the legend text.
    ``--legend-ncol``        int    The number of columns for the legend.
    ``--legend-alpha``       float  The alpha value of the legend.
    ``--title-fontsize``     int    The font size of the title.
    ``--label-fontsize``     int    The font size of the X and Y labels.
    ``--axis-fontsize``      int    The font size of the X and Y axis.
    ``--adjust-left``        float  Adjust the left margin.
    ``--adjust-right``       float  Adjust the right margin.
    ``--adjust-top``         float  Adjust the top margin.
    ``--adjust-bottom``      float  Adjust the bottom margin.
    ``--out``                string The prefix of the output files.
    ======================== ====== =========================================

    .. note::
        No option check is done here (except for the one automatically done by
        argparse). Those need to be done elsewhere (see :py:func:`checkArgs`).

    """
    args = parser.parse_args()

    return args


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
desc = "Creates a MDS plot."
long_desc = None
parser = argparse.ArgumentParser(description=desc)
parser.add_argument("-v", "--version", action="version",
                    version="pyGenClean version {}".format(__version__))

# The INPUT files
group = parser.add_argument_group("Input File")
group.add_argument("--file", type=str, metavar="FILE", required=True,
                   help="The MBS file.")
group.add_argument("--population-file", type=str, metavar="FORMAT",
                   required=True,
                   help=("A file containing population information. There "
                         "must be three columns: famID, indID and population "
                         "information."))

group = parser.add_argument_group("Population Properties")
group.add_argument("--population-order", type=str, metavar="STRING",
                   default="CEU,YRI,JPT-CHB,SOURCE,OUTLIER",
                   help=("The order to print the different populations. "
                         "[default: %(default)s]"))
group.add_argument("--population-colors", type=str, metavar="STRING",
                   default="377eb8,4daf4a,984ea3,e41a1c,ff7f00",
                   help=("The population point color in the plot [default: "
                         "%(default)s]"))
group.add_argument("--population-sizes", type=str, metavar="STRING",
                   default="12,12,12,8,3",
                   help=("The population point size in the plot. [default: "
                         "%(default)s]"))
group.add_argument("--population-markers", type=str, metavar="STRING",
                   default=".,.,.,+,D",
                   help=("The population point marker in the plot. [default: "
                         "%(default)s]"))
group.add_argument("--population-alpha", type=str, metavar="STRING",
                   default="1.0,1.0,1.0,1.0,1.0",
                   help=("The population alpha value in the plot. [default: "
                         "%(default)s]"))

# The graphical options
group = parser.add_argument_group("Graphical Properties")
group.add_argument("--format", type=str, metavar="FORMAT", default="png",
                   choices=["png", "ps", "pdf", "X11"],
                   help=("The output file format (png, ps, pdf, or X11 "
                         "formats are available). [default: %(default)s]"))
group.add_argument("--title", type=str, metavar="STRING",
                   default="C2 in function of C1 - MDS",
                   help="The title of the MDS plot. [default: %(default)s]")
group.add_argument("--xaxis", type=str, metavar="STRING", default="C1",
                   help=("The component to print on the X axis. [default: "
                         "%(default)s]"))
group.add_argument("--xlabel", type=str, metavar="STRING", default="C1",
                   help="The label of the X axis. [default: %(default)s]")
group.add_argument("--yaxis", type=str, metavar="STRING", default="C2",
                   help=("The component to print on the Y axis. [default: "
                         "%(default)s]"))
group.add_argument("--ylabel", type=str, metavar="STRING", default="C2",
                   help="The label of the Y axis. [default: %(default)s]")
group.add_argument("--legend-position", type=str, metavar="STRING",
                   default="best",
                   choices=["upper-left", "upper-right", "lower-left",
                            "lower-right", "best"],
                   help="The position of the legend. [default: %(default)s]")
group.add_argument("--legend-size", type=int, metavar="INT", default=10,
                   help="The size of the legend. [default: %(default)d]")
group.add_argument("--legend-ncol", type=int, metavar="INT", default=1,
                   help=("The number of column for the legend. [default: "
                         "%(default)d]"))
group.add_argument("--legend-alpha", type=float, metavar="FLOAT", default=1.0,
                   help=("The alpha value of the legend frame. [default: "
                         "%(default).1f]"))
group.add_argument("--title-fontsize", type=int, metavar="INT", default=15,
                   help="The font size of the title. [default: %(default)d]")
group.add_argument("--label-fontsize", type=int, metavar="INT", default=12,
                   help=("The font size of the X and Y labels. [default: "
                         "%(default)d]"))
group.add_argument("--axis-fontsize", type=int, metavar="INT", default=12,
                   help=("The font size of the X and Y axis. [Default: "
                         "%(default)d]"))
group.add_argument("--adjust-left", type=float, metavar="FLOAT", default=0.12,
                   help="Adjust the left margin. [Default: %(default).2f]")
group.add_argument("--adjust-right", type=float, metavar="FLOAT", default=0.90,
                   help="Adjust the right margin. [Default: %(default).2f]")
group.add_argument("--adjust-top", type=float, metavar="FLOAT", default=0.90,
                   help="Adjust the top margin. [Default: %(default).2f]")
group.add_argument("--adjust-bottom", type=float, metavar="FLOAT",
                   default=0.10, help="Adjust the bottom margin. "
                                      "[Default: %(default).2f]")

# The OUTPUT files
group = parser.add_argument_group("Output File")
group.add_argument("--out", type=str, metavar="FILE", default="mds",
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
