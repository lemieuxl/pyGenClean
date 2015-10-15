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

from . import createRowFromPlinkSpacedOutput


logger = logging.getLogger("plot_MDS")


desc = "Creates a MDS plot"
long_desc = None
parser = argparse.ArgumentParser(description=desc)


def main():
    """The main function of the module.

    These are the steps:

    1. Reads the population file (:py:func:`readPopulations`).
    2. Extract the MDS data (:py:func:`extractData`).
    3. Plots the MDS data (:py:func:`plotMDS`).

    """
    # Getting and checking the options
    args = parseArgs()
    checkArgs(args)

    # Reads the population file
    populations = readPopulations(args.population_file)

    # Acquire the data
    theData, theLabels = extractData(args.file, populations)
    theOrders = range(len(theLabels))
    theColors = [(0.215686275, 0.000494118, 0.721568627),
                 (0.301960784, 0.68627451, 0.290196078),
                 (0.596078431, 0.305882353, 0.639215686),
                 (0.894117647, 0.101960784, 0.109803922),
                 (0.596078431, 0.305882353, 0.639215686)]
    theSizes = [12, 12, 12, 8, 12]
    theMarkers = [".", ".", ".", "+", "D"]

    # Plot the data
    plotMDS(theData, theOrders, theLabels, theColors, theSizes, theMarkers,
            args)


def readPopulations(inputFileName):
    """Reads a population file.

    :param inputFileName: the name of the population file.

    :type inputFileName: str

    :returns: a :py:class:`dict` of population for each of the samples.

    """
    populations = {}
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

            # Save the population
            populations[(famID, indID)] = pop

    return populations


def plotMDS(data, theOrders, theLabels, theColors, theSizes, theMarkers,
            options):
    """Plot the MDS data.

    :param data: the data to plot (MDS values).
    :param theOrders: the order of the populations to plot.
    :param theLabels: the names of populations to plot.
    :param theColors: the colors of the populations to plot.
    :param theSizes: the sizes of the markers for each population to plot.
    :param theMarkers: the type of markers for each population to plot.
    :param options: the options.

    :type data: list of numpy.array
    :type theOrders: list
    :type theLabels: list
    :type theColors: list
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
    ax = fig.add_subplot(111)

    # The plot
    plotObject = []
    labels = []
    for i, index in enumerate(theOrders):
        tmp, = ax.plot(data[0][index], data[1][index], theMarkers[i],
                       color=theColors[i], mec=theColors[i],
                       markersize=theSizes[i])
        plotObject.append(tmp)
        labels.append(theLabels[index])

    # The legend
    prop = mpl.font_manager.FontProperties(size=10)
    leg = ax.legend(plotObject, labels, loc="upper left", numpoints=1,
                    fancybox=True, prop=prop)

    # The title and XY labels
    ax.set_title(options.title)
    ax.set_xlabel(options.xlabel)
    ax.set_ylabel(options.ylabel)

    if options.format == "X11":
        # Show the plot
        plt.show()
    else:
        fileName = options.out + "." + options.format
        try:
            plt.savefig(fileName)
        except IOError:
            msg = "%(fileName)s: can't write file" % locals()
            raise ProgramError(msg)


def extractData(fileName, populations):
    """Extract the C1 and C2 columns for plotting.

    :param fileName: the name of the MDS file.
    :param populations: the population of each sample in the MDS file.

    :type fileName: str
    :type fileName: dict

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
    possibleLabels = list(set(populations.values()))
    nbPossibleLabels = len(possibleLabels)

    c1 = [[] for i in xrange(nbPossibleLabels)]
    c2 = [[] for i in xrange(nbPossibleLabels)]
    with open(fileName, 'r') as inputFile:
        headerIndex = None
        for i, line in enumerate(inputFile):
            row = createRowFromPlinkSpacedOutput(line)

            if i == 0:
                # This is the header
                headerIndex = dict([(row[j], j) for j in xrange(len(row))])

                for columnName in ["FID", "IID", "C1", "C2"]:
                    if columnName not in headerIndex:
                        msg = "%(fileName)s: no column named " \
                              "%(columnName)s" % locals()
                        raise ProgramError(msg)

            else:
                # Getting the component 1 and 2
                currC1 = row[headerIndex["C1"]]
                currC2 = row[headerIndex["C2"]]

                # Getting the individual informations
                famID = row[headerIndex["FID"]]
                indID = row[headerIndex["IID"]]

                curLabel = ""
                if (famID, indID) in populations:
                    curLabel = populations[(famID, indID)]
                else:
                    continue

                c1[possibleLabels.index(curLabel)].append(currC1)
                c2[possibleLabels.index(curLabel)].append(currC2)

    return (np.array(c1), np.array(c2)), possibleLabels


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

    return True


def parseArgs():  # pragma: no cover
    """Parses the command line options and arguments.

    :returns: A :py:class:`argparse.Namespace` object created by the
              :py:mod:`argparse` module. It contains the values of the
              different options.

    ===================== ====== =========================================
            Options        Type                Description
    ===================== ====== =========================================
    ``--file``            string The MBS file.
    ``--population-file`` string A file containing population information.
    ``--format``          string The output file format.
    ``--title``           string The title of the MDS plot.
    ``--xlabel``          string The label of the X axis.
    ``--ylabel``          string The label of the Y axis.
    ``--out``             string The prefix of the output files.
    ===================== ====== =========================================

    .. note::
        No option check is done here (except for the one automatically done by
        argparse). Those need to be done elsewhere (see :py:func:`checkArgs`).

    """
    # The INPUT files
    group = parser.add_argument_group("Input File")
    group.add_argument("--file", type=str, metavar="FILE", required=True,
                       help="The MBS file.")
    parser.add_argument("--population-file", type=str, metavar="FORMAT",
                        required=True,
                        help="A file containing population information. "
                             "There must be three columns: famID, indID "
                             "and population information.")

    # The graphical options
    group = parser.add_argument_group("Graphical Options")
    addCustomOptions(group)

    # The OUTPUT files
    group = parser.add_argument_group("Output File")
    group.add_argument("--out", type=str, metavar="FILE",
                       default="mds",
                       help="The prefix of the output files. [default: "
                            "%(default)s]")

    args = parser.parse_args()

    return args


def addCustomOptions(parser):
    """Adds custom options to a parser.

    :param parser: the parser.

    :type parser: argparse.parser

    """
    parser.add_argument("--format", type=str, metavar="FORMAT", default="png",
                        choices=["png", "ps", "pdf", "X11"],
                        help="The output file format (png, ps, pdf, or X11 "
                             "formats are available). [default: %(default)s]")
    parser.add_argument("--title", type=str, metavar="STRING",
                        default="C2 in function of C1 - MDS",
                        help="The title of the MDS plot. [default: "
                             "%(default)s]")
    parser.add_argument("--xlabel", type=str, metavar="STRING", default="C1",
                        help="The label of the X axis. [default: %(default)s]")
    parser.add_argument("--ylabel", type=str, metavar="STRING", default="C2",
                        help="The label of the Y axis. [default: %(default)s]")


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
