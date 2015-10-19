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

from .. import __version__


desc = "Compares BIM file."
long_desc = None
parser = argparse.ArgumentParser(description=desc)
parser.add_argument("-v", "--version", action="version",
                    version="pyGenClean version {}".format(__version__))


logger = logging.getLogger("compare_bim")


def main():
    """ The main function of the module.

    The purpose of this module is to find markers that were removed by Plink.
    When Plinks exclude some markers from binary files, there are no easy way
    to find the list of removed markers, except by comparing the two BIM files
    (before and after modification).

    Here are the steps of this module:

    1. Reads the BIM file before the modification (:py:func:`readBIM`).
    2. Reads the BIM file after the modification (:py:func:`readBIM`).
    3. Compares the list of markers before and after modification, and write
       the removed markers into a file (:py:func:`compareSNPs`).

    .. note::
        This module only finds marker that were removed (since adding markers
        to a BIM file usually includes a companion file to tell Plink which
        marker to add.

    """
    # Getting and checking the options
    args = parseArgs()
    checkArgs(args)

    # Reading the bim files
    beforeBIM = readBIM(args.before)
    afterBIM = readBIM(args.after)

    # Find differences...
    compareSNPs(beforeBIM, afterBIM, args.out)


def compareSNPs(before, after, outFileName):
    """Compares two set of SNPs.

    :param before: the names of the markers in the ``before`` file.
    :param after: the names of the markers in the ``after`` file.
    :param outFileName: the name of the output file.

    :type before: set
    :type after: set
    :type outFileName: str

    Finds the difference between two sets of markers, and write them in the
    ``outFileName`` file.

    .. note::
        A :py:class:`ProgramError` is raised if:

        1. There are more markers in the ``after`` set than in the  ``before``
           set.
        2. Some markers that are in the ``after`` set are not in the ``before``
           set.

    """
    # First, check that "before" is larger than "after"
    if len(after) > len(before):
        msg = "there are more SNPs after than before"
        raise ProgramError(msg)

    # Checks that all the SNPs "after" are in "before"
    if not (after <= before):
        msg = "some after SNPs are not in before"
        raise ProgramError(msg)

    # Printing the SNPs
    try:
        with open(outFileName, "w") as outputFile:
            differences = before - after
            if len(differences) > 0:
                print >>outputFile, "\n".join(differences)
    except IOError:
        msg = "%(outFileName)s: can't write to file" % locals()
        raise ProgramError(msg)


def readBIM(fileName):
    """Reads a BIM file.

    :param fileName: the name of the BIM file to read.

    :type fileName: str

    :returns: the set of markers in the BIM file.

    Reads a Plink BIM file and extract the name of the markers. There is one
    marker per line, and the name of the marker is in the second column. There
    is no header in the BIM file.

    """
    # Reading the first BIM file
    snps = set()
    with open(fileName, "r") as inputFile:
        for line in inputFile:
            row = line.rstrip("\r\n").split("\t")
            snpName = row[1]
            snps.add(snpName)

    return snps


def checkArgs(args):
    """Checks the arguments and options.

    :param args: an object containing the options of the program.

    :type args: argparse.Namespace

    :returns: ``True`` if everything was OK.

    If there is a problem with an option, an exception is raised using the
    :py:class:`ProgramError` class, a message is printed to the
    :class:`sys.stderr` and the program exists with code 1.

    """
    # Check the "before" file
    if not args.before.endswith(".bim"):
        msg = "%s: not a BIM file (extension must be .bim)" % args.before
        raise ProgramError(msg)
    elif not os.path.isfile(args.before):
        msg = "%s: no such file" % args.before
        raise ProgramError(msg)

    # Check the "after" file
    if not args.after.endswith(".bim"):
        msg = "%s: not a BIM file (extension must be .bim)" % args.after
        raise ProgramError(msg)
    elif not os.path.isfile(args.after):
        msg = "%s: no such file" % args.after
        raise ProgramError(msg)

    return True


def parseArgs():  # pragma: no cover
    """Parses the command line options and arguments.

    :returns: A :py:class:`argparse.Namespace` object created by the
              :py:mod:`argparse` module. It contains the values of the
              different options.

    ============ ====== =============================================
       Options    Type                 Description
    ============ ====== =============================================
    ``--before`` string The name of the BIM file before modification.
    ``--after``  string The name of the BIM file after modification.
    ``--out``    string The prefix of the output files
    ============ ====== =============================================

    .. note::
        No option check is done here (except for the one automatically done by
        argparse). Those need to be done elsewhere (see :py:func:`checkArgs`).

    """
    # The INPUT files
    group = parser.add_argument_group("Input File")
    group.add_argument("--before", type=str, metavar="FILE", required=True,
                       help="The name of the bim FILE before modification.")
    group.add_argument("--after", type=str, metavar="FILE", required=True,
                       help="The name of the bim FILE after modification.")

    # The OUTPUT files
    group = parser.add_argument_group("Output File")
    group.add_argument("--out", type=str, metavar="FILE",
                       default="snp_removed",
                       help="The prefix of the output files. [default: "
                            "%(default)s]")

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
