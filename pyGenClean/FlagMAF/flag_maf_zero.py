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
import subprocess

from .. import __version__
from ..PlinkUtils import createRowFromPlinkSpacedOutput


logger = logging.getLogger("flag_maf_zero")


def main(argString=None):
    """The main function.

    :param argString: the options.

    :type argString: list

    These are the steps:

    1. Prints the options.
    2. Computes the frequencies using Plinl (:py:func:`computeFrequency`).
    3. Finds markers with MAF of 0, and saves them in a file
       (:py:func:`findSnpWithMaf0`).

    """
    # Getting and checking the options
    args = parseArgs(argString)
    checkArgs(args)

    logger.info("Options used:")
    for key, value in vars(args).iteritems():
        logger.info("  --{} {}".format(key.replace("_", "-"), value))

    # Compute frequency using plink
    logger.info("Computing the frequencies using Plink")
    computeFrequency(args)

    # Read the freqency file
    logger.info("Flagging SNPs with MAF = 0")
    findSnpWithMaf0(args.out + ".frq", args.out)


def findSnpWithMaf0(freqFileName, prefix):
    """Finds SNPs with MAF of 0 and put them in a file.

    :param freqFileName: the name of the frequency file.
    :param prefix: the prefix of all the files.

    :type freqFileName: str
    :type prefix: str

    Reads a frequency file from Plink, and find markers with a minor allele
    frequency of zero.

    """
    maf_0_set = set()
    na_set = set()

    try:
        with open(freqFileName, "r") as inputFile:
            headerIndex = None
            for i, line in enumerate(inputFile):
                row = createRowFromPlinkSpacedOutput(line)
                if i == 0:
                    # We have the header
                    headerIndex = dict([
                        (row[i], i) for i in xrange(len(row))
                    ])
                    for columnName in ["SNP", "MAF"]:
                        if columnName not in headerIndex:
                            msg = "%(freqFileName)s: no column named " \
                                  "%(columnName)s" % locals()
                            raise ProgramError(msg)

                else:
                    # We have data
                    snpName = row[headerIndex["SNP"]]
                    snpMAF = row[headerIndex["MAF"]]

                    if snpMAF == "0":
                        # We want to flag this SNP
                        maf_0_set.add(snpName)

                    elif snpMAF == "NA":
                        # We want to flag this SNP, because the MAF est NA
                        na_set.add(snpName)

    except IOError:
        msg = "%(freqFileName)s: no such file" % locals()
        raise ProgramError(msg)

    # Creating the output files
    if len(maf_0_set) == 0:
        logger.info("  - There are no markers with MAF 0")
    else:
        logger.info("  - There are {} markers with MAF 0".format(
            len(maf_0_set),
        ))
    outputFile = None
    try:
        with open(prefix + ".list", "w") as output_file:
            for marker_name in maf_0_set:
                print >>output_file, marker_name
    except IOError:
        msg = "{}.list: can't write file".format(prefix)
        raise ProgramError(msg)

    if len(na_set) > 0:
        logger.info("  - There are {} markers with NA MAF".format(len(na_set)))
        try:
            with open(prefix + ".na_list", "w") as output_file:
                for marker_name in na_set:
                    print >>output_file, marker_name
        except IOError:
            msg = "{}.na_list: can't write file".format(prefix)
            raise ProgramError(msg)


def computeFrequency(options):
    """Compute the frequency of the SNPs.

    :param options: the options.

    :type options: argparse.Namespace

    """
    # The plink command
    plinkCommand = ["plink", "--noweb", "--bfile", options.bfile,
                    "--freq", "--out", options.out]

    output = None
    try:
        output = subprocess.check_output(plinkCommand,
                                         stderr=subprocess.STDOUT, shell=False)
    except subprocess.CalledProcessError:
        msg = "couldn't run command\n" + " ".join(plinkCommand)
        raise ProgramError(msg)


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
    # Check if we have the tped and the tfam files
    for fileName in [args.bfile + i for i in [".bed", ".bim", ".fam"]]:
        if not os.path.isfile(fileName):
            msg = "%(fileName)s: no such file" % locals()
            raise ProgramError(msg)

    return True


def parseArgs(argString=None):  # pragma: no cover
    """Parses the command line options and arguments.

    :param argString: the options.

    :type argString: list

    :returns: A :py:class:`argparse.Namespace` object created by the
              :py:mod:`argparse` module. It contains the values of the
              different options.

    =========== ====== ==========================================
      Options    Type                 Description
    =========== ====== ==========================================
    ``--bfile`` string The input file prefix (Plink binary file).
    ``--out``   string The prefix of the output files.
    =========== ====== ==========================================

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
pretty_name = "MAF flagging"
desc = "Flags SNPs with MAF of 0."
long_desc = ("The script flags markers with a minor allele frequency of 0 "
             r"\textit{i.e} monomorphic markers).")
parser = argparse.ArgumentParser(description=desc)
parser.add_argument("-v", "--version", action="version",
                    version="pyGenClean version {}".format(__version__))

# The INPUT files
group = parser.add_argument_group("Input File")
group.add_argument("--bfile", type=str, metavar="FILE", required=True,
                   help=("The input file prefix (will find the plink binary "
                         "files by appending the prefix to the .bim, .bed and "
                         ".fam files, respectively."))
# The OUTPUT files
group = parser.add_argument_group("Output File")
group.add_argument("--out", type=str, metavar="FILE", default="flag_maf_0",
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
