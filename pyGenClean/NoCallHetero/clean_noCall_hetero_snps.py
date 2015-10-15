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
import shutil
import logging
import argparse

import numpy as np

from .. import __version__


logger = logging.getLogger("noCall_hetero_snps")


def main(argString=None):
    """The main function of the module.

    :param argString: the options.

    :type argString: list

    These are the steps:

    1. Prints the options.
    2. Reads the ``tfam`` and ``tped`` files and find all heterozygous and all
       failed markers (:py:func:`processTPEDandTFAM`).

    """
    # Getting and checking the options
    args = parseArgs(argString)
    checkArgs(args)

    logger.info("Options used:")
    for key, value in vars(args).iteritems():
        logger.info("  --{} {}".format(key.replace("_", "-"), value))

    # Process the TPED and TFAM file
    logger.info("Processing the TPED and TFAM file")
    processTPEDandTFAM(args.tfile + ".tped", args.tfile + ".tfam", args.out)


def processTPEDandTFAM(tped, tfam, prefix):
    """Process the TPED and TFAM files.

    :param tped: the name of the ``tped`` file.
    :param tfam: the name of the ``tfam`` file.
    :param prefix: the prefix of the output files.

    :type tped: str
    :type tfam: str
    :type prefix: str

    Copies the original ``tfam`` file into ``prefix.tfam``. Then, it reads the
    ``tped`` file and keeps in memory two sets containing the markers which are
    all failed or which contains only heterozygous genotypes.

    It creates two output files, ``prefix.allFailed`` and ``prefix.allHetero``,
    containing the markers that are all failed and are all heterozygous,
    respectively.

    .. note::
        All heterozygous markers located on the mitochondrial chromosome are
        not remove.

    """
    # Copying the tfam file
    try:
        shutil.copy(tfam, prefix + ".tfam")
    except IOError:
        msg = "%s: can't write file" % prefix + ".tfam"
        raise ProgramError(msg)

    try:
        outputFile = open(prefix + ".tped", "w")
    except IOError:
        msg = "%s: can't write to file" % prefix + ".tped"
        raise ProgramError(msg)

    # The name of the bad SNPs
    allFailed = set()
    allHetero = set()

    with open(tped, 'r') as inputFile:
        for line in inputFile:
            row = line.rstrip("\r\n").split("\t")
            snpInfo = row[:4]
            chromosome = snpInfo[0]
            genotypes = np.array([i.upper() for i in row[4:]])

            # Testing the genotypes
            uniqueGenotypes = np.unique(genotypes)
            if len(uniqueGenotypes) == 1:
                # We have only one kind of genotype, either all homo, all
                # hetero or all no call
                if uniqueGenotypes[0] == "0 0":
                    # This one is a no call
                    allFailed.add(snpInfo[1])
                elif len(set(uniqueGenotypes[0].split(" "))) == 2:
                    # There are two different alleles, hence, hetero
                    if chromosome not in {"26", "MT"}:
                        # The SNP is not on a mitochondrial chromosome
                        # (because we want to keep those)
                        allHetero.add(snpInfo[1])

            # If the SNP is good (neither in allFailed or allHetero), we keep
            # the SNP
            if (snpInfo[1] not in allFailed) and (snpInfo[1] not in allHetero):
                print >>outputFile, "\t".join(row)

    outputFile.close()

    # Now printing the summary files
    # The SNPs with no calls
    fileName = prefix + ".allFailed"
    try:
        with open(fileName, 'w') as outputFile:
            if len(allFailed) > 0:
                print >>outputFile, "\n".join(allFailed)
    except IOError:
        msg = "%(fileName)s: can't write to file" % locals()
        raise ProgramError(msg)

    # The SNPs with only hetero calls
    fileName = prefix + ".allHetero"
    try:
        with open(fileName, 'w') as outputFile:
            if len(allHetero) > 0:
                print >>outputFile, "\n".join(allHetero)
    except IOError:
        msg = "%(fileName)s: can't write to file" % locals()
        raise ProgramError(msg)


def checkArgs(args):
    """Checks the arguments and options.

    :param args: an object containing the options of the program.

    :type args: argparse.Namespace

    :returns: ``True`` if everything was OK.

    If there is a problem with an option, an exception is raised using the
    :py:class:`ProgramError` class, a message is printed to the
    :class:`sys.stderr` and the program exists with code 1.

    """
    # Checking the input files
    for suffix in [".tped", ".tfam"]:
        fileName = args.tfile + suffix
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

    =========== ====== ====================================
      Options    Type              Description
    =========== ====== ====================================
    ``--tfile`` string The input file prefix (Plink tfile).
    ``--out``   string The prefix of the output files
    =========== ====== ====================================

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
pretty_name = "No calls and heterozygous only markers"
desc = 'Removes "no calls" only and heterozygous only markers.'
long_desc = ("The script removes completely failed markers (no calls) or "
             "markers with only heterozygous genotypes.")
parser = argparse.ArgumentParser(description=desc)
parser.add_argument("-v", "--version", action="version",
                    version="pyGenClean version {}".format(__version__))

# The INPUT files
group = parser.add_argument_group("Input File")
group.add_argument("--tfile", type=str, metavar="FILE", required=True,
                   help=("The input file prefix (will find the tped and tfam "
                         "file by appending the prefix to .tped and .tfam, "
                         "respectively."))
# The OUTPUT files
group = parser.add_argument_group("Output File")
group.add_argument("--out", type=str, metavar="FILE",
                   default="clean_noCall_hetero",
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
