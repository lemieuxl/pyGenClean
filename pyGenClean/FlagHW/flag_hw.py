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
from ..PlinkUtils import compare_bim as CompareBIM
from ..PlinkUtils import createRowFromPlinkSpacedOutput


logger = logging.getLogger("flag_hw")


class Dummy(object):
    pass


def main(argString=None):
    """The main function.

    :param argString: the options.

    :type argString: list

    These are the steps performed by this module:

    1. Prints the options of the module.
    2. Computes the number of markers in the input file
       (:py:func:`computeNumberOfMarkers`).
    3. If there are no markers, the module stops.
    4. Computes the Bonferroni therhold (:math:`0.05 / \\textrm{nbMarkers}`).
    5. Runs Plink to find failed markers with the Bonferroni threshold.
    6. Runs Plink to find failed markers with the default threshold.
    7. Compares the ``bim`` files for the Bonferroni threshold.
    8. Compares the ``bim`` files for the default threshold.
    9. Computes the "in between" marker list, which is the markers from the
       default threshold and the Bonferroni threshold.

    """
    # Getting and checking the options
    args = parseArgs(argString)
    checkArgs(args)

    logger.info("Options used:")
    for key, value in vars(args).iteritems():
        logger.info("  --{} {}".format(key.replace("_", "-"), value))

    # Compute the number of markers
    logger.info("Counting the number of markers")
    nbMarkers = computeNumberOfMarkers(args.bfile + ".bim")

    if nbMarkers <= 0:
        logger.info("  - There are no markers: STOPPING NOW!")
    else:
        logger.info("  - There are {} markers".format(nbMarkers))

        customThreshold = str(0.05 / nbMarkers)

        # Run the plink command
        logger.info("Computing the HW equilibrium for {}".format(
            customThreshold,
        ))
        computeHWE(args.bfile, customThreshold,
                   args.out + ".threshold_" + customThreshold)
        logger.info("Computing the HW equilibrium for {}".format(args.hwe))
        computeHWE(args.bfile, args.hwe, args.out + ".threshold_" + args.hwe)

        # Compare the BIM files
        logger.info("Creating the flagged SNP list for {}".format(
            customThreshold,
        ))
        custom_snps = compareBIMfiles(
            args.bfile + ".bim",
            args.out + ".threshold_" + customThreshold + ".bim",
            args.out + ".snp_flag_threshold_" + customThreshold,
        )
        logger.info("Creating the flagged SNP list for {}".format(args.hwe))
        hwe_snps = compareBIMfiles(
            args.bfile + ".bim",
            args.out + ".threshold_" + args.hwe + ".bim",
            args.out + ".snp_flag_threshold_" + args.hwe,
        )

        logger.info("Creating the in between SNP list ([{}, {}[)".format(
            args.hwe,
            customThreshold,
        ))
        file_name = args.out + ".snp_flag_threshold_between_{}-{}".format(
            args.hwe,
            customThreshold,
        )
        try:
            with open(file_name, 'w') as output_file:
                differences = hwe_snps - custom_snps
                if len(differences) > 0:
                    print >>output_file, "\n".join(differences)
        except IOError:
            msg = "{}: can't write file".format(file_name)
            raise ProgramError(msg)


def compareBIMfiles(beforeFileName, afterFileName, outputFileName):
    """Compare two BIM files for differences.

    :param beforeFileName: the name of the file before modification.
    :param afterFileName: the name of the file after modification.
    :param outputFileName: the name of the output file (containing the
                           differences between the ``before`` and the ``after``
                           files.

    :type beforeFileName: str
    :type afterFileName: str
    :type outputFileName: str

    :returns: the number of differences between the two files.

    The ``bim`` files contain the list of markers in a given dataset. The
    ``before`` file should have more markers than the ``after`` file. The
    ``after`` file should be a subset of the markers in the ``before`` file.

    """
    # Creating the options
    options = Dummy()
    options.before = beforeFileName
    options.after = afterFileName
    options.out = outputFileName

    # Checking the options
    CompareBIM.checkArgs(options)

    # Reading the BIM files
    beforeBIM = CompareBIM.readBIM(options.before)
    afterBIM = CompareBIM.readBIM(options.after)

    # Finding the differences
    CompareBIM.compareSNPs(beforeBIM, afterBIM, options.out)

    return beforeBIM - afterBIM


def computeNumberOfMarkers(inputFileName):
    """Count the number of marker (line) in a BIM file.

    :param inputFileName: the name of the ``bim`` file.

    :type inputFileName: str

    :returns: the number of marker in the ``bim`` file.

    """
    nbLine = 0
    with open(inputFileName, "r") as inputFile:
        nbLine = len(inputFile.readlines())

    return nbLine


def computeHWE(prefix, threshold, outPrefix):
    """Compute the Hardy Weinberg test using Plink.

    :param prefix: the prefix of all the files.
    :param threshold: the Hardy Weinberg threshold.
    :param outPrefix: the prefix of the output file.

    :type prefix: str
    :type threshold: str
    :type outPrefix: str

    Uses Plink to exclude markers that failed the Hardy-Weinberg test at a
    specified significance threshold.

    """
    plinkCommand = ["plink", "--noweb", "--bfile", prefix, "--hwe", threshold,
                    "--make-bed", "--out", outPrefix]
    runCommand(plinkCommand)


def runCommand(command):
    """Run a command.

    :param command: the command to run.

    :type command: list

    Tries to run a command. If it fails, raise a :py:class:`ProgramError`. This
    function uses the :py:mod:`subprocess` module.

    .. warning::
        The variable ``command`` should be a list of stings (no other type).

    """
    output = None
    try:
        output = subprocess.check_output(command,
                                         stderr=subprocess.STDOUT, shell=False)
    except subprocess.CalledProcessError:
        msg = "couldn't run command\n" + " ".join(command)
        raise ProgramError(msg)


def checkArgs(args):
    """Checks the arguments and options.

    :param args: a :py:class:`argparse.Namespace` object containing the options
                 of the program.

    :type args: :py:class:`argparse.Namespace`

    :returns: ``True`` if everything was OK.

    If there is a problem with an option, an exception is raised using the
    :py:class:`ProgramError` class, a message is printed to the
    :py:class:`sys.stderr` and the program exists with code 1.

    """
    # Check if we have the tped and the tfam files
    for fileName in [args.bfile + i for i in [".bed", ".bim", ".fam"]]:
        if not os.path.isfile(fileName):
            msg = "%(fileName)s: no such file" % locals()
            raise ProgramError(msg)

    # Check the HW pvalue
    hwValue = args.hwe
    try:
        hwValue = float(hwValue)
    except ValueError:
        msg = "%(hwValue)s: not a valid HW value" % locals()
        raise ProgramError(msg)
    if (hwValue < 0) or (hwValue > 1):
        msg = "%s: not a valid HW value" % args.hwe
        raise ProgramError(msg)

    return True


def parseArgs(argString=None):  # pragma: no cover
    """Parses the command line options and arguments.

    :param argString: the options.

    :type argString: list

    :returns: A :py:class:`argparse.Namespace` object created by
              the :py:mod:`argparse` module. It contains the values of the
              different options.

    =========== ====== ==========================================
      Options    Type                     Description
    =========== ====== ==========================================
    ``--bfile`` string The input file prefix (binary Plink file).
    ``--hwe``   float  The Hardy-Weinberg equilibrium threshold.
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
pretty_name = "Hardy-Weinberg flagging"
desc = "Flags SNPs with Hardy-Weinberg disequilibrium."
long_desc = ("The script tests for Hardy-Weinberg equilibrium for each marker "
             "(using an exact test). It adjusts for multiple testing using "
             "Bonferroni.")
parser = argparse.ArgumentParser(description=desc)
parser.add_argument("-v", "--version", action="version",
                    version="pyGenClean version {}".format(__version__))

# The INPUT files
group = parser.add_argument_group("Input File")
group.add_argument("--bfile", type=str, metavar="FILE", required=True,
                   help=("The input file prefix (will find the plink binary "
                         "files by appending the prefix to the .bim, .bed and "
                         ".fam files, respectively."))
# The options
group = parser.add_argument_group("Options")
group.add_argument("--hwe", type=str, metavar="FLOAT", default="1e-4",
                   help=("The Hardy-Weinberg equilibrium threshold. "
                         "[default: %(default)s]"))
# The OUTPUT files
group = parser.add_argument_group("Output File")
group.add_argument("--out", type=str, metavar="FILE", default="flag_hw",
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
