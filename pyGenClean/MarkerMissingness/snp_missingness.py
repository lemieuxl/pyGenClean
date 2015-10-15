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


logger = logging.getLogger("snp_missingness")


def main(argString=None):
    """The main function of the module.

    :param argString: the options.

    :type argString: list

    These are the steps:

    1. Prints the options.
    2. Runs Plink with the ``geno`` option (:py:func:`runPlink`).
    3. Compares the two ``bim`` files (before and after the Plink ``geno``
       analysis) (:py:func:`compareBIM`).

    """
    # Getting and checking the options
    args = parseArgs(argString)
    checkArgs(args)

    logger.info("Options used:")
    for key, value in vars(args).iteritems():
        logger.info("  --{} {}".format(key.replace("_", "-"), value))

    # Run plink
    logger.info("Running Plink")
    runPlink(args)

    # Comparing the bim
    logger.info("Comparing BIM files")
    compareBIM(args)


def compareBIM(args):
    """Compare two BIM file.

    :param args: the options.

    :type args: argparse.Namespace

    Creates a *Dummy* object to mimic an :py:class:`argparse.Namespace` class
    containing the options for the :py:mod:`pyGenClean.PlinkUtils.compare_bim`
    module.

    """
    # Creating the CompareBIM options
    class Dummy(object):
        pass
    compareBIM_args = Dummy()
    compareBIM_args.before = args.bfile + ".bim"
    compareBIM_args.after = args.out + ".bim"
    compareBIM_args.out = args.out + ".removed_snps"

    try:
        # Checking the arguments
        CompareBIM.checkArgs(compareBIM_args)

        # Reading the BIM files
        beforeBIM = CompareBIM.readBIM(compareBIM_args.before)
        afterBIM = CompareBIM.readBIM(compareBIM_args.after)

        # Comparing the BIM files
        CompareBIM.compareSNPs(beforeBIM, afterBIM, compareBIM_args.out)
    except CompareBIM.ProgramError as e:
        raise ProgramError("CompareBIM: " + e.message)


def runPlink(options):
    """Runs Plink with the geno option.

    :param options: the options.

    :type options: argparse.Namespace

    """
    # The plink command
    plinkCommand = ["plink", "--noweb", "--bfile", options.bfile, "--geno",
                    str(options.geno), "--make-bed", "--out", options.out]

    output = None
    try:
        output = subprocess.check_output(plinkCommand,
                                         stderr=subprocess.STDOUT, shell=False)
    except subprocess.CalledProcessError:
        msg = "plink: couldn't run plink"
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
    # Check if we have the tped and the tfam files
    for fileName in [args.bfile + i for i in [".bed", ".bim", ".fam"]]:
        if not os.path.isfile(fileName):
            msg = "%(fileName)s: no such file" % locals()
            raise ProgramError(msg)

    # Check the mind option (between 0 and 1, inclusive)
    if (args.geno < 0) or (args.geno > 1):
        msg = "mind=%f: must be between 0 and 1 (inclusive)" % args.geno
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
      Options    Type                Description
    =========== ====== ==========================================
    ``--bfile`` string The input file prefix (Plink binary file).
    ``--geno``  float  The missingness threshold.
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
pretty_name = "Marker missingness"
desc = "Computes marker missingness using Plink."
long_desc = ("The script identifies poorly performing markers with genotyping "
             r"succes rate $<{success_rate}$\%.")
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
group.add_argument("--geno", type=float, metavar="FLOAT", default=0.02,
                   help=("The missingness threshold (remove SNPs with more "
                         "than x percent missing genotypes). [Default: "
                         "%(default).3f]"))
# The OUTPUT files
group = parser.add_argument_group("Output File")
group.add_argument("--out", type=str, metavar="FILE",
                   default="clean_geno",
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
