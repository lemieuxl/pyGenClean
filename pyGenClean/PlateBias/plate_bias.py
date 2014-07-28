#!/usr/bin/env python2.7
## This file is part of pyGenClean.
## 
## pyGenClean is free software: you can redistribute it and/or modify it under
## the terms of the GNU General Public License as published by the Free Software
## Foundation, either version 3 of the License, or (at your option) any later
## version.
## 
## pyGenClean is distributed in the hope that it will be useful, but WITHOUT ANY
## WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
## A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
## 
## You should have received a copy of the GNU General Public License along with
## pyGenClean.  If not, see <http://www.gnu.org/licenses/>.

import os
import sys
import glob
import argparse
import subprocess

from PlinkUtils import createRowFromPlinkSpacedOutput

def main(argString=None):
    """The main function of this module.

    :param argString: the options.

    :type argString: list of strings

    These are the steps:

    1. Runs a plate bias analysis using Plink
       (:py:func:`executePlateBiasAnalysis`).
    2. Extracts the list of significant markers after plate bias analysis
       (:py:func:`extractSignificantSNPs`).
    3. Computes the frequency of all significant markers after plate bias
       analysis (:py:func:`computeFrequencyOfSignificantSNPs`).

    """
    # Getting and checking the options
    args = parseArgs(argString)
    checkArgs(args)

    print "   - Options used:"
    for key, value in vars(args).iteritems():
        print "      --{} {}".format(key, value)
    
    # Run plink
    print "   - Running Plink to check the plate bias"
    executePlateBiasAnalysis(args)

    # Extract significant SNPs
    print "   - Extracting significant SNPs"
    extractSignificantSNPs(args.out)

    # Remove significant SNPs using plink
    print "   - Computing frequency of significant SNPs"
    computeFrequencyOfSignificantSNPs(args)


def extractSignificantSNPs(prefix):
    """Extract significant SNPs in the fisher file.

    :param prefix: the prefix of the input file.

    :type prefix: string

    Reads a list of significant markers (``prefix.assoc.fisher``) after plate
    bias analysis with Plink. Writes a file (``prefix.significant_SNPs.txt``)
    containing those significant markers.

    """
    fileNames = glob.glob(prefix + ".*.assoc.fisher")

    snpNames = set()
    for fileName in fileNames:
        try:
            with open(fileName, 'r') as inputFile:
                headerIndex = None
                for line in inputFile:
                    row = createRowFromPlinkSpacedOutput(line)
                    
                    if headerIndex is None:
                        # This is the header line
                        headerIndex = dict([(row[i], i) \
                                                for i in xrange(len(row))])
                    else:
                        snpName = row[headerIndex["SNP"]]
                        snpNames.add(snpName)
        except IOError:
            msg = "%(fileName)s: no such file" % locals()
            raise ProgramError(msg)

    # The output file
    outputFileName = prefix + ".significant_SNPs.txt"
    outputFile = None
    try:
        outputFile = open(outputFileName, "w")
    except IOError:
        msg = "%(outputFileName)s: can't write file" % locals()
        raise ProgramError(msg)

    print >>outputFile, "\n".join(snpNames)

    # Closing the file
    outputFile.close()


def computeFrequencyOfSignificantSNPs(options):
    """Computes the frequency of the significant markers.

    :param options: the options.

    :type options: argparse.Namespace

    Extract a list of markers (significant after plate bias analysis) and
    computes their frequencies.

    """
    # The plink command
    plinkCommand = ["plink", "--noweb", "--bfile", options.bfile, "--extract",
                    options.out + ".significant_SNPs.txt", "--freq", "--out",
                    options.out + ".significant_SNPs"]
    runCommand(plinkCommand)
    

def executePlateBiasAnalysis(options):
    """Execute the plate bias analysis with Plink.

    :param options: the options.

    :type options: argparse.Namespace

    """
    # The plink command
    plinkCommand = ["plink", "--noweb", "--bfile", options.bfile,
                    "--loop-assoc", options.loop_assoc, "--fisher",
                    "--pfilter", str(options.pfilter), "--out", options.out]
    runCommand(plinkCommand)


def runCommand(command):
    """Run a command.

    :param command: the command to run.

    :type command: list of strings

    Tries to run a command. If it fails, raise a :py:class:`ProgramError`. This
    function uses the :py:mod:`subprocess` module.

    .. warning::
        The variable ``command`` should be a list of strings (no other type).

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

    # Check the plate bias file
    if not os.path.isfile(args.loop_assoc):
        msg = "%s: no such file" % args.loop_assoc
        raise ProgramError(msg)

    return True


def parseArgs(argString=None): # pragma: no cover
    """Parses the command line options and arguments.

    :param argString: the options.

    :type argString: list of strings

    :returns: A :py:class:`argparse.Namespace` object created by the
              :py:mod:`argparse` module. It contains the values of the different
              options.

    ================ ====== ==================================================
        Options       Type                    Description
    ================ ====== ==================================================
    ``--bfile``      string The input file prefix (Plink binary).
    ``--loop-assoc`` string The file containing the plate organization of each
                            samples.
    ``--pfilter``    float  The significance threshold used for the plate
                            effect.
    ``--out``        string The prefix of the output files.
    ================ ====== ==================================================

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
desc = """Check for plate bias."""
parser = argparse.ArgumentParser(description=desc)

# The INPUT files
group = parser.add_argument_group("Input File")
group.add_argument("--bfile", type=str, metavar="FILE", required=True,
                   help=("The input file prefix (will find the plink binary "
                         "files by appending the prefix to the .bim, .bed "
                         "and .fam files, respectively."))
group.add_argument("--loop-assoc", type=str, metavar="FILE", required=True,
                    help=("The file containing the plate organization of "
                          "each samples. Must contains three column (with no "
                          "header): famID, indID and plateName."))
# The options
group = parser.add_argument_group("Options")
group.add_argument("--pfilter", type=float, metavar="FLOAT", default=1e-7,
                    help=("The significance threshold used for the plate "
                          "effect. [default: %(default).1e]"))
# The OUTPUT files
group = parser.add_argument_group("Output File")
group.add_argument("--out", type=str, metavar="FILE",
                    default="plate_bias",
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