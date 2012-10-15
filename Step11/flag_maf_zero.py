#!/usr/bin/env python2.7

import os
import sys
import argparse
import subprocess

from PlinkUtils import createRowFromPlinkSpacedOutput

def main(argString=None):
    # Getting and checking the options
    args = parseArgs(argString)
    checkArgs(args)

    print "   - Options used:"
    for key, value in vars(args).iteritems():
        print "      --{} {}".format(key, value)

    # Compute frequency using plink
    print "   - Computing the frequencies using Plink"
    computeFrequency(args)

    # Read the freqency file
    print "   - Flagging SNPs with MAF = 0"
    findSnpWithMaf0(args.out + ".frq", args.out)


def findSnpWithMaf0(freqFileName, prefix):
    """Finds SNPs with MAF of 0 and put them in a file."""
    maf_0_set = set()
    na_set = set()

    try:
        with open(freqFileName, "r") as inputFile:
            headerIndex = None
            for i, line in enumerate(inputFile):
                row = createRowFromPlinkSpacedOutput(line)
                if i == 0:
                    # We have the header
                    headerIndex = dict([(row[i], i) \
                                            for i in xrange(len(row))])
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
        print "      - There are no markers with MAF 0"
    else:
        print "      - There are {} markers with MAF 0".format(len(maf_0_set))
    outputFile = None
    try:
        with open(prefix + ".list", "w") as output_file:
            for marker_name in maf_0_set:
                print >>output_file, marker_name
    except IOError:
        msg = "{}.list: can't write file".format(prefix)
        raise ProgramError(msg)

    if len(na_set) > 0:
        print "      - There are {} markers with NA MAF".format(len(na_set))
        try:
            with open(prefix + ".na_list", "w") as output_file:
                for marker_name in na_set:
                    print >>output_file, marker_name
        except IOError:
            msg = "{}.na_list: can't write file".format(prefix)
            raise ProgramError(msg)



def computeFrequency(options):
    """Compute the frequency of the SNPs."""
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

    :param args: a :py:class:`Namespace` object containing the options of the
                 program.
    :type args: :py:class:`argparse.Namespace`

    :returns: ``True`` if everything was OK.

    If there is a problem with an option, an exception is raised using the
    :py:class:`ProgramError` class, a message is printed
    to the :class:`sys.stderr` and the program exists with code 1.

    """
    # Check if we have the tped and the tfam files
    for fileName in [args.bfile + i for i in [".bed", ".bim", ".fam"]]:
        if not os.path.isfile(fileName):
            msg = "%(fileName)s: no such file" % locals()
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
prog = "flag_maf_zero"
desc = """Flag SNPs with MAF of 0."""
parser = argparse.ArgumentParser(description=desc, prog=prog)

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

if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print >>sys.stderr, "Cancelled by user"
        sys.exit(0)
    except ProgramError as e:
        parser.error(e.message)
