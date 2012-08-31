#!/usr/bin/env python2.7

import os
import sys
import argparse
import subprocess

def main(argString=None):
    # Getting and checking the options
    args = parseArgs(argString)
    checkArgs(args)

    print "   - Options used:"
    for key, value in vars(args).iteritems():
        print "      --{} {}".format(key, value)

    # Running Plink
    print "   - Running Plink"
    runPlink(args)


def runPlink(options):
    """Run Plink with the following options."""
    # The plink command
    plinkCommand = ["plink", "--noweb",
                    "--bfile" if options.is_bfile else "--tfile", options.ifile,
                    "--mind", str(options.mind), "--make-bed", "--out",
                    options.out]

    output = None
    try:
        output = subprocess.check_output(plinkCommand,
                                         stderr=subprocess.STDOUT, shell=False)
    except subprocess.CalledProcessError:
        msg = "plink: couldn't run plink"
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
    required_file_extensions = {".tfam", ".tped"}
    if args.is_bfile:
        required_file_extensions = {".bed", ".bim", ".fam"}
    for fileName in [args.ifile + i for i in required_file_extensions]:
        if not os.path.isfile(fileName):
            msg = "{}: no such file".format(fileName)
            raise ProgramError(msg)

    # Check the mind option (between 0 and 1, inclusive)
    if (args.mind < 0) or (args.mind > 1):
        msg = "mind=%f: must be between 0 and 1 (inclusive)" % args.mind
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
prog = "sample_missingess"
desc = """Comput sample missingness using Plink"""
parser = argparse.ArgumentParser(description=desc, prog=prog)

# The INPUT files
group = parser.add_argument_group("Input File")
group.add_argument("--ifile", type=str, metavar="FILE", required=True,
                   help=("The input file prefix (by default, this input file "
                         "must be a tfile. If options --is-bfile is used, the "
                         "input file must be a bfile)."))
# The options
group = parser.add_argument_group("Options")
group.add_argument("--is-bfile", action="store_true",
                   help=("The input file (--ifile) is a bfile instead of a "
                         "tfile."))
group.add_argument("--mind", type=float, metavar="FLOAT", default=0.1,
                   help=("The missingness threshold (remove samples with "
                         "more than x percent missing genotypes). [Default: "
                         "%(default).3f]"))
# The OUTPUT files
group = parser.add_argument_group("Output File")
group.add_argument("--out", type=str, metavar="FILE", default="clean_mind",
                   help=("The prefix of the output files (wich will be a Plink "
                         "binary file). [default: %(default)s]"))

if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print >>sys.stderr, "Cancelled by user"
        sys.exit(0)
    except ProgramError as e:
        parser.error(e.message)
