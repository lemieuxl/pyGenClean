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
import argparse
import subprocess

def main(argString=None):
    """The main function of the module.

    :param argString: the options.

    :type argString: list of strings

    These are the steps:

    1. Prints the options.
    2. Runs Plink with the ``mind`` option (:py:func:`runPlink`).

    """
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
    """Run Plink with the ``mind`` option.

    :param options: the options.

    :type options: argparse.Namespace

    """
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

    :param args: a object containing the options of the program.

    :type args: argparse.Namespace

    :returns: ``True`` if everything was OK.

    If there is a problem with an option, an exception is raised using the
    :py:class:`ProgramError` class, a message is printed to the
    :class:`sys.stderr` and the program exists with code 1.

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

    :param argString: the options.

    :type argString: list of strings

    :returns: A :py:class:`argparse.Namespace` object created by the
              :py:mod:`argparse` module. It contains the values of the different
              options.

    ============== ====== ======================================================
        Options     Type                     Description
    ============== ====== ======================================================
    ``--ifile``    string The input file prefix (either a Plink binary file or a
                          tfile).
    ``--is-bfile`` bool   The input file (``--ifile``) is a bfile instead of a
                          tfile.
    ``--mind``     float  The missingness threshold.
    ``--out``      string The prefix of the output files.
    ============== ====== ======================================================

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
desc = """Computes sample missingness using Plink"""
parser = argparse.ArgumentParser(description=desc)

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
