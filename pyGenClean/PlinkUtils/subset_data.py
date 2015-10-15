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


logger = logging.getLogger("subset")


def main(argString=None):
    """The main function of the modile.

    :param argString: the options.

    :type argString: list

    Here are the steps:

    1. Prints the options.
    2. Subset the data (:py:func:`subset_data`).

    .. note::
        The type of the output files are determined by the type of the input
        files (*e.g.* if the input files are binary files, so will be the
        output ones).

    """
    # Getting and checking the options
    args = parseArgs(argString)
    checkArgs(args)

    logger.info("Options used:")
    for key, value in vars(args).iteritems():
        logger.info("  --{} {}".format(key.replace("_", "-"), value))

    # Subset the data
    logger.info("Subsetting the data using Plink")
    subset_data(args)


def subset_data(options):
    """Subset the data.

    :param options: the options.

    :type options: argparse.Namespace

    Subset the data using either ``--exclude`` or ``--extract``for markers or
    ``--remove`` or ``keep`` for samples.

    """
    # The plink command
    plinkCommand = ["plink", "--noweb"]

    # The input file prefix
    if options.is_bfile:
        plinkCommand += ["--bfile", options.ifile, "--make-bed"]
    elif options.is_tfile:
        plinkCommand += ["--tfile", options.ifile, "--recode", "--transpose",
                         "--tab"]
    elif options.is_file:
        plinkCommand += ["--file", options.ifile, "--recode", "--tab"]

    # The subset command for SNPs
    if options.exclude is not None:
        plinkCommand += ["--exclude", options.exclude]
    elif options.extract is not None:
        plinkCommand += ["--extract", options.extract]

    # The subset command for samples
    if options.keep is not None:
        plinkCommand += ["--keep", options.keep]
    elif options.remove is not None:
        plinkCommand += ["--remove", options.remove]

    # The output prefix
    plinkCommand += ["--out", options.out]

    # Running the command
    runCommand(plinkCommand)


def runCommand(command):
    """Runs a command.

    :param command: the command to run.

    :type command: list

    If there is a problem, a :py:class:`ProgramError` is raised.

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

    .. note::
        Only one operation for markers and one operation for samples can be
        done at a time. Hence, one of ``--exclude`` or ``--extract`` can be
        done for markers, and one of ``--remove`` or ``--keep`` can be done for
        samples.

    """
    # Check the input files
    if not args.is_bfile and not args.is_tfile and not args.is_file:
        msg = "needs one input file type (--is-bfile, --is-tfile or --is-file)"
        raise ProgramError(msg)

    if args.is_bfile and not args.is_tfile and not args.is_file:
        for fileName in [args.ifile + i for i in [".bed", ".bim", ".fam"]]:
            if not os.path.isfile(fileName):
                msg = "{}: no such file".format(fileName)
                raise ProgramError(msg)

    elif args.is_tfile and not args.is_bfile and not args.is_file:
        for fileName in [args.ifile + i for i in [".tped", ".tfam"]]:
            if not os.path.isfile(fileName):
                msg = "{}: no such file".format(fileName)
                raise ProgramError(msg)

    elif args.is_file and not args.is_bfile and not args.is_tfile:
        for fileName in [args.ifile + i for i in [".ped", ".map"]]:
            if not os.path.isfile(fileName):
                msg = "{}: no such file". format(fileName)
                raise ProgramError(msg)

    else:
        msg = ("needs only one input file type (--is-bfile, --is-tfile or "
               "--is-file)")
        raise ProgramError(msg)

    # Check that we have at least one of exclude, extract remove or keep
    if args.exclude is None and args.extract is None and \
            args.remove is None and args.keep is None:
        msg = "needs at least one of --exclude, --extract, --remove or --keep"
        raise ProgramError(msg)

    # Check for SNPs
    if args.exclude is not None and args.extract is None:
        if not os.path.isfile(args.exclude):
            msg = "{}: no such file".format(args.exclude)
            raise ProgramError(msg)
    elif args.extract is not None and args.exclude is None:
        if not os.path.isfile(args.extract):
            msg = "{}: no such file".format(args.extract)
            raise ProgramError(msg)
    elif args.exclude is not None and args.extract is not None:
        msg = "use only one of --extract or --exclude"
        raise ProgramError(msg)

    # Check for samples
    if args.remove is not None and args.keep is None:
        if not os.path.isfile(args.remove):
            msg = "{}: no such file".format(args.remove)
            raise ProgramError(msg)
    elif args.keep is not None and args.remove is None:
        if not os.path.isfile(args.keep):
            msg = "{}: no such file".format(args.keep)
            raise ProgramError(msg)
    elif args.remove is not None and args.keep is not None:
        msg = "use only one of --keep or --remove"
        raise ProgramError(msg)

    return True


def parseArgs(argString=None):  # pragma: no cover
    """Parses the command line options and arguments.

    :param argString: the parameters.

    :type argString: list

    :returns: A :py:class:`argparse.Namespace` object created by the
              :py:mod:`argparse` module. It contains the values of the
              different options.

    ============== ====== =====================================================
        Options     Type                      Description
    ============== ====== =====================================================
    ``--ifile``    string The input file prefix.
    ``--is-bfile`` bool   The input file is a bfile
    ``--is-tfile`` bool   The input file is a tfile
    ``--is-file``  bool   The input file is a file
    ``--exclude``  string A file containing SNPs to exclude from the data set.
    ``--extract``  string A file containing SNPs to extract from the data set.
    ``--remove``   string A file containing samples (FID and IID) to remove
                          from the data set.
    ``--keep``     string A file containing samples (FID and IID) to keep from
                          the data set.
    ``--out``      string The prefix of the output files.
    ============== ====== =====================================================

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
pretty_name = "Data subset"
desc = "Subsets genotype data using Plink."
long_desc = "The script excludes a set of markers and/or samples."
parser = argparse.ArgumentParser(description=desc)
parser.add_argument("-v", "--version", action="version",
                    version="pyGenClean version {}".format(__version__))

# The INPUT files
group = parser.add_argument_group("Input File")
group.add_argument("--ifile", type=str, metavar="FILE", required=True,
                   help=("The input file prefix. The format will be specified "
                         "by --is-bfile, --is-tfile or --is-file, for bfile, "
                         "tfile and file, respectively."))
group.add_argument("--is-bfile", action="store_true",
                   help="The file specified by --ifile is a bfile")
group.add_argument("--is-tfile", action="store_true",
                   help="The file specified by --ifile is a tfile")
group.add_argument("--is-file", action="store_true",
                   help="The file specified by --ifile is a file")
# The options
group = parser.add_argument_group("Options")
group.add_argument("--exclude", type=str, metavar="FILE",
                   help=("A file containing SNPs to exclude from the data "
                         "set."))
group.add_argument("--extract", type=str, metavar="FILE",
                   help=("A file containing SNPs to extract from the data "
                         "set."))
group.add_argument("--remove", type=str, metavar="FILE",
                   help=("A file containing samples (FID and IID) to "
                         "remove from the data set."))
group.add_argument("--keep", type=str, metavar="FILE",
                   help=("A file containing samples (FID and IID) to keep "
                         "from the data set."))
# The OUTPUT files
group = parser.add_argument_group("Output File")
group.add_argument("--out", type=str, metavar="FILE", default="subset",
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
