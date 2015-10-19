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
import re
import sys
import time
import shutil
import logging
import argparse
from glob import glob

from . import auto_report
from ..PlinkUtils import get_plink_version


logger = logging.getLogger("merge_reports")


def main(argString=None):
    """The main function of this module.

    :param argString: the options.

    :type argString: list of strings

    """
    # Getting and checking the options
    args = parseArgs(argString)
    checkArgs(args)

    logger.info("Options used:")
    for key, value in vars(args).iteritems():
        logger.info("  --{} {}".format(key.replace("_", "-"), value))

    # Checking if the output directory exists, creating it otherwise
    if not os.path.isdir(args.out_dir):
        os.mkdir(args.out_dir)

    # Ordering the directories according to their name
    qc_dir = order_qc_dir(args.qc_dir)

    # First, we want to merge the required files
    merge_required_files(qc_dir, args.out_dir)

    # Then, we want to copy the initial_files file
    copy_initial_files(os.path.join(qc_dir[0], "initial_files.txt"),
                       args.out_dir)

    # Get the final number of markers and samples
    final_nb_markers, final_nb_samples = get_final_numbers(
        os.path.join(qc_dir[-1], "final_files.txt"),
        args.out_dir,
    )

    # Getting the steps summary file (TeX)
    summary_files = get_summary_files(qc_dir)

    # Generating the report
    generate_report(args.out_dir, summary_files, final_nb_markers,
                    final_nb_samples, args)


def order_qc_dir(dirnames):
    """Order the QC directory names according to their date.

    :param dirnames: the list of directories to merge data from.

    :type dirnames: list

    :returns: the sorted list of directories
    :rtype: list

    """
    return sorted(
        dirnames, key=lambda dn: time.strptime(
            os.path.basename(dn.rstrip("/"))[14:],
            "%Y-%m-%d_%H.%M.%S",
        )
    )


def merge_required_files(dirnames, out_dir):
    """Merges the required files from each of the directories.

    :param dirnames: the list of directories to merge data from.
    :param out_dir: the name of the output directory.

    :type dirnames: list
    :type out_dir: str

    """
    # The list of files to merge
    fn_to_merge = ("steps_summary.tex", "excluded_markers.txt",
                   "excluded_samples.txt")

    # Merging the files
    for fn in fn_to_merge:
        o_fn = os.path.join(out_dir, fn)
        with open(o_fn, "w") as o_file:
            for dn in dirnames:
                i_fn = os.path.join(dn, fn)
                with open(i_fn, "r") as i_file:
                    o_file.write(i_file.read())

    # Merging the result summary file
    o_fn = os.path.join(out_dir, "results_summary.txt")
    with open(o_fn, "w") as o_file:
        for i, dn in enumerate(dirnames):
            i_fn = os.path.join(dn, "results_summary.txt")
            with open(i_fn, "r") as i_file:
                if i != 0:
                    # We skip the first 4 lines (file descriptions)
                    [i_file.readline() for i in range(4)]
                o_file.write(i_file.read())

    # Merging the graphic paths file
    graphic_paths = set()
    for dn in dirnames:
        fn = os.path.join(dn, "graphic_paths.txt")
        if os.path.isfile(fn):
            with open(fn, "r") as i_file:
                graphic_paths.update({
                    os.path.join(dn, path)
                    for path in i_file.read().splitlines()
                })
    if len(graphic_paths) > 0:
        with open(os.path.join(out_dir, "graphic_paths.txt"), "w") as o_file:
            for path in sorted(graphic_paths):
                print >>o_file, os.path.relpath(path, out_dir)


def copy_initial_files(filename, out_dir):
    """Copy the initial_files file to the final directory.

    :param filename: the name of the file.
    :param out_dir: the name of the output directory

    :type dirname: str
    :type out_dir: str

    """
    shutil.copy(filename, out_dir)


def get_final_numbers(filename, out_dir):
    """Copy the final_files file and get the number of markers and samples.

    :param filename: the name of the file.
    :param out_dir: the output directory.

    :type filename: str
    :type out_dir: str

    :returns: the final number of markers and samples
    :rtype: tuple

    """
    # Copying the file
    shutil.copy(filename, out_dir)

    # Reading the number of markers and samples
    nb_samples = None
    nb_markers = None
    with open(filename, "r") as i_file:
        for line in i_file:
            row = line.rstrip("\r\n").split("\t")
            if len(row) == 1:
                continue
            path, ext = os.path.splitext(row[0])
            if ext in {".bim", ".tped", ".map"}:
                nb_markers = row[1]
            elif ext in {".fam", ".ped", ".tfam"}:
                nb_samples = row[1]

    assert nb_samples
    assert nb_markers

    return nb_markers, nb_samples


def get_summary_files(dirnames):
    """Gets the TeX summary files for each test.

    :param dirnames: the list of directories to merge data from.

    :type dirnames: list

    :returns: a list of summary file names.
    :rtype: list

    """
    # A useful regular expression to get step number in the current directory
    step_nb_re = re.compile(r"^([0-9]+)_\S+")

    # The final list of summary files
    final_summary_files = []

    # For each of the directory
    for dn in dirnames:
        # Getting the step directories
        step_dir = [
            n for n in os.listdir(dn)
            if os.path.isdir(os.path.join(dn, n)) and step_nb_re.match(n)
        ]

        # Sorting the step directories
        step_dir.sort(key=lambda x: int(step_nb_re.match(x).group(1)))

        # Getting the name of the summary file for each of the step directory
        step_summary_files = [
            glob(os.path.join(dn, sn, "*.summary.tex")) for sn in step_dir
        ]

        # Checking we have only one summary file
        for summary_file in step_summary_files:
            if len(summary_file) > 1:
                raise ProgramError("{}: multiple summary files".format(
                    os.path.join(dn, sn),
                ))

            if not summary_file:
                raise ProgramError("{}: missing summary file".format(
                    os.apth.join(dn, sn),
                ))

        final_summary_files.extend(i[0] for i in step_summary_files)

    return [os.path.abspath(fn) for fn in final_summary_files]


def generate_report(out_dir, latex_summaries, nb_markers, nb_samples, options):
    """Generates the report.

    :param out_dir: the output directory.
    :param latex_summaries: the list of LaTeX summaries.
    :param nb_markers: the final number of markers.
    :param nb_samples: the final number of samples.
    :param options: the list of options.

    :type out_dir: str
    :type latex_summaries: list
    :type nb_markers: str
    :type nb_samples: str
    :type options: argparse.Namespace

    """
    # Getting the graphic paths file
    graphic_paths_fn = None
    if os.path.isfile(os.path.join(out_dir, "graphic_paths.txt")):
        graphic_paths_fn = os.path.join(out_dir, "graphic_paths.txt")

    # We create the automatic report
    report_name = os.path.join(out_dir, "merged_report.tex")
    auto_report.create_report(
        out_dir,
        report_name,
        project_name=options.report_number,
        steps_filename=os.path.join(out_dir, "steps_summary.tex"),
        summaries=latex_summaries,
        background=options.report_background,
        summary_fn=os.path.join(out_dir, "results_summary.txt"),
        report_title=options.report_title,
        report_author=options.report_author,
        initial_files=os.path.join(out_dir, "initial_files.txt"),
        final_files=os.path.join(out_dir, "final_files.txt"),
        final_nb_markers=nb_markers,
        final_nb_samples=nb_samples,
        plink_version=get_plink_version(),
        graphic_paths_fn=graphic_paths_fn,
    )


def checkArgs(args):
    """Checks the arguments and options.

    :param args: an object containing the options of the program.

    :type args: argparse.Namespace

    :returns: ``True`` if everything was OK.

    If there is a problem with an option, an exception is raised using the
    :py:class:`ProgramError` class, a message is printed to the
    :class:`sys.stderr` and the program exists with code 1.

    """
    # For all input directories
    for dn in args.qc_dir:
        # Checking that all the directories exists
        if not os.path.isdir(dn):
            raise ProgramError("{}: no such directory".format(dn))

        # Checking that this is a directory created from pyGenClean
        if not os.path.basename(dn.rstrip("/")).startswith("data_clean_up."):
            raise ProgramError("{}: not a pyGenClean directory".format(dn))

        # Checking that each directory contains the required files
        for fn in ("excluded_markers.txt", "excluded_samples.txt",
                   "results_summary.txt", "steps_summary.tex",
                   "initial_files.txt", "final_files.txt"):
            required_fn = os.path.join(dn, fn)
            if not os.path.isfile(required_fn):
                raise ProgramError("{}: missing required "
                                   "file".format(required_fn))

    return True


def parseArgs(argString=None):  # pragma: no cover
    """Parses the command line options and arguments.

    :param argString: the options.

    :type argString: list of strings

    :returns: A :py:class:`argparse.Namespace` object created by the
              :py:mod:`argparse` module. It contains the values of the
              different options.

    ======================= ====== ============================================
             Options         Type                    Description
    ======================= ====== ============================================
    ``--report-author``     String The current project number.
    ``--report-number``     String The current project author.
    ``--report-background`` String Text of file containing the background
                                   section of the report.
    ======================= ====== ============================================

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


def add_custom_options(parser):
    """Adds custom options to a parser.

    :param parser: the parser to which to add options.

    :type parser: argparse.ArgumentParser

    """
    parser.add_argument("--report-title", type=str, metavar="TITLE",
                        default="Genetic Data Clean Up",
                        help="The report title. [default: %(default)s]")
    parser.add_argument("--report-author", type=str, metavar="AUTHOR",
                        default="pyGenClean",
                        help="The current project number. "
                             "[default: %(default)s]")
    parser.add_argument("--report-number", type=str, metavar="NUMBER",
                        default="Simple Project",
                        help="The current project author. "
                             "[default: %(default)s]")
    parser.add_argument("--report-background", type=str, metavar="BACKGROUND",
                        default="The aim of this project is to perform data "
                                "QC prior to genetic analysis.",
                        help="Text of file containing the background section "
                             "of the report.")


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
pretty_name = "Merge reports"
desc = """Merges automatic reports from other pyGenClean runs."""
parser = argparse.ArgumentParser(description=desc)

# The INPUT files
group = parser.add_argument_group("Input")
group.add_argument("--qc-dir", nargs="+", required=True, metavar="DIR",
                   help="A list of directory containing pyGenClean runs.")

# The options
group = parser.add_argument_group("Report Options")
add_custom_options(group)

# The OUTPUT files
group = parser.add_argument_group("Output Directory")
group.add_argument("--out-dir", type=str, metavar="FILE",
                   default="pyGenClean_report",
                   help=("The name of the directory that will contain the "
                         "final report. [default: %(default)s]"))


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
