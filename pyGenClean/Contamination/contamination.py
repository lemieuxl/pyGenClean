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
import glob
import pipes
import logging
import argparse
import subprocess

from .. import __version__


logger = logging.getLogger("contamination")


def main(argString=None):
    """The main function of the module.

    :param argString: the options.

    :type argString: list

    These are the steps:

    1. Prints the options.
    2. Compute frequency using Plink.
    3. Runs bafRegress.

    """
    # Getting and checking the options
    args = parseArgs(argString)
    checkArgs(args)

    logger.info("Options used:")
    for key, value in vars(args).iteritems():
        logger.info("  --{} {}".format(key.replace("_", "-"), value))

    # Checks the sample raw data
    logger.info("Checking the raw data files")
    sample_files = check_sample_files(args.bfile + ".fam", args.raw_dir)

    # Finds the markers to extract
    logger.info("Creating extraction list (autosome only)")
    create_extraction_file(args.bfile + ".bim", args.out)

    # Run plink
    logger.info("Computing frequency using Plink")
    run_plink(args.bfile, args.out, args.out + ".to_extract")

    # Run bafRegress
    logger.info("Running bafRegress")
    if args.sge:
        run_bafRegress_sge(sample_files, args.out, args.out + ".to_extract",
                           args.out + ".frq", args)
    else:
        run_bafRegress(sample_files, args.out, args.out + ".to_extract",
                       args.out + ".frq", args)


def check_sample_files(fam_filename, raw_dirname):
    """Checks the raw sample files.

    :param fam_filename: the name of the FAM file.
    :param raw_dirname: the name of the directory containing the raw file.

    :type fam_filename: str
    :type raw_dirname: str

    :returns: the set of all the sample files that are compatible with the FAM
              file.
    :rtype: set

    """
    # Reading the sample identification number from the FAM file
    fam_samples = None
    with open(fam_filename, "r") as i_file:
        fam_samples = {line.split()[1] for line in i_file.read().splitlines()}

    # Checking the files in the raw directory
    sample_files = set()
    all_samples = set()
    for filename in glob.glob(os.path.join(raw_dirname, "*")):
        sample = os.path.splitext(os.path.basename(filename))[0]
        all_samples.add(sample)
        if sample not in fam_samples:
            logger.warning("{}: sample not in FAM file".format(sample))
        else:
            sample_files.add(filename)

    for sample in fam_samples - all_samples:
        logger.warning("{}: sample not in raw directory".format(sample))

    if len(sample_files) == 0:
        raise ProgramError("no sample left for analysis")

    return sample_files


def create_extraction_file(bim_filename, out_prefix):
    """Creates an extraction file (keeping only markers on autosomes).

    :param bim_filename: the name of the BIM file.
    :param out_prefix: the prefix for the output file.

    :type bim_filename: str
    :type out_prefix: str

    """
    o_file = None
    try:
        o_file = open(out_prefix + ".to_extract", "w")
    except IOError:
        raise ProgramError("{}: cannot write file".format(
            out_prefix + ".to_extract"
        ))

    # Reading the BIM file and extracts only the markers on autosome
    autosomes = set(map(str, range(1, 23)))
    nb_markers = 0
    header = dict(zip(["chrom", "name", "cm", "pos", "a1", "a2"], range(6)))
    with open(bim_filename, "r") as i_file:
        for line in i_file:
            row = line.rstrip("\r\n").split()
            if row[header["chrom"]] in autosomes:
                print >>o_file, row[header["name"]]
                nb_markers += 1

    # Closing the file
    o_file.close()

    logger.info("  - {:,d} markers will be used for contamination "
                "estimation".format(nb_markers))


def run_bafRegress(filenames, out_prefix, extract_filename, freq_filename,
                   options):
    """Runs the bafRegress function.

    :param filenames: the set of all sample files.
    :param out_prefix: the output prefix.
    :param extract_filename: the name of the markers to extract.
    :param freq_filename: the name of the file containing the frequency.
    :param options: the other options.

    :type filenames: set
    :type out_prefix: str
    :type extract_filename: str
    :type freq_filename: str
    :type options: argparse.Namespace

    """
    # The command
    command = [
        "bafRegress.py",
        "estimate",
        "--freqfile", freq_filename,
        "--freqcol", "2,5",
        "--extract", extract_filename,
        "--colsample", options.colsample,
        "--colmarker", options.colmarker,
        "--colbaf", options.colbaf,
        "--colab1", options.colab1,
        "--colab2", options.colab2,
    ]
    command.extend(filenames)

    output = None
    try:
        output = subprocess.check_output(command, stderr=subprocess.STDOUT,
                                         shell=False)
    except subprocess.CalledProcessError as exc:
        raise ProgramError("bafRegress.py: couldn't run "
                           "bafRegress.py\n{}".format(exc.output))

    # Saving the output
    try:
        with open(out_prefix + ".bafRegress", "w") as o_file:
            o_file.write(output)
    except IOError:
        raise ProgramError("{}: cannot write file".format(
            out_prefix + ".bafRegress",
        ))


def run_bafRegress_sge(filenames, out_prefix, extract_filename, freq_filename,
                       options):
    """Runs the bafRegress function using SGE.

    :param filenames: the set of all sample files.
    :param out_prefix: the output prefix.
    :param extract_filename: the name of the markers to extract.
    :param freq_filename: the name of the file containing the frequency.
    :param options: the other options.

    :type filenames: set
    :type out_prefix: str
    :type extract_filename: str
    :type freq_filename: str
    :type options: argparse.Namespace

    """
    # Checks the environment variable for DRMAA package
    if "DRMAA_LIBRARY_PATH" not in os.environ:
        msg = "could not load drmaa: set DRMAA_LIBRARY_PATH"
        raise ProgramError(msg)

    # Import the drmaa package
    try:
        import drmaa
    except ImportError:
        raise ProgramError("drmaa is not install, install drmaa")

    # Initializing a session
    s = drmaa.Session()
    s.initialize()

    # The base command
    base_command = [
        "bafRegress.py",
        "estimate",
        "--freqfile", freq_filename,
        "--freqcol", "2,5",
        "--extract", extract_filename,
        "--colsample", options.colsample,
        "--colmarker", options.colmarker,
        "--colbaf", options.colbaf,
        "--colab1", options.colab1,
        "--colab2", options.colab2,
    ]

    # Creating chunks
    chunks = [
        list(filenames)[i: i+options.sample_per_run_for_sge]
        for i in range(0, len(filenames), options.sample_per_run_for_sge)
    ]

    # Run for each sub task...
    job_ids = []
    job_templates = []
    for i, chunk in enumerate(chunks):
        # Creating the final command
        command = [pipes.quote(token) for token in base_command + chunk]
        command.append(
            "> {}.bafRegress_{}".format(pipes.quote(out_prefix), i+1),
        )
        command = " ".join(command)

        # Creating the job template
        jt = s.createJobTemplate()
        jt.remoteCommand = "bash"
        jt.workingDirectory = os.getcwd()
        jt.jobEnvironment = os.environ
        jt.args = ["-c", command]
        jt.jobName = "_bafRegress_{}".format(i+1)

        # Cluster specifics
        if options.sge_walltime is not None:
            jt.hardWallclockTimeLimit = options.sge_walltime
        if options.sge_nodes is not None:
            native_spec = "-l nodes={}:ppn={}".format(options.sge_nodes[0],
                                                      options.sge_nodes[1])
            jt.nativeSpecification = native_spec

        job_ids.append(s.runJob(jt))
        job_templates.append(jt)

    # Waiting for the jobs to finish
    had_problems = []
    for job_id in job_ids:
        ret_val = s.wait(job_id, drmaa.Session.TIMEOUT_WAIT_FOREVER)
        had_problems.append(ret_val.exitStatus == 0)

    # Deleting the jobs
    for jt in job_templates:
        s.deleteJobTemplate(jt)

    # Closing the session
    s.exit()

    # Checking for problems
    if not all(had_problems):
        raise ProgramError("Some SGE jobs had errors...")

    # Merging the output files
    o_file = None
    try:
        o_file = open(out_prefix + ".bafRegress", "w")
    except IOError:
        msg = "{}: cannot write file".format(out_prefix + ".bafRegress")
        raise ProgramError(msg)

    to_remove = set()
    for i in range(len(chunks)):
        filename = out_prefix + ".bafRegress_{}".format(i+1)
        if not os.path.isfile(filename):
            raise ProgramError("{}: no such file".format(filename))
        to_remove.add(filename)

        with open(filename, "r") as i_file:
            if i == 0:
                # First file, we write everything
                o_file.write(i_file.read())
                continue

            for j, line in enumerate(i_file):
                # Skipping first line
                if j == 0:
                    continue
                o_file.write(line)

    for filename in to_remove:
        # Deleting the chunk
        os.remove(filename)

    # Closing
    o_file.close()


def run_plink(in_prefix, out_prefix, extract_filename):
    """Runs Plink with the geno option.

    :param in_prefix: the input prefix.
    :param out_prefix: the output prefix.
    :param extract_filename: the name of the file containing markers to
                             extract.

    :type in_prefix: str
    :type out_prefix: str
    :param extract_filename: str

    """
    # The plink command
    plink_command = [
        "plink",
        "--noweb",
        "--bfile", in_prefix,
        "--extract", extract_filename,
        "--freq",
        "--out", out_prefix,
    ]
    output = None
    try:
        output = subprocess.check_output(plink_command,
                                         stderr=subprocess.STDOUT, shell=False)
    except subprocess.CalledProcessError as exc:
        msg = "plink: couldn't run plink\n{}".format(exc.output)
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
    for filename in [args.bfile + i for i in [".bed", ".bim", ".fam"]]:
        if not os.path.isfile(filename):
            raise ProgramError("{}: no such file".format(filename))

    # Checking that the raw directory exists
    if not os.path.isdir(args.raw_dir):
        raise ProgramError("{}: no such directory".format(args.raw_dir))

    return True


def parseArgs(argString=None):  # pragma: no cover
    """Parses the command line options and arguments.

    :param argString: the options.

    :type argString: list

    :returns: A :py:class:`argparse.Namespace` object created by the
              :py:mod:`argparse` module. It contains the values of the
              different options.

    ============================ ====== =======================================
               Options            Type                Description
    ============================ ====== =======================================
    ``--bfile``                  string The input file prefix (will find the
                                        plink binary files by appending the
                                        prefix to the .bim, .bed and .fam
                                        files, respectively).
    ``--raw-dir``                string Directory containing the raw data (one
                                        file per sample, where the name of the
                                        file (minus the extension) is the
                                        sample identification number.
    ``--colsample``              string The sample column.
    ``--colmarker``              string The marker column.
    ``--colbaf``                 string The B allele frequency column.
    ``--colab1``                 string The AB Allele 1 column.
    ``--colab2``                 string The AB Allele 2 column.
    ``--out``                    string The prefix of the output files.
    ``--sge``                    bool   Use SGE for parallelization.
    ``--sge-walltime``           string The walltime for the job to run on the
                                        cluster. Do not use if you are not
                                        required to specify a walltime for your
                                        jobs on your cluster (e.g. 'qsub
                                        -lwalltime=1:0:0' on the cluster).
    ``--sge-nodes``              int    The number of nodes and the number of
                                        processor per nodes to use (e.g. 'qsub
                                        -lnodes=X:ppn=Y' on the cluster, where
                                        X is the number of nodes and Y is  the
                                        number of processor to use. Do not use
                                        if you are not required to specify the
                                        number of nodes for your jobs on the
                                        cluster.
    ``--sample-per-run-for-sge`` int    The number of sample to run for a
                                        single SGE job.
    ============================ ====== =======================================

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
pretty_name = "Contamination"
desc = "Check BAF and LogR ratio for data contamination."
long_desc = ("The script search for sample contamination using the "
             r"\path{bafRegress.py} software.")
parser = argparse.ArgumentParser(description=desc)
parser.add_argument("-v", "--version", action="version",
                    version="pyGenClean version {}".format(__version__))

# The INPUT files
group = parser.add_argument_group("Input File")
group.add_argument("--bfile", type=str, metavar="FILE", required=True,
                   help=("The input file prefix (will find the plink binary "
                         "files by appending the prefix to the .bim, .bed and "
                         ".fam files, respectively)."))

group = parser.add_argument_group("Raw Data")
group.add_argument("--raw-dir", type=str, metavar="DIR", required=True,
                   help=("Directory containing the raw data (one file per "
                         "sample, where the name of the file (minus the "
                         "extension) is the sample identification number."))

group.add_argument("--colsample", type=str, metavar="COL",
                   default="Sample Name",
                   help="The sample column. [default: %(default)s]")
group.add_argument("--colmarker", type=str, metavar="COL",
                   default="SNP Name",
                   help="The marker column. [default: %(default)s]")
group.add_argument("--colbaf", type=str, metavar="COL",
                   default="B Allele Freq",
                   help=("The B allele frequency column. "
                         "[default: %(default)s]"))
group.add_argument("--colab1", type=str, metavar="COL",
                   default="Allele1 - AB",
                   help="The AB Allele 1 column. [default: %(default)s]")
group.add_argument("--colab2", type=str, metavar="COL",
                   default="Allele2 - AB",
                   help="The AB Allele 2 column. [default: %(default)s]")

group = parser.add_argument_group("SGE Options")
group.add_argument("--sge", action="store_true",
                   help="Use SGE for parallelization.")
group.add_argument("--sge-walltime", type=str, metavar="TIME",
                   help=("The walltime for the job to run on the cluster. Do "
                         "not use if you are not required to specify a "
                         "walltime for your jobs on your cluster (e.g. 'qsub "
                         "-lwalltime=1:0:0' on the cluster)."))
group.add_argument("--sge-nodes", type=int, metavar="INT", nargs=2,
                   help=("The number of nodes and the number of processor per "
                         "nodes to use (e.g. 'qsub -lnodes=X:ppn=Y' on the "
                         "cluster, where X is the number of nodes and Y is "
                         " the number of processor to use. Do not use if you "
                         "are not required to specify the number of nodes for "
                         "your jobs on the cluster."))
group.add_argument("--sample-per-run-for-sge", type=int, metavar="INT",
                   default=30,
                   help=("The number of sample to run for a single SGE job. "
                         "[default: %(default)d]"))

# The OUTPUT files
group = parser.add_argument_group("Output File")
group.add_argument("--out", type=str, metavar="FILE",
                   default="contamination",
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
