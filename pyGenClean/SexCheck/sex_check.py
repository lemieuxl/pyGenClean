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

import numpy as np

from .import gender_plot
from . import baf_lrr_plot

from .. import __version__
from ..PlinkUtils import createRowFromPlinkSpacedOutput


logger = logging.getLogger("sex_check")


def main(argString=None):
    """The main function of the module.

    :param argString: the options.

    :type argString: list

    These are the following steps:

    1.  Prints the options.
    2.  Checks if there are enough markers on the chromosome ``23``
        (:py:func:`checkBim`). If not, quits here.
    3.  Runs the sex check analysis using Plink (:py:func:`runPlinkSexCheck`).
    4.  If there are no sex problems, then quits (:py:func:`readCheckSexFile`).
    5.  Creates the recoded file for the chromosome ``23``
        (:py:func:`createPedChr23UsingPlink`).
    6.  Computes the heterozygosity percentage on the chromosome ``23``
        (:py:func:`computeHeteroPercentage`).
    7.  If there are enough markers on chromosome ``24`` (at least 1), creates
        the recoded file for this chromosome
        (:py:func:`createPedChr24UsingPlink`).
    8.  Computes the number of no call on the chromosome ``24``
        (:py:func:`computeNoCall`).
    9.  If required, plots the gender plot (:py:func:`createGenderPlot`).
    10. If required, plots the BAF and LRR plot (:py:func:`createLrrBafPlot`).

    """
    # Getting and checking the options
    args = parseArgs(argString)
    checkArgs(args)

    logger.info("Options used:")
    for key, value in vars(args).iteritems():
        logger.info("  --{} {}".format(key.replace("_", "-"), value))

    # Reads the bim file to see if chromosome 23 is there
    hasSexProblems = None
    if not checkBim("{}.bim".format(args.bfile), args.nbChr23, "23"):
        logger.info("There are not enough markers on chromosome 23: "
                    "STOPPING NOW!")
    else:
        # Run plink "check-sex"
        logger.info("Running Plink for sex check")
        runPlinkSexCheck(args)

        # Reading plink "check-sex" output file
        logger.info("Reading Plink's sex check output to find sex problems")
        hasSexProblems = readCheckSexFile(args.out + ".sexcheck",
                                          args.out + ".list_problem_sex",
                                          args.out + ".list_problem_sex_ids",
                                          args.femaleF, args.maleF)

    if hasSexProblems is not None and hasSexProblems:
        # Run plink to recode chr 23 in a ped format
        logger.info("Creating recoded file for chr23 using Plink")
        createPedChr23UsingPlink(args)

        # Compute the hetero percentage
        logger.info("Computing the heterozygous percentage")
        computeHeteroPercentage(args.out + ".chr23_recodeA.raw")

        # Run plink to get chr 24
        if checkBim("{}.bim".format(args.bfile), 1, "24"):
            logger.info("Creating recoded file for chr24 using Plink")
            createPedChr24UsingPlink(args)

            # Compute the number of no call
            logger.info("Computing the number of no calls")
            computeNoCall(args.out + ".chr24_recodeA.raw")
        else:
            logger.info("Not enough markers on chr24")

        # If required, let's plot the gender plot
        if args.gender_plot:
            logger.info("Creating the gender plot")
            createGenderPlot(args.bfile, args.sex_chr_intensities,
                             args.out + ".list_problem_sex",
                             args.gender_plot_format, args.out)

        # If required, let's plot the LRR and BAF plot
        if args.lrr_baf:
            logger.info("Creating the LRR and BAF plot")
            createLrrBafPlot(
                raw_dir=args.lrr_baf_raw_dir,
                problematic_samples=args.out + ".list_problem_sex_ids",
                format=args.lrr_baf_format,
                dpi=args.lrr_baf_dpi,
                out_prefix=args.out,
            )


def createGenderPlot(bfile, intensities, problematic_samples, format,
                     out_prefix):
    """Creates the gender plot.

    :param bfile: the prefix of the input binary file.
    :param intensities: the file containing the intensities.
    :param problematic_samples: the file containing the problematic samples.
    :param format: the format of the output plot.
    :param out_prefix: the prefix of the output file.

    :type bfile: str
    :type intensities: str
    :type problematic_samples: str
    :type format: str
    :type out_prefix: str

    Creates the gender plot of the samples using the
    :py:mod:`pyGenClean.SexCheck.gender_plot` module.

    """
    gender_plot_options = ["--bfile", bfile, "--intensities", intensities,
                           "--sex-problems", problematic_samples, "--format",
                           format, "--out", out_prefix]
    try:
        gender_plot.main(gender_plot_options)
    except gender_plot.ProgramError as e:
        msg = "gender plot: {}".format(e)
        raise ProgramError(msg)


def createLrrBafPlot(raw_dir, problematic_samples, format, dpi, out_prefix):
    """Creates the LRR and BAF plot.

    :param raw_dir: the directory containing the intensities.
    :param problematic_samples: the file containing the problematic samples.
    :param format: the format of the plot.
    :param dpi: the DPI of the resulting images.
    :param out_prefix: the prefix of the output file.

    :type raw_dir: str
    :type problematic_samples: str
    :type format: str
    :type out_prefix: str

    Creates the LRR (Log R Ratio) and BAF (B Allele Frequency) of the
    problematic samples using the :py:mod:`pyGenClean.SexCheck.baf_lrr_plot`
    module.

    """
    # First, we create an output directory
    dir_name = out_prefix + ".LRR_BAF"
    if not os.path.isdir(dir_name):
        os.mkdir(dir_name)

    # The options
    baf_lrr_plot_options = ["--problematic-samples", problematic_samples,
                            "--raw-dir", raw_dir, "--format", format,
                            "--dpi", str(dpi),
                            "--out", os.path.join(dir_name, "baf_lrr")]
    try:
        baf_lrr_plot.main(baf_lrr_plot_options)
    except baf_lrr_plot.ProgramError as e:
        msg = "BAF LRR plot: {}".format(e)
        raise ProgramError(msg)


def checkBim(fileName, minNumber, chromosome):
    """Checks the BIM file for chrN markers.

    :param fileName:
    :param minNumber:
    :param chromosome:

    :type fileName: str
    :type minNumber: int
    :type chromosome: str

    :returns: ``True`` if there are at least ``minNumber`` markers on
              chromosome ``chromosome``, ``False`` otherwise.

    """
    nbMarkers = 0
    with open(fileName, 'r') as inputFile:
        for line in inputFile:
            row = line.rstrip("\r\n").split("\t")
            if row[0] == chromosome:
                nbMarkers += 1
    if nbMarkers < minNumber:
        return False
    return True


def computeNoCall(fileName):
    """Computes the number of no call.

    :param fileName: the name of the file

    :type fileName: str

    Reads the ``ped`` file created by Plink using the ``recodeA`` options (see
    :py:func:`createPedChr24UsingPlink`) and computes the number and percentage
    of no calls on the chromosome ``24``.
    """
    outputFile = None
    try:
        outputFile = open(fileName + ".noCall", "w")
    except IOError:
        msg = "%s: can't write file" % fileName + ".noCall"
        raise ProgramError(msg)
    print >>outputFile, "\t".join(["PED", "ID", "SEX", "nbGeno", "nbNoCall"])

    try:
        toPrint = []
        with open(fileName, "r") as inputFile:
            for i, line in enumerate(inputFile):
                row = line.rstrip("\r\n").split(" ")
                if i != 0:
                    # This is data
                    genotypes = np.array(row[6:])

                    nbMarker = len(genotypes)
                    nbNA = len(np.where(genotypes == "NA")[0])

                    toPrint.append((row[0], row[1], row[4], str(nbMarker),
                                    str(nbNA)))

        toPrint.sort(reverse=True, key=lambda values: int(values[4]))

        for row in toPrint:
            print >>outputFile, "\t".join(row)

    except IOError:
        msg = "%(fileName)s: no such file" % locals()
        raise ProgramError(msg)

    # Closing the output file
    outputFile.close()


def computeHeteroPercentage(fileName):
    """Computes the heterozygosity percentage.

    :param fileName: the name of the input file.

    :type fileName: str

    Reads the ``ped`` file created by Plink using the ``recodeA`` options (see
    :py:func:`createPedChr23UsingPlink`) and computes the heterozygosity
    percentage on the chromosome ``23``.

    """
    outputFile = None
    try:
        outputFile = open(fileName + ".hetero", "w")
    except IOError:
        msg = "%s: can't write file" % fileName + ".hetero"
        raise ProgramError(msg)
    print >>outputFile, "\t".join(["PED", "ID", "SEX", "HETERO"])

    try:
        toPrint = []
        with open(fileName, "r") as inputFile:
            for i, line in enumerate(inputFile):
                row = line.rstrip("\r\n").split(" ")
                if i != 0:
                    # This is data
                    genotypes = np.array(row[6:])
                    genotypes = genotypes[np.where(genotypes != "NA")]

                    nbMarker = len(genotypes)
                    nbHetero = len(np.where(genotypes == "1")[0])
                    percentHetero = -9999
                    if nbMarker != 0:
                        percentHetero = nbHetero / float(nbMarker)

                    toPrint.append((row[0], row[1], row[4], percentHetero))

        # Sorting the data
        toPrint.sort(reverse=True, key=lambda values: values[3])

        # Printing the data
        for row in toPrint:
            value = list(row)
            if value[3] == -9999:
                value[3] = "ALL_NA"
            else:
                value[3] = str(value[3])
            print >>outputFile, "\t".join(value)

    except IOError:
        msg = "%(fileName)s: no such file" % locals()
        raise ProgramError(msg)

    # Closing the output file
    outputFile.close()


def readCheckSexFile(fileName, allProblemsFileName, idsFileName, femaleF,
                     maleF):
    """Reads the Plink check-sex output file.

    :param fileName: the name of the input file.
    :param allProblemsFileName: the name of the output file that will contain
                                all the problems.
    :param idsFileName: the name of the output file what will contain samples
                        with sex problems.
    :param femaleF: the F threshold for females.
    :param maleF: the F threshold for males.

    :type fileName: str
    :type allProblemsFileName: str
    :type idsFileName: str
    :type femaleF: float
    :type maleF: float

    :returns: ``True`` if there are sex problems, ``False`` otherwise.

    Reads sex check file provided by :py:func:`runPlinkSexCheck` (Plink) and
    extract the samples that have sex problems.

    """
    allProblemsFile = None
    try:
        allProblemsFile = open(allProblemsFileName, 'w')
    except IOError:
        msg = "%(allProblemsFileName)s: can't write file" % locals()
        raise ProgramError(msg)

    idsFile = None
    try:
        idsFile = open(idsFileName, 'w')
    except IOError:
        msg = "%(idsFileName)s: can't write file" % locals()
        raise ProgramError(msg)

    try:
        with open(fileName, 'r') as inputFile:
            headerIndex = None
            nbProblems = 0
            nbTotalProblems = 0
            nbSexUnknown = 0
            nbFemaleThreshold = 0
            nbMaleThreshold = 0
            for line in inputFile:
                row = createRowFromPlinkSpacedOutput(line)

                if headerIndex is None:
                    # This is the header
                    headerIndex = dict([
                        (row[i], i) for i in xrange(len(row))
                    ])
                    for columnName in ["STATUS", "PEDSEX", "SNPSEX", "F",
                                       "FID", "IID"]:
                        if columnName not in headerIndex:
                            msg = "%(fileName)s: no column named " \
                                  "%(columnName)s" % locals()
                            raise ProgramError(msg)
                    print >>allProblemsFile, "\t".join(row)
                    continue

                # We have data
                status = row[headerIndex["STATUS"]]

                if status == "PROBLEM":
                    # We have a sex problem
                    nbTotalProblems += 1
                    pedsex = row[headerIndex["PEDSEX"]]

                    if pedsex == "0":
                        # The individual was "0", so we skip
                        nbSexUnknown += 1
                        continue

                    snpsex = row[headerIndex["SNPSEX"]]
                    if snpsex == "0":
                        # The new sex is unknown
                        f = None
                        try:
                            f = float(row[headerIndex["F"]])
                        except ValueError:
                            msg = "F=%s: not a float" % row[headerIndex["F"]]
                            raise ProgramError(msg)

                        if pedsex == "2":
                            # We have a female
                            if f < femaleF:
                                nbFemaleThreshold += 1
                                continue

                        if pedsex == "1":
                            # We have a male
                            if f > maleF:
                                nbMaleThreshold += 1
                                continue

                    print >>allProblemsFile, "\t".join(row)

                    famID = row[headerIndex["FID"]]
                    indID = row[headerIndex["IID"]]
                    print >>idsFile, "\t".join([famID, indID])

                    nbProblems += 1

            logger.info("Sex Check Summary")
            logger.info("  - {:,d} total problems".format(nbTotalProblems))
            logger.info("  - {:,d} pedsex unknown".format(nbSexUnknown))
            logger.info("  - {:,d} female F < {}".format(nbFemaleThreshold,
                                                         femaleF))
            logger.info("  - {:,d} male F > {}".format(nbMaleThreshold, maleF))
            logger.info("  - {:,d} problems kept".format(nbProblems))

    except IOError:
        msg = "%(fileName)s: no such file"

    # Closing the output files
    idsFile.close()
    allProblemsFile.close()

    if nbProblems == 0:
        # There are no sex problems to investigate
        logger.info("There are no sex problem to investigate...")
        logger.info("  - Nothing else to do...")
        return False
    return True


def runPlinkSexCheck(options):
    """Runs Plink to perform a sex check analysis.

    :param options: the options.

    :type options: argparse.Namespace

    Uses Plink to perform a sex check analysis.

    """
    # The plink command
    plinkCommand = ["plink", "--noweb", "--bfile", options.bfile,
                    "--check-sex", "--out", options.out]
    runCommand(plinkCommand)


def createPedChr23UsingPlink(options):
    """Run Plink to create a ped format.

    :param options: the options.

    :type options: argparse.Namespace

    Uses Plink to create a ``ped`` file of markers on the chromosome ``23``. It
    uses the ``recodeA`` options to use additive coding. It also subsets the
    data to keep only samples with sex problems.

    """
    plinkCommand = ["plink", "--noweb", "--bfile", options.bfile, "--chr",
                    "23", "--recodeA", "--keep",
                    options.out + ".list_problem_sex_ids", "--out",
                    options.out + ".chr23_recodeA"]
    runCommand(plinkCommand)


def createPedChr24UsingPlink(options):
    """Run plink to create a ped format.

    :param options: the options.

    :type options: argparse.Namespace

    Uses Plink to create a ``ped`` file of markers on the chromosome ``24``. It
    uses the ``recodeA`` options to use additive coding. It also subsets the
    data to keep only samples with sex problems.

    """
    plinkCommand = ["plink", "--noweb", "--bfile", options.bfile, "--chr",
                    "24", "--recodeA", "--keep",
                    options.out + ".list_problem_sex_ids", "--out",
                    options.out + ".chr24_recodeA"]
    runCommand(plinkCommand)


def runCommand(command):
    """Run a command.

    :param command: the command to run.

    :type command: list

    Tries to run a command. If it fails, raise a :py:class:`ProgramError`. This
    function uses the :py:mod:`subprocess` module.

    .. warning::
        The variable command should be a list of strings (no other type).

    """
    output = None
    try:
        output = subprocess.check_output(command,
                                         stderr=subprocess.STDOUT, shell=False)
    except subprocess.CalledProcessError:
        msg = "plink: couldn't run plink\n"
        msg += " ".join(command)
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

    # Ceck the number of markers on chromosome 23
    if args.nbChr23 < 0:
        msg = ("{}: number of markers on chr 23 must be "
               "positive".format(args.nbChr23))
        raise ProgramError(msg)

    # If we ask for LRR and BAF, we need a directory
    if args.lrr_baf:
        if not os.path.isdir(args.lrr_baf_raw_dir):
            msg = "{}: no such directory".format(args.lrr_baf_raw_dir)
            raise ProgramError(msg)
        if args.lrr_baf_dpi < 10:
            msg = "{}: DPI too low".format(args.dpi)
            raise ProgramError(msg)

    return True


def parseArgs(argString=None):  # pragma: no cover
    """Parses the command line options and arguments.

    :param argString: the options.

    :type argString: list

    :returns: A :py:class:`argparse.Namespace` object created by the
              :py:mod:`argparse` module. It contains the values of the
              different options.

    ========================= ====== ==========================================
             Options           Type                  Description
    ========================= ====== ==========================================
    ``--bfile``               string The input file prefix (Plink binary).
    ``--femaleF``             float  The female F threshold.
    ``--maleF``               float  The male F threshold.
    ``--nbChr23``             int    The minimum number of markers on
                                     chromosome 23 before computing Plink's sex
                                     check.
    ``--gender-plot``         bool   Create the gender plot.
    ``--sex-chr-intensities`` string A file containing alleles intensities for
                                     each of the markers located on the X and Y
                                     chromosome.
    ``--gender-plot-format``  string The output file format for the gender
                                     plot.
    ``--lrr-baf``             bool   Create the LRR and BAF plot.
    ``--lrr-baf-raw-dir``     string Directory containing information about
                                     every samples (BAF and LRR).
    ``--lrr-baf-format``      string The output file format.
    ``--lrr-baf-dpi``         int    The pixel density of the figure(s) (DPI).
    ``--out``                 string The prefix of the output files.
    ========================= ====== ==========================================

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

        :param msg: the message to print to the user.

        :type msg: str

        """
        self.message = str(msg)

    def __str__(self):
        return self.message


# The parser object
pretty_name = "Gender check"
desc = "Check sample's gender using Plink."
long_desc = ("The script identifies any individual with discrepancies between "
             "phenotype and genotype data for sex. Individuals with sex error "
             "are to be investigated.")
parser = argparse.ArgumentParser(description=desc)
parser.add_argument("-v", "--version", action="version",
                    version="pyGenClean version {}".format(__version__))

# The INPUT files
group = parser.add_argument_group("Input File")
group.add_argument("--bfile", type=str, metavar="FILE", required=True,
                   help=("The input file prefix (will find the Plink binary "
                         "files by appending the prefix to the .bed, .bim, "
                         "and .fam files, respectively."))
# The options
group = parser.add_argument_group("Options")
group.add_argument("--femaleF", type=float, metavar="FLOAT", default=0.3,
                   help="The female F threshold. [default: < %(default)f]")
group.add_argument("--maleF", type=float, metavar="FLOAT", default=0.7,
                   help="The male F threshold. [default: > %(default)f]")
group.add_argument("--nbChr23", type=int, metavar="INT", default=50,
                   help=("The minimum number of markers on chromosome 23 "
                         "before computing Plink's sex check [default: "
                         "%(default)d]"))
group = parser.add_argument_group("Gender Plot")
group.add_argument("--gender-plot", action="store_true",
                   help=("Create the gender plot (summarized chr Y "
                         "intensities in function of summarized chr X "
                         "intensities) for problematic samples."))
group.add_argument("--sex-chr-intensities", type=str, metavar="FILE",
                   help=("A file containing alleles intensities for each of "
                         "the markers located on the X and Y chromosome for "
                         "the gender plot."))
group.add_argument("--gender-plot-format", type=str, metavar="FORMAT",
                   default="png", choices=["png", "ps", "pdf", "X11"],
                   help=("The output file format for the gender plot (png, "
                         "ps, pdf, or X11 formats are available). "
                         "[default: %(default)s]"))
group = parser.add_argument_group("LRR and BAF Plot")
group.add_argument("--lrr-baf", action="store_true",
                   help=("Create the LRR and BAF plot for problematic "
                         "samples"))
group.add_argument("--lrr-baf-raw-dir", type=str, metavar="DIR",
                   help=("Directory or list of directories containing "
                         "information about every samples (BAF and LRR)."))
group.add_argument("--lrr-baf-format", type=str, metavar="FORMAT",
                   default="png", choices=["png", "ps", "pdf", "X11"],
                   help=("The output file format for the LRR and BAF plot "
                         "(png, ps, pdf, or X11 formats are available). "
                         "[default: %(default)s]"))
group.add_argument("--lrr-baf-dpi", type=int, metavar="DPI", default=300,
                   help=("The pixel density of the figure(s) (DPI). "
                         "[default: %(default)d]"))
# The OUTPUT files
group = parser.add_argument_group("Output File")
group.add_argument("--out", type=str, metavar="FILE", default="sexcheck",
                   help=("The prefix of the output files (which will be a "
                         "Plink binary file. [default: %(default)s]"))


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
