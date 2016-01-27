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
import gzip
import logging
import argparse
import subprocess
from collections import defaultdict

import numpy as np

from .. import __version__
from ..PlinkUtils import createRowFromPlinkSpacedOutput
from ..RelatedSamples.merge_related_samples import merge_related_samples


logger = logging.getLogger("find_related_samples")


def main(argString=None):
    """The main function of this module.

    :param argString: the options.

    :type argString: list

    Here are the steps for this function:

    1. Prints the options.
    2. Uses Plink to extract markers according to LD
       (:py:func:`selectSNPsAccordingToLD`).
    3. Checks if there is enough markers after pruning
       (:py:func:`checkNumberOfSNP`). If not, then quits.
    4. Extract markers according to LD (:py:func:`extractSNPs`).
    5. Runs Plink with the ``genome`` option (:py:func:`runGenome`). Quits here
       if the user asker only for the ``genome`` file.
    6. Finds related individuals and gets values for plotting
       (:py:func:`extractRelatedIndividuals`).
    7. Plots ``Z1`` in function of ``IBS2 ratio`` for related individuals
       (:py:func:`plot_related_data`).
    8. Plots ``Z2`` in function of ``IBS2 ratio`` for related individuals
       (:py:func:`plot_related_data`).

    """
    # Getting and checking the options
    args = parseArgs(argString)
    checkArgs(args)

    logger.info("Options used:")
    for key, value in vars(args).iteritems():
        logger.info("  --{} {}".format(key.replace("_", "-"), value))

    # Run plink
    logger.info("Running Plink to extract SNPs according to LD")
    snpsToExtract = selectSNPsAccordingToLD(args)

    # Check there is enough SNP in the output file
    logger.info("Checking if there are enough extracted SNP")
    if not checkNumberOfSNP(snpsToExtract, args.min_nb_snp):
        # There are not enough markers
        logger.info("There are not enough SNPs: STOPPING NOW!")

    else:
        # Extract the SNPs
        logger.info("Extracting the SNPs using Plink")
        newBfile = extractSNPs(snpsToExtract, args)

        # Run the genome command from plink
        logger.info("Creating the genome file using Plink")
        genomeFileName = runGenome(newBfile, args)

        if args.genome_only:
            # We just want the genome file
            return newBfile

        # Extract related individuals
        logger.info("Finding related individuals from genome file")
        related_data = extractRelatedIndividuals(genomeFileName, args.out,
                                                 args.ibs2_ratio)
        # Are there related samples?
        if related_data is None:
            logger.info("There are no related samples in the dataset")

        else:
            # Plot the related data
            logger.info("Plotting related individuals")
            plot_related_data(related_data["IBS2_RATIO"], related_data["Z1"],
                              related_data["CODE"], r"$Z_1$",
                              args.out + ".related_individuals_z1.png", args)
            plot_related_data(related_data["IBS2_RATIO"], related_data["Z2"],
                              related_data["CODE"], r"$Z_2$",
                              args.out + ".related_individuals_z2.png", args)


def plot_related_data(x, y, code, ylabel, fileName, options):
    """Plot Z1 and Z2 in function of IBS2* ratio.

    :param x: the x axis of the plot (``IBS2 ratio``).
    :param y: the y axis of the plot (either ``z1`` or ``z2``).
    :param code: the code of the relatedness of each sample pair.
    :param ylabel: the label of the y axis (either ``z1`` or ``z2``).
    :param fileName: the name of the output file.
    :param options: the options.

    :type x: numpy.array of floats
    :type y: numpy.array of floats
    :type code: numpy.array
    :type ylabel: str
    :type fileName: str
    :type options: argparse.Namespace

    There are four different relation codes (represented by 4 different color
    in the plots:

    ==== =============================================== ===========
    Code                    Relation                        Color
    ==== =============================================== ===========
    1    Full-sbis                                       ``#CC0000``
    2    Half-sibs or Grand-parent-Child or Uncle-Nephew ``#0099CC``
    3    Parent-Child                                    ``#FF8800``
    4    Twins or Duplicated samples                     ``#9933CC``
    ==== =============================================== ===========

    Sample pairs with unknown relation are plotted using ``#669900`` as color.

    """
    import matplotlib as mpl
    if mpl.get_backend() != "agg":
        mpl.use("Agg")
    import matplotlib.pyplot as plt
    plt.ioff()

    fig = plt.figure()
    ax = fig.add_subplot(111)

    # Setting the title, the X and Y label
    ax.set_title((r"%d pairs with $IBS2^\ast_{ratio} >$ "
                  r"%f" % (len(code), options.ibs2_ratio)))
    ax.set_xlabel(r"$IBS2^\ast_{ratio}$")
    ax.set_ylabel(ylabel)

    # Plotting the data (there are 5 error codes)
    c5, = ax.plot(x[code == "5"], y[code == "5"], "o", ms=3, mec="#669900",
                  mfc="#669900")
    c1, = ax.plot(x[code == "1"], y[code == "1"], "o", ms=3, mec="#CC0000",
                  mfc="#CC0000")
    c2, = ax.plot(x[code == "2"], y[code == "2"], "o", ms=3, mec="#0099CC",
                  mfc="#0099CC")
    c3, = ax.plot(x[code == "3"], y[code == "3"], "o", ms=3, mec="#FF8800",
                  mfc="#FF8800")
    c4, = ax.plot(x[code == "4"], y[code == "4"], "o", ms=3, mec="#9933CC",
                  mfc="#9933CC")

    # The legend
    prop = mpl.font_manager.FontProperties(size=8)
    leg = ax.legend([c1, c2, c3, c4, c5],
                    ["Full sibs (n={})".format(np.sum(code == "1")),
                     ("Half sibs, grand-parent-child or uncle-nephew "
                      "(n={})".format(np.sum(code == "2"))),
                     "Parent-child (n={})".format(np.sum(code == "3")),
                     ("Twins or duplicated samples "
                      "(n={})".format(np.sum(code == "4"))),
                     "Unknown (n={})".format(np.sum(code == "5"))],
                    loc="best", numpoints=1, fancybox=True, prop=prop)
    leg.get_frame().set_alpha(0.5)

    # Setting the limits
    ax.set_xlim((options.ibs2_ratio - 0.01, 1.01))
    ax.set_ylim((-0.01, 1.01))

    # Modifying the spines
    ax.xaxis.set_ticks_position("bottom")
    ax.yaxis.set_ticks_position("left")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    # Saving the figure
    plt.savefig(fileName)


def extractRelatedIndividuals(fileName, outPrefix, ibs2_ratio_threshold):
    """Extract related individuals according IBS2* ratio.

    :param fileName: the name of the input file.
    :param outPrefix: the prefix of the output files.
    :param ibs2_ratio_threshold: the ibs2 ratio threshold (tells if sample pair
                                 is related or not).

    :type fileName: str
    :type outPrefix: str
    :type ibs2_ratio_threshold: float

    :returns: a :py:class:`numpy.recarray` data set containing (for each
              related sample pair) the ``ibs2 ratio``, ``Z1``, ``Z2`` and the
              type of relatedness.

    Reads a ``genome`` file (provided by :py:func:`runGenome`) and extract
    related sample pairs according to ``IBS2 ratio``.

    A ``genome`` file contains at least the following information for each
    sample pair:

    * **FID1:** the family ID of the first sample in the pair.
    * **IID1:** the individual ID of the first sample in the pair.
    * **FID2:** the family ID of the second sample in the pair.
    * **IID2:** the individual ID of the second sample in the pair.
    * **Z0:** the probability that :math:`IBD = 0`.
    * **Z1:** the probability that :math:`IBD = 1`.
    * **Z2:** the probability that :math:`IBD = 2`.
    * **HOMHOM:** the number of :math:`IBS = 0` SNP pairs used in ``PPC`` test.
    * **HETHET:** the number of :math:`IBS = 2` het/het SNP pairs in ``PPC``
      test.

    The ``IBS2 ratio`` is computed using the following formula:

    .. math::
        \\textrm{IBS2 ratio} = \\frac{\\textrm{HETHET}}
                                     {\\textrm{HOMHOM} + \\textrm{HETHET}}

    If the ``IBS2 ratio`` is higher than the threshold, the samples in the pair
    are related. The following values help in finding the relatedness of the
    sample pair.

    +---------------------------------------+-----------------------+------+
    |           Values                      |        Relation       | Code |
    +=======================================+=======================+======+
    | :math:`0.17 \\leq z_0 \\leq 0.33` and   | Full-sibs             | 1    |
    | :math:`0.40 \\leq z_1 \\leq 0.60`       |                       |      |
    +---------------------------------------+-----------------------+------+
    | :math:`0.40 \\leq z_0 \\leq 0.60` and   | Half-sibs or          | 2    |
    | :math:`0.40 \\leq z_1 \\leq 0.60`       | Grand-parent-Child or |      |
    |                                       | Uncle-Nephew          |      |
    +---------------------------------------+-----------------------+------+
    | :math:`z_0 \\leq 0.05` and             | Parent-Child          | 3    |
    | :math:`z_1 \\geq 0.95` and             |                       |      |
    | :math:`z_2 \\leq 0.05`                 |                       |      |
    +---------------------------------------+-----------------------+------+
    | :math:`z_0 \\leq 0.05` and             | Twins or Duplicated   | 4    |
    | :math:`z_1 \\leq 0.05` and             | samples               |      |
    | :math:`z_2 \\geq 0.95`                 |                       |      |
    +---------------------------------------+-----------------------+------+

    """
    # The output file
    outputFile = None
    try:
        outputFile = open(outPrefix + ".related_individuals", "w")
    except IOError:
        msg = "%(outPrefix)s.related_individuals: can't write file" % locals()
        raise ProgramError(msg)

    # The input file
    inputFile = None
    try:
        if fileName.endswith(".gz"):
            inputFile = gzip.open(fileName, 'rb')
        else:
            inputFile = open(fileName, 'r')
    except IOError:
        msg = "{}: no such file".format(fileName)
        raise ProgramError(msg)

    headerIndex = None
    data = []
    for i, line in enumerate(inputFile):
        row = createRowFromPlinkSpacedOutput(line)

        if i == 0:
            # This is the header
            headerIndex = dict([(colName, j) for j, colName in enumerate(row)])

            # Checking columns
            for columnName in ["FID1", "IID1", "FID2", "IID2", "Z0", "Z1",
                               "Z2", "HOMHOM", "HETHET"]:
                if columnName not in headerIndex:
                    msg = "{}: no culumn named {}".format(fileName, columnName)
                    raise ProgramError(msg)

            # Writing header
            print >>outputFile, "\t".join(row + ["IBS2_ratio", "status",
                                                 "code"])

        else:
            # This is data
            homhom = row[headerIndex["HOMHOM"]]
            hethet = row[headerIndex["HETHET"]]
            try:
                homhom = float(homhom)
                hethet = float(hethet)
            except ValueError:
                msg = "{}: invalid HOMHOM or HETHET".format(fileName)
                raise ProgramError(float)

            # Computing IBS2* ratio
            ibs2_ratio = hethet / (homhom + hethet)

            if ibs2_ratio > ibs2_ratio_threshold:
                # Those pairs might be related
                # Finding the status
                status = "unknown"
                code = "5"
                z0 = row[headerIndex["Z0"]]
                z1 = row[headerIndex["Z1"]]
                z2 = row[headerIndex["Z2"]]
                try:
                    z0 = float(z0)
                    z1 = float(z1)
                    z2 = float(z2)
                except ValueError:
                    msg = "{}: invalid value for Z0, Z1 or Z2".format(fileName)
                    raise ProgramError(msg)

                if (z0 >= 0.17 and z0 <= 0.33) and (z1 >= 0.40 and z1 <= 0.60):
                    # Full sibs
                    status = "full-sibs"
                    code = "1"

                elif (z0 >= 0.4 and z0 <= 0.6) and (z1 >= 0.4 and z1 <= 0.6):
                    # half sibs, grand-parent child, uncle nephew
                    status = ";".join(["half-sibs", "grand-parent-child",
                                       "uncle-nephew"])
                    code = "2"

                elif (z0 <= 0.05) and (z1 >= 0.95) and (z2 <= 0.05):
                    # parent child
                    status = "parent-child"
                    code = "3"

                elif (z0 <= 0.05) and (z1 <= 0.05) and (z2 >= 0.95):
                    # twin
                    status = "twins"
                    code = "4"

                # Printing to file
                print >>outputFile, "\t".join(row + [str(ibs2_ratio), status,
                                                     code])
                data.append((ibs2_ratio, z1, z2, code))

    # Closing the output and input files
    inputFile.close()
    outputFile.close()

    # Merging the related individuals
    merge_related_samples(outPrefix + ".related_individuals", outPrefix, False)

    # If there are no related samples, we return nothing
    if len(data) == 0:
        return None

    # Creating the numpy array if there are related samples
    data = np.array(data, dtype=[
        ("IBS2_RATIO", float),
        ("Z1", float),
        ("Z2", float),
        ("CODE", "S{}".format(max([len(i[3]) for i in data]))),
    ])

    return data


def checkNumberOfSNP(fileName, minimumNumber):
    """Check there is enough SNPs in the file (with minimum).

    :param fileName: the name of the file.
    :param minimumNumber: the minimum number of markers that needs to be in the
                          file.

    :type fileName: str
    :type minimumNumber: int

    :returns: ``True`` if there is enough markers in the file, ``False``
              otherwise.

    Reads the number of markers (number of lines) in a file.

    """
    nbSNP = 0
    try:
        with open(fileName, 'r') as inputFile:
            for line in inputFile:
                nbSNP += 1
    except IOError:
        msg = "{}: no such file".format(fileName)
        raise ProgramError(msg)

    if nbSNP < minimumNumber:
        return False
    return True


def splitFile(inputFileName, linePerFile, outPrefix):
    """Split a file.

    :param inputFileName: the name of the input file.
    :param linePerFile: the number of line per file (after splitting).
    :param outPrefix: the prefix of the output files.

    :type inputFileName: str
    :type linePerFile: int
    :type outPrefix: str

    :returns: the number of created temporary files.

    Splits a file (``inputFileName`` into multiple files containing at most
    ``linePerFile`` lines.

    """
    nbTmpFile = 1
    nbLine = 0
    tmpFile = None
    try:
        with open(inputFileName, "r") as inputFile:
            for line in inputFile:
                row = line.rstrip("\r\n").split(" ")
                nbLine += 1

                if tmpFile is None:
                    try:
                        tmpFile = open(
                            outPrefix + "_tmp.list%d" % nbTmpFile,
                            "w",
                        )
                    except IOError:
                        msg = "tmp.list%d: can't write file" % nbTmpFile
                        raise ProgramError(msg)

                print >>tmpFile, " ".join(row[:2])

                if nbLine == linePerFile:
                    nbLine = 0
                    nbTmpFile += 1
                    tmpFile.close()
                    try:
                        tmpFile = open(
                            outPrefix + "_tmp.list%d" % nbTmpFile,
                            "w",
                        )
                    except IOError:
                        msg = "tmp.list%d: can't write file" % nbTmpFile
                        raise ProgramError(msg)
        tmpFile.close()

        # Check if the number of line is zero (hence the last file is empty)
        if nbLine == 0:
            # We delete the last file
            file_name = outPrefix + "_tmp.list{}".format(nbTmpFile)
            if os.path.isfile(file_name):
                os.remove(file_name)
            nbTmpFile -= 1

    except IOError:
        msg = "%s: no such file" % inputFileName
        raise ProgramError(msg)

    return nbTmpFile


def runGenome(bfile, options):
    """Runs the genome command from plink.

    :param bfile: the input file prefix.
    :param options: the options.

    :type bfile: str
    :type options: argparse.Namespace

    :returns: the name of the ``genome`` file.

    Runs Plink with the ``genome`` option. If the user asks for SGE
    (``options.sge`` is True), a frequency file is first created by plink.
    Then, the input files are split in ``options.line_per_file_for_sge`` and
    Plink is called (using the ``genome`` option) on the cluster using SGE
    (:py:func:`runGenomeSGE`). After the analysis, Plink's output files and
    logs are merged using :py:func:`mergeGenomeLogFiles`.

    """
    outPrefix = options.out + ".genome"
    if options.sge:
        # We run genome using SGE
        # We need to create a frequency file using plink
        plinkCommand = ["plink", "--noweb", "--bfile", bfile, "--freq",
                        "--out", options.out + ".frequency"]
        runCommand(plinkCommand)

        # We need to split the .fam file
        nbJob = splitFile(bfile + ".fam", options.line_per_file_for_sge,
                          outPrefix)

        runGenomeSGE(bfile, options.out + ".frequency.frq", nbJob,
                     outPrefix, options)

        # Merging genome files
        mergeGenomeLogFiles(outPrefix, nbJob)

    else:
        plinkCommand = ["plink", "--noweb", "--bfile", bfile, "--genome",
                        "--genome-full", "--out", outPrefix]
        runCommand(plinkCommand)

    return outPrefix + ".genome"


def mergeGenomeLogFiles(outPrefix, nbSet):
    """Merge genome and log files together.

    :param outPrefix: the prefix of the output files.
    :param nbSet: The number of set of files to merge together.

    :type outPrefix: str
    :type nbSet: int

    :returns: the name of the output file (the ``genome`` file).

    After merging, the files are deleted to save space.

    """
    outputFile = None
    try:
        outputFile = open(outPrefix + ".genome", "w")
        outputLog = open(outPrefix + ".log", "w")
    except IOError:
        msg = "%s or %s: can't write file" % (outPrefix + ".genome",
                                              outPrefix + ".log")
        raise ProgramError(msg)

    for i in xrange(1, nbSet + 1):
        for j in xrange(i, nbSet + 1):
            fileName = outPrefix + "_output.sub.%(i)d.%(j)d.genome" % locals()

            printHeader = False
            if (i == 1) and (j == 1):
                # This is the first file we open
                printHeader = True

            # Read file here
            try:
                with open(fileName, 'r') as inputFile:
                    for nbLine, line in enumerate(inputFile):
                        if nbLine == 0:
                            if printHeader:
                                outputFile.write(line)
                        else:
                            outputFile.write(line)
            except IOError:
                msg = "%(fileName)s: no such file" % locals()
                raise ProgramError(msg)

            # Deleting the file
            try:
                os.remove(fileName)
            except IOError:
                msg = "%(fileName)s: can't delete the file" % locals()
                raise ProgramError(msg)

            # Read file here
            fileName = outPrefix + "_output.sub.%(i)d.%(j)d.log" % locals()
            try:
                with open(fileName, 'r') as inputFile:
                    for line in inputFile:
                        outputLog.write(line)
            except IOError:
                msg = "%(fileName)s: no such file" % locals()
                raise ProgramError(msg)

            # Deleting the file
            try:
                os.remove(fileName)
            except IOError:
                msg = "%(fileName)s: can't delete the file" % locals()
                raise ProgramError(msg)

    # Removing the tmp.list* files
    try:
        for fileName in glob.glob(outPrefix + "_tmp.list*"):
            os.remove(fileName)
    except IOError:
        msg = "can't delete the tmp.list* files"
        raise ProgramError(msg)

    # Removing the output.sub.*
    try:
        for fileName in glob.glob(outPrefix + "_output.sub.*"):
            os.remove(fileName)
    except IOError:
        msg = "can't delete the output.sub.* files"
        raise ProgramError(msg)

    # Closing the output files
    outputFile.close()
    outputLog.close()

    return outPrefix + ".genome"


def runGenomeSGE(bfile, freqFile, nbJob, outPrefix, options):
    """Runs the genome command from plink, on SGE.

    :param bfile: the prefix of the input file.
    :param freqFile: the name of the frequency file (from Plink).
    :param nbJob: the number of jobs to launch.
    :param outPrefix: the prefix of all the output files.
    :param options: the options.

    :type bfile: str
    :type freqFile: str
    :type nbJob: int
    :type outPrefix: str
    :type options: argparse.Namespace

    Runs Plink with the ``genome`` options on the cluster (using SGE).

    """
    # Add the environment variable for DRMAA package
    if "DRMAA_LIBRARY_PATH" not in os.environ:
        msg = "could not load drmaa: set DRMAA_LIBRARY_PATH"
        raise ProgramError(msg)

    # Import the python drmaa library
    try:
        import drmaa
    except ImportError:
        raise ProgramError("drmaa is not install, install drmaa")

    # Initializing a session
    s = drmaa.Session()
    s.initialize()

    # Run for each sub task...
    jobIDs = []
    jobTemplates = []
    for i in xrange(1, nbJob + 1):
        for j in xrange(i, nbJob + 1):
            # The command to run
            plinkCommand = ["plink", "--noweb", "--bfile", bfile,
                            "--read-freq", freqFile, "--genome",
                            "--genome-full", "--genome-lists",
                            "{}_tmp.list{}".format(outPrefix, i),
                            "{}_tmp.list{}".format(outPrefix, j), "--out",
                            "{}_output.sub.{}.{}".format(outPrefix, i, j)]

            # Creating the job template
            jt = s.createJobTemplate()
            jt.remoteCommand = plinkCommand[0]
            jt.workingDirectory = os.getcwd()
            jt.jobEnvironment = os.environ
            jt.args = plinkCommand[1:]
            jt.jobName = "_plink_genome_{}_{}".format(i, j)

            # Cluster specifics
            if options.sge_walltime is not None:
                jt.hardWallclockTimeLimit = options.sge_walltime
            if options.sge_nodes is not None:
                native_spec = "-l nodes={}:ppn={}".format(options.sge_nodes[0],
                                                          options.sge_nodes[1])
                jt.nativeSpecification = native_spec

            jobIDs.append(s.runJob(jt))
            jobTemplates.append(jt)

    # Waiting for the jobs to finish
    hadProblems = []
    for jobID in jobIDs:
        retVal = s.wait(jobID, drmaa.Session.TIMEOUT_WAIT_FOREVER)
        hadProblems.append(retVal.exitStatus == 0)

    # Deleting the jobs
    for jt in jobTemplates:
        s.deleteJobTemplate(jt)

    # Closing the session
    s.exit()

    # Checking for problems
    for hadProblem in hadProblems:
        if not hadProblem:
            msg = "Some SGE jobs had errors..."
            raise ProgramError(msg)


def extractSNPs(snpsToExtract, options):
    """Extract markers using Plink.

    :param snpsToExtract: the name of the file containing markers to extract.
    :param options: the options

    :type snpsToExtract: str
    :type options: argparse.Namespace

    :returns: the prefix of the output files.

    """
    outPrefix = options.out + ".pruned_data"
    plinkCommand = ["plink", "--noweb", "--bfile", options.bfile, "--extract",
                    snpsToExtract, "--make-bed", "--out", outPrefix]
    runCommand(plinkCommand)
    return outPrefix


def selectSNPsAccordingToLD(options):
    """Compute LD using Plink.

    :param options: the options.

    :type options: argparse.Namespace

    :returns: the name of the output file (from Plink).

    """
    # The plink command
    outPrefix = options.out + ".pruning_" + options.indep_pairwise[2]
    plinkCommand = [
        "plink",
        "--noweb",
        "--bfile", options.bfile,
        "--maf", options.maf,
        "--indep-pairwise",
    ] + options.indep_pairwise + ["--out", outPrefix]

    runCommand(plinkCommand)

    # Finding the autosomal markers
    autosomes = {str(i) for i in range(1, 23)}
    autosomal_snps = set()
    with open(options.bfile + ".bim", "r") as i_file:
        for line in i_file:
            chrom, snp = line.rstrip("\r\n").split("\t")[:2]
            if chrom in autosomes:
                autosomal_snps.add(snp)

    # Reading the pruned markers
    pruned_snps = None
    with open(outPrefix + ".prune.in", "r") as i_file:
        pruned_snps = set(i_file.read().splitlines())

    # Writing the pruned markers located on an autosome
    with open(outPrefix + ".prune.in.autosomal", "w") as o_file:
        for snp in autosomal_snps & pruned_snps:
            print >>o_file, snp

    return outPrefix + ".prune.in.autosomal"


def runCommand(command):
    """Run a command.

    :param command: the command to run.

    :type command: list

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

    # Check the indep-pairwise option
    # The two first must be int, the last one float
    try:
        for i in xrange(2):
            tmp = int(args.indep_pairwise[i])
        tmp = float(args.indep_pairwise[2])
    except ValueError:
        msg = "indep-pairwise: need INT INT FLOAT"
        raise ProgramError(msg)

    # Check the maf value
    tmpMAF = None
    try:
        tmpMAF = float(args.maf)
    except ValueError:
        msg = "maf: must be a float, not %s" % args.maf
        raise ProgramError(msg)
    if (tmpMAF > 0.5) or (tmpMAF < 0.0):
        msg = "maf: must be between 0.0 and 0.5, not %s" % args.maf
        raise ProgramError(msg)

    # Check the number of line per file
    if args.line_per_file_for_sge < 1:
        msg = "line-per-file-for-sge: must be above 0, not " \
              "%d" % args.line_per_file_for_sge
        raise ProgramError(msg)

    # Check the minimum number of SNPs
    if args.min_nb_snp < 1:
        msg = "min-nb-snp: must be above 1"
        raise ProgramError(msg)

    return True


def parseArgs(argString=None):  # pragma: no cover
    """Parses the command line options and arguments.

    :param argString: the options.

    :type argString: list

    :returns: A :py:class:`argparse.Namespace` object created by
              the :py:mod:`argparse` module. It contains the values of the
              different options.

    =========================== ====== ========================================
        Options       Type                     Description
    =========================== ====== ========================================
    ``--bfile``                 string The input file prefix (Plink binary
                                       file).
    ``--genome-only``           bool   Only create the genome file.
    ``--min-nb-snp``            int    The minimum number of markers needed to
                                       compute IBS values.
    ``--indep-pairwise``        string Three numbers: window size, window shift
                                       and the r2 threshold.
    ``--maf``                   string Restrict to SNPs with MAF >= threshold.
    ``--ibs2-ratio``            float  The initial IBS2* ratio (the minimum
                                       value to show in the plot.
    ``--sge``                   bool   Use SGE for parallelization.
    ``--sge-walltime``          int    The time limit (for clusters).
    ``--sge-nodes``             int    Two INTs (number of nodes and number of
                                int    processor per nodes).
    ``--line-per-file-for-sge`` int    The number of line per file for SGE task
                                       array.
    ``--out``                   string The prefix of the output files.
    =========================== ====== ========================================

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


# Default values
_ibs2_ratio_default = 0.8
_indep_pairwise_r2_default = "0.1"

# The parser object
pretty_name = "Related samples"
desc = "Finds related samples according to IBS values."
long_desc = ("The script conducts close familial relationship checks with "
             "pairwise IBD. It flags and removes all but one pair-member of "
             r"samples duplicates ($IBS2^\ast_{ratio}>{ratio_value}$) based "
             "on a selection of uncorrelated markers ($r^2 < {r_squared}$).")
parser = argparse.ArgumentParser(description=desc)
parser.add_argument("-v", "--version", action="version",
                    version="pyGenClean version {}".format(__version__))

# The INPUT files
group = parser.add_argument_group("Input File")
group.add_argument("--bfile", type=str, metavar="FILE", required=True,
                   help=("The input file prefix (will find the plink binary "
                         "files by appending the prefix to the .bim, .bed and "
                         ".fam files, respectively.)"))
# The options
group = parser.add_argument_group("Options")
group.add_argument("--genome-only", action="store_true",
                   help="Only create the genome file")
group.add_argument("--min-nb-snp", type=int, metavar="INT", default=10000,
                   help=("The minimum number of markers needed to compute IBS "
                         "values. [Default: %(default)d]"))
group.add_argument("--indep-pairwise", type=str, metavar="STR", nargs=3,
                   default=["50", "5", _indep_pairwise_r2_default],
                   help=("Three numbers: window size, window shift and the r2 "
                         "threshold. [default: %(default)s]"))
group.add_argument("--maf", type=str, metavar="FLOAT", default="0.05",
                   help=("Restrict to SNPs with MAF >= threshold. [default: "
                         "%(default)s]"))
group.add_argument("--ibs2-ratio", type=float, metavar="FLOAT",
                   default=_ibs2_ratio_default,
                   help=("The initial IBS2* ratio (the minimum value to show "
                         "in the plot. [default: %(default).1f]"))
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
                         "the number of processor to use. Do not use if you "
                         "are not required to specify the number of nodes for "
                         "your jobs on the cluster."))
group.add_argument("--line-per-file-for-sge", type=int, metavar="INT",
                   default=100, help=("The number of line per file for SGE "
                                      "task array. [default: " "%(default)d]"))
# The OUTPUT files
group = parser.add_argument_group("Output File")
group.add_argument("--out", type=str, metavar="FILE", default="ibs",
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
