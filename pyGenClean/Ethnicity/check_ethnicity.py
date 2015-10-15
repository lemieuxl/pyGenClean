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
from itertools import izip
from collections import defaultdict

from ..PlinkUtils import plot_MDS as PlotMDS
from ..PlinkUtils import createRowFromPlinkSpacedOutput

from .. import __version__
from . import find_outliers
from . import plot_eigenvalues as PlotEigenvalues
from ..DupSNPs.duplicated_snps import flipGenotype
from ..RelatedSamples import find_related_samples as Relatedness


logger = logging.getLogger("check_ethnicity")


class Dummy(object):
    pass


def main(argString=None):
    """The main function.

    :param argString: the options.

    :type argString: list

    These are the steps of this module:

    1.  Prints the options.
    2.  Finds the overlapping markers between the three reference panels and
        the source panel (:py:func:`findOverlappingSNPsWithReference`).
    3.  Extract the required markers from all the data sets
        (:py:func:`extractSNPs`).
    4.  Renames the reference panel's marker names to that they are the same as
        the source panel (for all populations) (:py:func:`renameSNPs`).
    5.  Combines the three reference panels together
        (:py:func:`combinePlinkBinaryFiles`).
    6.  Compute the frequency of all the markers from the reference and the
        source panels (:py:func:`computeFrequency`).
    7.  Finds the markers to flip in the reference panel (when compared to the
        source panel) (:py:func:`findFlippedSNPs`).
    8.  Excludes the unflippable markers from the reference and the source
        panels (:py:func:`excludeSNPs`).
    9.  Flips the markers that need flipping in their reference panel
        (:py:func:`flipSNPs`).
    10. Combines the reference and the source panels
        (:py:func:`combinePlinkBinaryFiles`).
    11. Runs part of :py:mod:`pyGenClean.RelatedSamples.find_related_samples`
        on the combined data set (:py:func:`runRelatedness`).
    12. Creates the ``mds`` file from the combined data set and the result of
        previous step (:py:func:`createMDSFile`).
    13. Creates the population file (:py:func:`createPopulationFile`).
    14. Plots the ``mds`` values (:py:func:`plotMDS`).
    15. Finds the outliers of a given reference population
        (:py:func:`find_the_outliers`).
    16. If required, computes the Eigenvalues using smartpca
        (:py:func:`compute_eigenvalues`).
    17. If required, creates a scree plot from smartpca resutls
        (:py:func:`create_scree_plot`).

    """
    # Getting and checking the options
    args = parseArgs(argString)
    checkArgs(args)

    logger.info("Options used:")
    for key, value in vars(args).iteritems():
        logger.info("  --{} {}".format(key.replace("_", "-"), value))

    newBfile = None
    popNames = ["CEU", "YRI", "JPT-CHB"]
    referencePrefixes = [args.ceu_bfile, args.yri_bfile, args.jpt_chb_bfile]
    if not args.skip_ref_pops:
        # Find overlap with the reference file
        logger.info("Finding overlapping SNPs between reference and "
                    "source panels")
        findOverlappingSNPsWithReference(
            prefix=args.bfile,
            referencePrefixes=referencePrefixes,
            referencePopulations=popNames,
            outPrefix=args.out,
        )

        # Extract the required SNPs using Plink (reference panels)
        logger.info("Extracting overlapping SNPs from the reference panels")
        extractSNPs(
            snpToExtractFileNames=[
                args.out + ".{}_snp_to_extract".format(popName)
                for popName in popNames
            ],
            referencePrefixes=referencePrefixes,
            popNames=popNames,
            outPrefix=args.out + ".reference_panel",
            runSGE=args.sge,
            options=args,
        )

        # Extract the required SNPs using Plink (source panel)
        logger.info("Extracting overlapping SNPs from the source panel")
        extractSNPs(
            snpToExtractFileNames=[args.out + ".source_snp_to_extract"],
            referencePrefixes=[args.bfile],
            popNames=["ALL"],
            outPrefix=args.out + ".source_panel",
            runSGE=False,
            options=args,
        )

        # Renaming the reference file, so that the SNP names are the same
        for pop in popNames:
            logger.info("Renaming reference panel's SNPs {} to match source "
                        "panel".format(pop))
            renameSNPs(
                inPrefix=args.out + ".reference_panel.{}".format(pop),
                updateFileName=args.out + ".{}_update_names".format(pop),
                outPrefix=args.out + ".reference_panel.{}.rename".format(pop),
            )

        # Combining the reference panel
        logger.info("Combining the reference panels")
        combinePlinkBinaryFiles(
            prefixes=[
                args.out + ".reference_panel.{}.rename".format(pop)
                for pop in popNames
            ],
            outPrefix=args.out + ".reference_panel.ALL.rename",
        )

        # Computing the frequency (reference panel)
        logger.info("Computing reference panel frequencies")
        computeFrequency(
            prefix=args.out + ".reference_panel.ALL.rename",
            outPrefix=args.out + ".reference_panel.ALL.rename.frequency",
        )

        # Computing the frequency (source panel)
        logger.info("Computing source panel frequencies")
        computeFrequency(
            prefix=args.out + ".source_panel.ALL",
            outPrefix=args.out + ".source_panel.ALL.frequency",
        )

        # Finding the SNPs to flip and flip them in the reference panel
        logger.info("Finding SNPs to flip or to exclude from reference panel")
        findFlippedSNPs(
            frqFile1=args.out + ".reference_panel.ALL.rename.frequency.frq",
            frqFile2=args.out + ".source_panel.ALL.frequency.frq",
            outPrefix=args.out,
        )

        # Excluding SNPs (reference panel)
        logger.info("Excluding SNPs from reference panel")
        excludeSNPs(
            inPrefix=args.out + ".reference_panel.ALL.rename",
            outPrefix=args.out + ".reference_panel.ALL.rename.cleaned",
            exclusionFileName=args.out + ".snp_to_remove",
        )

        # Excluding SNPs (source panel)
        logger.info("Excluding SNPs from source panel")
        excludeSNPs(args.out + ".source_panel.ALL",
                    args.out + ".source_panel.ALL.cleaned",
                    args.out + ".snp_to_remove")

        # Flipping the SNP that need to be flip in the reference
        logger.info("Flipping SNPs in reference panel")
        flipSNPs(
            inPrefix=args.out + ".reference_panel.ALL.rename.cleaned",
            outPrefix=args.out + ".reference_panel.ALL.rename.cleaned.flipped",
            flipFileName=args.out + ".snp_to_flip_in_reference",
        )

        # Combining the reference panel
        logger.info("Combining reference and source panels")
        combinePlinkBinaryFiles(
            prefixes=[args.out + ".reference_panel.ALL.rename.cleaned.flipped",
                      args.out + ".source_panel.ALL.cleaned"],
            outPrefix=args.out + ".final_dataset_for_genome",
        )

        # Runing the relatedness step
        logger.info("Creating the genome file using Plink")
        newBfile = runRelatedness(
            inputPrefix=args.out + ".final_dataset_for_genome",
            outPrefix=args.out,
            options=args,
        )

    else:
        # Just run relatedness on the dataset
        newBfile = runRelatedness(
            inputPrefix=args.bfile,
            outPrefix=args.out,
            options=args,
        )

    # Creating the MDS file
    logger.info("Creating the MDS file using Plink")
    createMDSFile(
        nb_components=args.nb_components,
        inPrefix=newBfile,
        outPrefix=args.out + ".mds",
        genomeFileName=args.out + ".ibs.genome.genome",
    )

    if not args.skip_ref_pops:
        # Creating the population files
        logger.info("Creating a population file")
        famFiles = [
            args.out + ".reference_panel." + i + ".fam" for i in popNames
        ]
        famFiles.append(args.out + ".source_panel.ALL.fam")
        labels = popNames + ["SOURCE"]
        createPopulationFile(
            inputFiles=famFiles,
            labels=labels,
            outputFileName=args.out + ".population_file",
        )

        # Plot the MDS value
        logger.info("Creating the MDS plot")
        plotMDS(
            inputFileName=args.out + ".mds.mds",
            outPrefix=args.out + ".mds",
            populationFileName=args.out + ".population_file",
            options=args,
        )

        # Finding the outliers
        logger.info("Finding the outliers")
        find_the_outliers(
            mds_file_name=args.out + ".mds.mds",
            population_file_name=args.out + ".population_file",
            ref_pop_name=args.outliers_of,
            multiplier=args.multiplier,
            out_prefix=args.out,
        )

    # De we need to create a scree plot?
    if args.create_scree_plot:
        # Computing the eigenvalues using smartpca
        logger.info("Computing eigenvalues")
        compute_eigenvalues(
            in_prefix=args.out + ".ibs.pruned_data",
            out_prefix=args.out + ".smartpca",
        )

        logger.info("Creating scree plot")
        create_scree_plot(
            in_filename=args.out + ".smartpca.evec.txt",
            out_filename=args.out + ".smartpca.scree_plot.png",
            plot_title=args.scree_plot_title,
        )


def create_scree_plot(in_filename, out_filename, plot_title):
    """Creates a scree plot using smartpca results.

    :param in_filename: the name of the input file.
    :param out_filename: the name of the output file.
    :param plot_title: the title of the scree plot.

    :type in_filename: str
    :type out_filename: str
    :type plot_title: str

    """
    # Constructing the arguments
    scree_plot_args = ("--evec", in_filename, "--out", out_filename,
                       "--scree-plot-title", plot_title)
    try:
        # Executing the main method
        PlotEigenvalues.main(argString=scree_plot_args)

    except PlotEigenvalues.ProgramError as e:
        msg = "PlotEigenvalues: {}".format(e)
        raise ProgramError(msg)


def compute_eigenvalues(in_prefix, out_prefix):
    """Computes the Eigenvalues using smartpca from Eigensoft.

    :param in_prefix: the prefix of the input files.
    :param out_prefix: the prefix of the output files.

    :type in_prefix: str
    :type out_prefix: str

    Creates a "parameter file" used by smartpca and runs it.

    """
    # First, we create the parameter file
    with open(out_prefix + ".parameters", "w") as o_file:
        print >>o_file, "genotypename:    " + in_prefix + ".bed"
        print >>o_file, "snpname:         " + in_prefix + ".bim"
        print >>o_file, "indivname:       " + in_prefix + ".fam"
        print >>o_file, "evecoutname:     " + out_prefix + ".evec.txt"
        print >>o_file, "evaloutname:     " + out_prefix + ".eval.txt"
        print >>o_file, "numoutlieriter:  0"
        print >>o_file, "altnormstyle:    NO"

    # Executing smartpca
    command = ["smartpca", "-p", out_prefix + ".parameters"]
    runCommand(command)


def find_the_outliers(mds_file_name, population_file_name, ref_pop_name,
                      multiplier, out_prefix):
    """Finds the outliers of a given population.

    :param mds_file_name: the name of the ``mds`` file.
    :param population_file_name: the name of the population file.
    :param ref_pop_name: the name of the reference population for which to find
                         outliers from.
    :param multiplier: the multiplier of the cluster standard deviation to
                       modify the strictness of the outlier removal procedure.
    :param out_prefix: the prefix of the output file.

    :type mds_file_name: str
    :type population_file_name: str
    :type ref_pop_name: str
    :type multiplier: float
    :type out_prefix: str

    Uses the :py:mod:`pyGenClean.Ethnicity.find_outliers` modules to find
    outliers. It requires the ``mds`` file created by :py:func:`createMDSFile`
    and the population file created by :py:func:`createPopulationFile`.

    """
    options = ["--mds", mds_file_name, "--population-file",
               population_file_name, "--outliers-of", ref_pop_name,
               "--multiplier", str(multiplier), "--out", out_prefix]

    try:
        find_outliers.main(options)
    except find_outliers.ProgramError as e:
        msg = "find_outliers: {}".format(e)
        raise ProgramError(msg)


def createPopulationFile(inputFiles, labels, outputFileName):
    """Creates a population file.

    :param inputFiles: the list of input files.
    :param labels: the list of labels (corresponding to the input files).
    :param outputFileName: the name of the output file.

    :type inputFiles: list
    :type labels: list
    :type outputFileName: str

    The ``inputFiles`` is in reality a list of ``tfam`` files composed of
    samples. For each of those ``tfam`` files, there is a label associated with
    it (representing the name of the population).

    The output file consists of one row per sample, with the following three
    columns: the family ID, the individual ID and the population of each
    sample.

    """
    outputFile = None
    try:
        outputFile = open(outputFileName, 'w')
    except IOError:
        msg = "%(outputFileName)s: can't write file"
        raise ProgramError(msg)

    for i in xrange(len(inputFiles)):
        # For each file
        fileName = inputFiles[i]
        label = labels[i]

        try:
            with open(fileName, 'r') as inputFile:
                for line in inputFile:
                    row = line.rstrip("\r\n").split(" ")

                    # Getting the informations
                    famID = row[0]
                    indID = row[1]

                    # Printing to file
                    print >>outputFile, "\t".join([famID, indID, label])
        except IOError:
            msg = "%(fileName)s: no such file" % locals()
            raise ProgramError(msg)

    # Closing the output file
    outputFile.close()


def plotMDS(inputFileName, outPrefix, populationFileName, options):
    """Plots the MDS value.

    :param inputFileName: the name of the ``mds`` file.
    :param outPrefix: the prefix of the output files.
    :param populationFileName: the name of the population file.
    :param options: the options

    :type inputFileName: str
    :type outPrefix: str
    :type populationFileName: str
    :type options: argparse.Namespace

    Plots the ``mds`` value according to the ``inputFileName`` file (``mds``)
    and the ``populationFileName`` (the population file).

    """
    # The options
    plotMDS_options = Dummy()
    plotMDS_options.file = inputFileName
    plotMDS_options.out = outPrefix
    plotMDS_options.format = options.format
    plotMDS_options.title = options.title
    plotMDS_options.xlabel = options.xlabel
    plotMDS_options.ylabel = options.ylabel
    plotMDS_options.population_file = populationFileName

    try:
        # Checking the options
        PlotMDS.checkArgs(plotMDS_options)

        # Reading the population file
        populations = PlotMDS.readPopulations(plotMDS_options.population_file)

        # Getting the data
        data, labels = PlotMDS.extractData(plotMDS_options.file, populations)
        order = [labels.index("CEU"), labels.index("YRI"),
                 labels.index("JPT-CHB"), labels.index("SOURCE")]
        if "HIGHLIGHT" in labels:
            order.append(labels.index("HIGHLIGHT"))
        color = [(0.215686275, 0.000494118, 0.721568627),
                 (0.301960784, 0.68627451, 0.290196078),
                 (0.596078431, 0.305882353, 0.639215686),
                 (0.894117647, 0.101960784, 0.109803922),
                 (0.596078431, 0.305882353, 0.639215686)]
        sizes = [12, 12, 12, 8, 12]
        markers = [".", ".", ".", "+", "D"]

        # Plotting the data
        PlotMDS.plotMDS(data, order, labels, color, sizes, markers,
                        plotMDS_options)

    except PlotMDS.ProgramError as e:
        msg = "PlotMDS: %s" % e
        raise ProgramError(msg)


def createMDSFile(nb_components, inPrefix, outPrefix, genomeFileName):
    """Creates a MDS file using Plink.

    :param nb_components: the number of component.
    :param inPrefix: the prefix of the input file.
    :param outPrefix: the prefix of the output file.
    :param genomeFileName: the name of the ``genome`` file.

    :type nb_components: int
    :type inPrefix: str
    :type outPrefix: str
    :type genomeFileName: str

    Using Plink, computes the MDS values for each individual using the
    ``inPrefix``, ``genomeFileName`` and the number of components. The results
    are save using the ``outPrefix`` prefix.

    """
    plinkCommand = ["plink", "--noweb", "--bfile", inPrefix, "--read-genome",
                    genomeFileName, "--cluster", "--mds-plot",
                    str(nb_components), "--out", outPrefix]
    runCommand(plinkCommand)


def runRelatedness(inputPrefix, outPrefix, options):
    """Run the relatedness step of the data clean up.

    :param inputPrefix: the prefix of the input file.
    :param outPrefix: the prefix of the output file.
    :param options: the options

    :type inputPrefix: str
    :type outPrefix: str
    :type options: argparse.Namespace

    :returns: the prefix of the new bfile.

    Runs :py:mod:`pyGenClean.RelatedSamples.find_related_samples` using the
    ``inputPrefix`` files and ``options`` options, and saves the results using
    the ``outPrefix`` prefix.

    """
    # The options
    new_options = ["--bfile", inputPrefix, "--genome-only",
                   "--min-nb-snp", str(options.min_nb_snp),
                   "--maf", options.maf,
                   "--out", "{}.ibs".format(outPrefix)]
    new_options += ["--indep-pairwise"] + options.indep_pairwise
    if options.sge:
        new_options.append("--sge")
        new_options += ["--line-per-file-for-sge",
                        str(options.line_per_file_for_sge)]
        if options.ibs_sge_walltime is not None:
            new_options += ["--sge-walltime", options.ibs_sge_walltime]
        if options.ibs_sge_nodes is not None:
            new_options += ["--sge-nodes"] + map(str, options.ibs_sge_nodes)

    # Checking the input file
    if not allFileExists([inputPrefix + i for i in [".bed", ".bim", ".fam"]]):
        msg = "{}: not a valid binary prefix".format(inputPrefix)
        raise ProgramError(msg)

    newBfile = None
    try:
        newBfile = Relatedness.main(new_options)
    except Relatedness.ProgramError as e:
        msg = "compute genome: {}".format(e)
        raise ProgramError(msg)

    return newBfile


def allFileExists(fileList):
    """Check that all file exists.

    :param fileList: the list of file to check.

    :type fileList: list

    Check if all the files in ``fileList`` exists.

    """
    allExists = True
    for fileName in fileList:
        allExists = allExists and os.path.isfile(fileName)
    return allExists


def flipSNPs(inPrefix, outPrefix, flipFileName):
    """Flip SNPs using Plink.

    :param inPrefix: the prefix of the input file.
    :param outPrefix: the prefix of the output file.
    :param flipFileName: the name of the file containing the markers to flip.

    :type inPrefix: str
    :type outPrefix: str
    :type flipFileName: str

    Using Plink, flip a set of markers in ``inPrefix``, and saves the results
    in ``outPrefix``. The list of markers to be flipped is in ``flipFileName``.

    """
    plinkCommand = ["plink", "--noweb", "--bfile", inPrefix, "--flip",
                    flipFileName, "--make-bed", "--out", outPrefix]
    runCommand(plinkCommand)


def excludeSNPs(inPrefix, outPrefix, exclusionFileName):
    """Exclude some SNPs using Plink.

    :param inPrefix: the prefix of the input file.
    :param outPrefix: the prefix of the output file.
    :param exclusionFileName: the name of the file containing the markers to be
                              excluded.

    :type inPrefix: str
    :type outPrefix: str
    :type exclusionFileName: str

    Using Plink, exclude a list of markers from ``inPrefix``, and saves the
    results in ``outPrefix``. The list of markers are in ``exclusionFileName``.

    """
    plinkCommand = ["plink", "--noweb", "--bfile", inPrefix, "--exclude",
                    exclusionFileName, "--make-bed", "--out", outPrefix]
    runCommand(plinkCommand)


def renameSNPs(inPrefix, updateFileName, outPrefix):
    """Updates the name of the SNPs using Plink.

    :param inPrefix: the prefix of the input file.
    :param updateFileName: the name of the file containing the updated marker
                           names.
    :param outPrefix: the prefix of the output file.

    :type inPrefix: str
    :type updateFileName: str
    :type outPrefix: str

    Using Plink, changes the name of the markers in ``inPrefix`` using
    ``updateFileName``. It saves the results in ``outPrefix``.

    """
    plinkCommand = ["plink", "--noweb", "--bfile", inPrefix, "--update-map",
                    updateFileName, "--update-name", "--make-bed", "--out",
                    outPrefix]
    runCommand(plinkCommand)


def findFlippedSNPs(frqFile1, frqFile2, outPrefix):
    """Find flipped SNPs and flip them in the data.

    :param frqFile1: the name of the first frequency file.
    :param frqFile2: the name of the second frequency file.
    :param outPrefix: the prefix of the output files.

    :type frqFile1: str
    :type frqFile2: str
    :type outPrefix: str

    By reading two frequency files (``frqFile1`` and ``frqFile2``), it finds a
    list of markers that need to be flipped so that the first file becomes
    comparable with the second one. Also finds marker that need to be removed.

    A marker needs to be flipped in one of the two data set if the two markers
    are not comparable (same minor allele), but become comparable if we flip
    one of them.

    A marker will be removed if it is all homozygous in at least one data set.
    It will also be removed if it's impossible to determine the phase of the
    marker (*e.g.* if the two alleles are ``A`` and ``T`` or ``C`` and ``G``).

    """
    allelesFiles = [{}, {}]
    files = [frqFile1, frqFile2]
    for k, fileName in enumerate(files):
        try:
            with open(fileName, "r") as inputFile:
                headerIndex = None
                for i, line in enumerate(inputFile):
                    row = createRowFromPlinkSpacedOutput(line)

                    if i == 0:
                        # This is the header
                        headerIndex = dict([
                            (row[j], j) for j in xrange(len(row))
                        ])

                        # Checking the columns
                        for columnName in ["SNP", "A1", "A2"]:
                            if columnName not in headerIndex:
                                msg = "%(fileName)s: no column named " \
                                      "%(columnName)s" % locals()
                                raise ProgramError(msg)
                    else:
                        snpName = row[headerIndex["SNP"]]
                        allele1 = row[headerIndex["A1"]]
                        allele2 = row[headerIndex["A2"]]

                        alleles = set([allele1, allele2])
                        if "0" in alleles:
                            alleles.remove("0")

                        allelesFiles[k][snpName] = alleles
        except IOError:
            msg = "%(fileName)s: no such file" % locals()
            raise ProgramError(msg)

    allelesFile1, allelesFile2 = allelesFiles

    # Finding the SNPs to flip
    toFlipOutputFile = None
    try:
        toFlipOutputFile = open(outPrefix + ".snp_to_flip_in_reference", "w")
    except IOError:
        msg = "%(outPrefix)s.snp_to_flip_in_reference: can't write " \
              "file" % locals()
        raise ProgramError(msg)

    toRemoveOutputFile = None
    try:
        toRemoveOutputFile = open(outPrefix + ".snp_to_remove", "w")
    except IOError:
        msg = "%(outPrefix)s.snp_to_remove: can't write file" % locals()
        raise ProgramError(msg)

    for snpName in allelesFile1.iterkeys():
        alleles1 = allelesFile1[snpName]
        alleles2 = allelesFile2[snpName]

        if (len(alleles1) == 2) and (len(alleles2) == 2):
            # Both are heterozygous
            if ({"A", "T"} == alleles1) or ({"C", "G"} == alleles1) or \
                    ({"A", "T"} == alleles2) or ({"C", "G"} == alleles2):
                # We can't flip those..., so we remove them
                print >>toRemoveOutputFile, snpName
            else:
                if alleles1 != alleles2:
                    # Let's try the flip one
                    if flipGenotype(alleles1) == alleles2:
                        # We need to flip it
                        print >>toFlipOutputFile, snpName
                    else:
                        # Those SNP are discordant...
                        print >>toRemoveOutputFile, snpName
        else:
            # We want to remove this SNP, because there is at least one
            # homozygous individual
            print >>toRemoveOutputFile, snpName

    # Closing output files
    toFlipOutputFile.close()
    toRemoveOutputFile.close()


def computeFrequency(prefix, outPrefix):
    """Compute the frequency using Plink.

    :param prefix: the prefix of the file binary file for which we need to
                   compute frequencies.
    :param outPrefix: the prefix of the output files.

    :type prefix: str
    :type outPrefix: str

    Uses Plink to compute the frequency of all the markers in the ``prefix``
    binary file.

    """
    plinkCommand = ["plink", "--noweb", "--bfile", prefix, "--freq", "--out",
                    outPrefix]
    runCommand(plinkCommand)


def combinePlinkBinaryFiles(prefixes, outPrefix):
    """Combine Plink binary files.

    :param prefixes: a list of the prefix of the files that need to be
                     combined.
    :param outPrefix: the prefix of the output file (the combined file).

    :type prefixes: list
    :type outPrefix: str

    It uses Plink to merge a list of binary files (which is a list of prefixes
    (strings)), and create the final data set which as ``outPrefix`` as the
    prefix.

    """
    # The first file is the bfile, the others are the ones to merge
    outputFile = None
    try:
        outputFile = open(outPrefix + ".files_to_merge", "w")
    except IOError:
        msg = "%(outPrefix)s.filesToMerge: can't write file" % locals()
        raise ProgramError(msg)

    for prefix in prefixes[1:]:
        print >>outputFile, " ".join([
            prefix + i for i in [".bed", ".bim", ".fam"]
        ])

    # Closing the output files
    outputFile.close()

    # Runing plink
    plinkCommand = ["plink", "--noweb", "--bfile", prefixes[0],
                    "--merge-list", outPrefix + ".files_to_merge",
                    "--make-bed", "--out", outPrefix]
    runCommand(plinkCommand)


def findOverlappingSNPsWithReference(prefix, referencePrefixes,
                                     referencePopulations, outPrefix):
    """Find the overlapping SNPs in 4 different data sets.

    :param prefix: the prefix of all the files.
    :param referencePrefixes: the prefix of the reference population files.
    :param referencePopulations: the name of the reference population (same
                                 order as ``referencePrefixes``)
    :param outPrefix: the prefix of the output files.

    :type prefix: str
    :type referencePrefixes: list
    :type referencePopulations: list
    :type outPrefix: str

    It starts by reading the ``bim`` file of the source data set
    (``prefix.bim``). It finds all the markers (excluding the duplicated ones).
    Then it reads all of the reference population ``bim`` files
    (``referencePrefixes.bim``) and find all the markers that were found in the
    source data set.

    It creates three output files:

    * ``outPrefix.ref_snp_to_extract``: the name of the markers that needs to
      be extracted from the three reference panels.
    * ``outPrefix.source_snp_to_extract``: the name of the markers that needs
      to be extracted from the source panel.
    * ``outPrefix.update_names``: a file (readable by Plink) that will help in
      changing the names of the selected markers in the reference panels, so
      that they become comparable with the source panel.

    """
    # Reading the main file
    sourceSnpToExtract = {}
    duplicates = set()
    try:
        with open(prefix + ".bim", "r") as inputFile:
            for line in inputFile:
                row = line.rstrip("\r\n").split("\t")
                chromosome = row[0]
                position = row[3]
                snpName = row[1]

                if (chromosome, position) not in sourceSnpToExtract:
                    sourceSnpToExtract[(chromosome, position)] = snpName
                else:
                    # It's a duplicate
                    duplicates.add((chromosome, position))
    except IOError:
        msg = "%s.bim: no such file" % prefix
        raise ProgramError(msg)

    # Removing duplicates from the list
    for snpID in duplicates:
        del sourceSnpToExtract[snpID]

    # Reading each of the reference markers
    refSnpToExtract = defaultdict(dict)
    refSnp = defaultdict(set)
    for refPop, refPrefix in izip(referencePopulations, referencePrefixes):
        try:
            with open(refPrefix + ".bim", "r") as inputFile:
                for line in inputFile:
                    row = line.rstrip("\r\n").split("\t")
                    chromosome = row[0]
                    position = row[3]
                    snpName = row[1]

                    key = (chromosome, position)
                    if key in sourceSnpToExtract:
                        # We want this SNP
                        refSnpToExtract[refPop][key] = snpName
                        refSnp[refPop].add(key)
                logger.info("  - {:,d} overlaps with {}".format(
                    len(refSnp[refPop]),
                    refPrefix,
                ))

        except IOError:
            msg = "%(refPrefix)s.bim: no such file" % locals()
            raise ProgramError(msg)

    # Creating the intersect of the reference SNP
    refSnpIntersect = refSnp[referencePopulations[0]]
    for refPop in referencePopulations[1:]:
        refSnpIntersect &= refSnp[refPop]
    logger.info("  - {:,d} in common between reference "
                "panels".format(len(refSnpIntersect)))

    # Printing the names of the SNPs to extract in the reference
    refOutputFiles = {}
    try:
        for refPop in referencePopulations:
            refOutputFiles[refPop] = open(
                outPrefix + ".{}_snp_to_extract".format(refPop),
                mode="w",
            )
    except IOError:
        msg = "{}.POP_snp_to_extract: can't write file".format(outPrefix)
        raise ProgramError(msg)

    sourceOutputFile = None
    try:
        sourceOutputFile = open(outPrefix + ".source_snp_to_extract", "w")
    except IOError:
        msg = "{}.source_snp_to_extract: can't write file".format(outPrefix)
        raise ProgramError(msg)

    changeNameOutputFiles = {}
    try:
        for refPop in referencePopulations:
            changeNameOutputFiles[refPop] = open(
                outPrefix + ".{}_update_names".format(refPop),
                mode="w",
            )
    except IOError:
        msg = "%(outPrefix)s.updateNames: can't write file" % locals()
        raise ProgramError(msg)

    # Writing the file containing the SNPs to extract in the source
    for snpID in refSnpIntersect:
        print >>sourceOutputFile, sourceSnpToExtract[snpID]
        for refPop in referencePopulations:
            print >>refOutputFiles[refPop], refSnpToExtract[refPop][snpID]
            print >>changeNameOutputFiles[refPop], "\t".join([
                refSnpToExtract[refPop][snpID],
                sourceSnpToExtract[snpID],
            ])

    # Closing the output file
    sourceOutputFile.close()
    for refOutputFile in refOutputFiles.itervalues():
        refOutputFile.close()
    for changeNameOutputFile in changeNameOutputFiles.itervalues():
        changeNameOutputFile.close()


def extractSNPs(snpToExtractFileNames, referencePrefixes, popNames, outPrefix,
                runSGE, options):
    """Extract a list of SNPs using Plink.

    :param snpToExtractFileNames: the name of the files which contains the
                                  markers to extract from the original data
                                  set.
    :param referencePrefixes: a list containing the three reference population
                              prefixes (the original data sets).
    :param popNames: a list containing the three reference population names.
    :param outPrefix: the prefix of the output file.
    :param runSGE: Whether using SGE or not.
    :param options: the options.

    :type snpToExtractFileNames: list
    :type referencePrefixes: list
    :type popNames: list
    :type outPrefix: str
    :type runSGE: boolean
    :type options: argparse.Namespace

    Using Plink, extract a set of markers from a list of prefixes.

    """
    s = None
    jobIDs = []
    jobTemplates = []
    if runSGE:
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

    zipped = izip(popNames, referencePrefixes, snpToExtractFileNames)
    for popName, refPrefix, snpToExtractFileName in zipped:
        plinkCommand = ["plink", "--noweb", "--bfile", refPrefix, "--extract",
                        snpToExtractFileName, "--make-bed", "--out",
                        outPrefix + "." + popName]

        if runSGE:
            # We run using SGE
            # Creating the job template
            jt = s.createJobTemplate()
            jt.remoteCommand = plinkCommand[0]
            jt.workingDirectory = os.getcwd()
            jt.jobEnvironment = os.environ
            jt.args = plinkCommand[1:]
            jt.jobName = "_plink_extract_snps"

            # Cluster specifics
            if options.sge_walltime is not None:
                jt.hardWallclockTimeLimit = options.sge_walltime
            if options.sge_nodes is not None:
                native_spec = "-l nodes={}:ppn={}".format(options.sge_nodes[0],
                                                          options.sge_nodes[1])
                jt.nativeSpecification = native_spec

            # Running the job
            jobID = s.runJob(jt)

            # Storing the job template and the job ID
            jobTemplates.append(jt)
            jobIDs.append(jobID)

        else:
            # We run normal
            runCommand(plinkCommand)

    if runSGE:
        # We wait for all the jobs to be over
        hadProblems = []
        for jobID in jobIDs:
            retVal = s.wait(jobID, drmaa.Session.TIMEOUT_WAIT_FOREVER)
            hadProblems.append(retVal.exitStatus == 0)

        # The jobs should be finished, so we clean everything
        # Deleting the job template, and exiting the session
        for jt in jobTemplates:
            s.deleteJobTemplate(jt)

        # Closing the connection
        s.exit()

        for hadProblem in hadProblems:
            if not hadProblem:
                msg = "Some SGE jobs had errors..."
                raise ProgramError(msg)


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
        output = subprocess.check_output(
            command,
            stderr=subprocess.STDOUT,
            shell=False,
        )
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
    :py:class:`sys.stderr` and the program exists with code 1.

    """
    # Check if we have the tped and the tfam files
    prefixes_to_check = [args.bfile]
    if not args.skip_ref_pops:
        for pop in ("ceu", "yri", "jpt_chb"):
            if vars(args)["{}_bfile".format(pop)] is None:
                raise ProgramError("argument --{}-bfile is "
                                   "required".format(pop.replace("_", "-")))
        prefixes_to_check += [
            args.ceu_bfile,
            args.yri_bfile,
            args.jpt_chb_bfile,
        ]
    for prefix in prefixes_to_check:
        if prefix is None:
            msg = "no input file"
            raise ProgramError(msg)
        for fileName in [prefix + i for i in [".bed", ".bim", ".fam"]]:
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

    # Check the number of component to compute
    if args.nb_components < 2 or args.nb_components > 10:
        msg = ("nb-components: must be between 2 and 10 (inclusive), "
               "not {}".format(args.nb_components))
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
              Options            Type               Description
    =========================== ====== ========================================
    ``--bfile``                 string The input file prefix (Plink binary
                                       file).
    ``--skip-ref-pops``         bool   Perform the MDS computation, but skip
                                       the three reference panels.
    ``--ceu-bfile``             string The input file prefix for the CEU
                                       population (Plink binary file).
    ``--yri-bfile``             string The input file prefix for the YRI
                                       population (Plink binary file).
    ``--jpt-chb-bfile``         string The input file prefix for the JPT-CHB
                                       population (Plink binary file).
    ``--min-nb-snp``            int    The minimum number of markers needed to
                                       compute IBS.
    ``--indep-pairwise``        string Three numbers: window size, window shift
                                       and the r2 threshold.
    ``--maf``                   string Restrict to SNPs with MAF >= threshold.
    ``--sge``                   bool   Use SGE for parallelization.
    ``--sge-walltime``          int    The time limit (for clusters).
    ``--sge-nodes``             int    Two INTs (number of nodes and number of
                                int    processor per nodes).
    ``--ibs-sge-walltime``      int    The time limit (for clusters) (for IBS)
    ``--ibs-sge-nodes``         int    Two INTs (number of nodes and number of
                                int    processor per nodes) (for IBS).
    ``--line-per-file-for-sge`` int    The number of line per file for SGE task
                                       array.
    ``--nb-components``         int    The number of component to compute.
    ``--outliers-of``           string Finds the ouliers of this population.
    ``--multiplier``            float  To find the outliers, we look for more
                                       than x times the cluster standard
                                       deviation.
    ``--xaxis``                 string The component to use for the X axis.
    ``--yaxis``                 string The component to use for the Y axis.
    ``--format``                string The output file format.
    ``--title``                 string The title of the MDS plot.
    ``--xlabel``                string The label of the X axis.
    ``--ylabel``                string The label of the Y axis.
    ``--create-scree-plot``     bool   Computes Eigenvalues and creates a scree
                                       plot.
    ``--scree-plot-title``      string The main title of the scree plot
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


# The parser object
pretty_name = "Ethnicity"
desc = "Checks sample's ethnicity using reference populations and IBS."
long_desc = ("The script uses pairwise IBS matrix as a distance metric to "
             "identify cryptic relatedness among samples and sample outliers "
             "by multi-dimensional scaling (MDS).")
parser = argparse.ArgumentParser(description=desc)
parser.add_argument("-v", "--version", action="version",
                    version="pyGenClean version {}".format(__version__))

# The INPUT files
group = parser.add_argument_group("Input File")
group.add_argument("--bfile", type=str, metavar="FILE", required=True,
                   help=("The input file prefix (will find the plink binary "
                         "files by appending the prefix to the .bim, .bed and "
                         ".fam files, respectively."))
group.add_argument("--skip-ref-pops", action="store_true",
                   help=("Perform the MDS computation, but skip the three "
                         "reference panels."))
group.add_argument("--ceu-bfile", type=str, metavar="FILE",
                   help=("The input file prefix (will find the plink binary "
                         "files by appending the prefix to the .bim, .bed and "
                         ".fam files, respectively.) for the CEU population"))
group.add_argument("--yri-bfile", type=str, metavar="FILE",
                   help=("The input file prefix (will find the plink binary "
                         "files by appending the prefix to the .bim, .bed and "
                         ".fam files, respectively.) for the YRI population"))
group.add_argument("--jpt-chb-bfile", type=str, metavar="FILE",
                   help=("The input file prefix (will find the plink binary "
                         "files by appending the prefix to the .bim, .bed and "
                         ".fam files, respectively.) for the JPT-CHB "
                         "population"))
# The options
group = parser.add_argument_group("Options")
group.add_argument("--min-nb-snp", type=int, metavar="INT", default=8000,
                   help=("The minimum number of markers needed to compute IBS "
                         "values. [Default: %(default)d]"))
group.add_argument("--indep-pairwise", type=str, metavar="STR",
                   nargs=3, default=["50", "5", "0.1"],
                   help=("Three numbers: window size, window shift and "
                         "the r2 threshold. [default: %(default)s]"))
group.add_argument("--maf", type=str, metavar="FLOAT", default="0.05",
                   help=("Restrict to SNPs with MAF >= threshold. "
                         "[default: %(default)s]"))
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
group.add_argument("--ibs-sge-walltime", type=str, metavar="TIME",
                   help=("The walltime for the IBS jobs to run on the "
                         "cluster. Do not use if you are not required to "
                         "specify a walltime for your jobs on your cluster "
                         "(e.g. 'qsub -lwalltime=1:0:0' on the cluster)."))
group.add_argument("--ibs-sge-nodes", type=int, metavar="INT", nargs=2,
                   help=("The number of nodes and the number of processor per "
                         "nodes to use for the IBS jobs (e.g. 'qsub "
                         "-lnodes=X:ppn=Y' on the cluster, where X is the "
                         "number of nodes and Y is the number of processor to "
                         "use. Do not use if you are not required to specify "
                         "the number of nodes for your jobs on the cluster."))
group.add_argument("--line-per-file-for-sge", type=int, metavar="INT",
                   default=100, help=("The number of line per file for SGE "
                                      "task array for the IBS jobs. [default: "
                                      "%(default)d]"))
group.add_argument("--nb-components", type=int, metavar="INT", default=10,
                   help=("The number of component to compute. [default: "
                         "%(default)d]"))

# The outlier options
group = parser.add_argument_group("Outlier Options")
find_outliers.add_custom_options(group)

# The MDS plotting options
group = parser.add_argument_group("MDS Plot Options")
PlotMDS.addCustomOptions(group)

# The Scree Plot options
group = parser.add_argument_group("Scree Plot Options")
group.add_argument("--create-scree-plot", action="store_true",
                   help="Computes Eigenvalues and creates a scree plot.")
PlotEigenvalues.add_custom_options(group)

# The OUTPUT files
group = parser.add_argument_group("Output File")
group.add_argument("--out", type=str, metavar="FILE", default="ethnicity",
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
