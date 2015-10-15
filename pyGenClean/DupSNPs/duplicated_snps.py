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
import math
import shutil
import random
import logging
import StringIO
import argparse
import subprocess
from collections import defaultdict

import numpy as np

from .. import __version__
from ..PlinkUtils import createRowFromPlinkSpacedOutput


logger = logging.getLogger("duplicated_snps")


def main(argString=None):
    """The main function of the module..

    Here are the steps for duplicated samples:

    1.  Prints the options.
    2.  Reads the ``map`` file to gather marker's position
        (:py:func:`readMAP`).
    3.  Reads the ``tfam`` file (:py:func:`readTFAM`).
    4.  Finds the unique markers using the ``map`` file
        (:py:func:`findUniques`).
    5.  Process the ``tped`` file. Write a file containing unique markers in
        ``prefix.unique_snps`` (``tfam`` and ``tped``). Keep in memory
        information about the duplicated markers (``tped``)
        (:py:func:`processTPED`).
    6.  If there are no duplicated markers, the file ``prefix.unique_snps``
        (``tped`` and ``tfam``) are copied to ``prefix.final``.
    7.  If there are duplicated markers, print a ``tped`` and ``tfam`` file
        containing the duplicated markers
        (:py:func:`printDuplicatedTPEDandTFAM`).
    8.  Computes the frequency of the duplicated markers (using Plink) and read
        the output file (:py:func:`computeFrequency`).
    9.  Computes the concordance and pairwise completion of each of the
        duplicated markers (:py:func:`computeStatistics`).
    10. Prints the problematic duplicated markers with a file containing the
        summary of the statistics (completion and pairwise concordance)
        (:py:func:`printProblems`).
    11. Print the pairwise concordance in a file (matrices)
        (:py:func:`printConcordance`).
    12. Choose the best duplicated markers using concordance and completion
        (:py:func:`chooseBestSnps`).
    13. Completes the chosen markers with the remaining duplicated markers
        (:py:func:`createAndCleanTPED`).
    14. Creates the final ``tped`` file, containing the unique markers, the
        chosen duplicated markers that were completed, and the problematic
        duplicated markers (for further analysis). This set excludes markers
        that were used for completing the chosen ones
        (:py:func:`createFinalTPEDandTFAM`).

    """
    # Getting and checking the options
    args = parseArgs(argString)
    checkArgs(args)

    logger.info("Options used:")
    for key, value in vars(args).iteritems():
        logger.info("  --{} {}".format(key.replace("_", "-"), value))

    # Reading the map file
    logger.info("Reading MAP file")
    mapF = readMAP(args.tfile + ".map", args.out)

    # Reading the tfam file
    logger.info("Reading TFAM file")
    tfam = readTFAM(args.tfile + ".tfam")

    # Find unique snps
    logger.info("Finding unique SNPs")
    uniqueSNPs = findUniques(mapF)

    # Process the TPED file
    logger.info("Reading TPED file")
    tped, duplicatedSNPs = processTPED(uniqueSNPs, mapF,
                                       args.tfile + ".tped",
                                       args.tfile + ".tfam", args.out)

    if len(tped) == 0:
        logger.info("There are no duplicated SNPs")
        logger.info("  - Creating final TFAM")
        # Copying the files
        # The TFAM
        try:
            shutil.copy(args.out + ".unique_snps.tfam",
                        args.out + ".final.tfam")
        except IOError:
            msg = "%s.unique_snps.tfam: can't copy file to " \
                  "%s.final.tfam" % (args.out, args.out)
            raise ProgramError(msg)

        # The TPED
        logger.info("  - Creating final TPED")
        try:
            shutil.copy(args.out + ".unique_snps.tped",
                        args.out + ".final.tped")
        except IOError:
            msg = "%s.unique_snps.tped: can't copy file to " \
                  "%s.final.tped" % (args.out, args.out)
            raise ProgramError(msg)

    else:
        # We print the TPED and TFAM for the duplicated SNPs
        logger.info("Printing duplicated SNPs TPED and TFAM files")
        printDuplicatedTPEDandTFAM(tped, args.tfile + ".tfam", args.out)

        # Computing the frequency of the duplicated SNPs
        logger.info("Computing duplicated SNPs' frequencies")
        dupSNPsFreq = computeFrequency(args.out + ".duplicated_snps",
                                       args.out)

        # Compute statistics
        logger.info("Computing concordance and completion of duplicated SNPs")
        completion, concordance = computeStatistics(tped, tfam, duplicatedSNPs)

        # Print the statistics
        logger.info("Printing duplicated SNPs summary file and finding errors "
                    "within duplicates")
        snpsToComplete = printProblems(completion, concordance, tped,
                                       duplicatedSNPs, dupSNPsFreq, args.out,
                                       args.frequency_difference)

        # Print the concordance file
        logger.info("Printing concordance file")
        printConcordance(concordance, args.out, tped, duplicatedSNPs)

        # Choose the best SNP
        logger.info("Choosing best SNP for each duplicates")
        chosenSNPs, comp, conc = chooseBestSnps(tped, duplicatedSNPs,
                                                completion, concordance,
                                                args.out)

        # Complete the SNPs
        logger.info("Completing chosen duplicates (removing discordant "
                    "genotypes)")
        newTPED, snpToRemove = createAndCleanTPED(
            tped,
            tfam,
            duplicatedSNPs,
            args.out,
            chosenSNPs,
            comp,
            conc,
            snpsToComplete,
            args.tfile + ".tfam",
            args.snp_completion_threshold,
            args.snp_concordance_threshold,
        )

        # Creates the final tped
        logger.info("Writing final TPED and TFAM file")
        createFinalTPEDandTFAM(newTPED, args.out + ".unique_snps",
                               args.out, snpToRemove)


def createFinalTPEDandTFAM(tped, toReadPrefix, prefix, snpToRemove):
    """Creates the final TPED and TFAM.

    :param tped: a representation of the ``tped`` of duplicated markers.
    :param toReadPrefix: the prefix of the unique files.
    :param prefix: the prefix of the output files.
    :param snpToRemove: the markers to remove.

    :type tped: numpy.array
    :type toReadPrefix: str
    :type prefix: str
    :type snpToRemove: set

    Starts by copying the unique markers' ``tfam`` file to
    ``prefix.final.tfam``. Then, it copies the unique markers' ``tped`` file,
    in which the chosen markers will be appended.

    The final data set will include the unique markers, the chosen markers
    which were completed, and the problematic duplicated markers (for further
    analysis). The markers that were used to complete the chosen ones are not
    present in the final data set.

    """
    # First, copying the tfam
    try:
        shutil.copy(toReadPrefix + ".tfam", prefix + ".final.tfam")
    except IOError:
        msg = "%(toReadPrefix)s.tfam: can't copy file to " \
              "%(prefix)s.final.tfam" % locals()
        raise ProgramError(msg)

    # Next, copy the tped, and append at the end
    try:
        shutil.copy(toReadPrefix + ".tped", prefix + ".final.tped")
    except IOError:
        msg = "%(toReadPrefix)s.tped: can't copy fil to " \
              "%(prefix)s.final.tped" % locals()
        raise ProgramError(msg)
    tpedFile = None
    try:
        tpedFile = open(prefix + ".final.tped", "a")
    except IOError:
        msg = "%(prefix)s.final.tped: can't append to file" % locals()
        raise ProgramError(msg)
    for i, row in enumerate(tped):
        if i not in snpToRemove:
            print >>tpedFile, "\t".join(row)
    tpedFile.close()


def createAndCleanTPED(tped, tfam, snps, prefix, chosenSNPs, completion,
                       concordance, snpsToComplete, tfamFileName, completionT,
                       concordanceT):
    """Complete a TPED for duplicated SNPs.

    :param tped: a representation of the ``tped`` of duplicated markers.
    :param tfam: a representation of the ``tfam``.
    :param snps: the position of duplicated markers in the ``tped``.
    :param prefix: the prefix of the output files.
    :param chosenSNPs: the markers that were chosen for completion (including
                       problems).
    :param completion: the completion of each of the duplicated markers.
    :param concordance: the pairwise concordance of the duplicated markers.
    :param snpsToComplete: the markers that will be completed (excluding
                           problems).
    :param tfamFileName: the name of the original ``tfam`` file.
    :param completionT: the completion threshold.
    :param concordanceT: the concordance threshold.

    :type tped: numpy.array
    :type tfam: list
    :type snps: dict
    :type prefix: str
    :type chosenSNPs: dict
    :type completion: numpy.array
    :type concordance: dict
    :type snpsToComplete: set
    :type tfamFileName: str
    :type completionT: float
    :type concordanceT: float

    :returns: a tuple containing the new ``tped`` after completion
              (:py:class:`numpy.array` as the first element, and the index of
              the markers that will need to be rid of (:py:class:`set`) as the
              last element.

    It creates three different files:

    * ``prefix.zeroed_out``: contains information about markers and samples
                             where the genotyped was zeroed out.
    * ``prefix.not_good_enough``: contains information about markers that were
                                  not good enough to help in completing the
                                  chosen markers (because of concordance or
                                  completion).
    * ``prefix.removed_duplicates``: the list of markers that where used for
                                     completing the chosen one, hence they will
                                     be removed from the final data
                                     set.

    Cycling through every genotypes of every samples of every duplicated
    markers, checks if the genotypes are all the same. If the chosen one was
    not called, but the other ones were, then we complete the chosen one with
    the genotypes for the others (assuming that they are all the same). If
    there is a difference between the genotypes, it is zeroed out for the
    chosen marker.

    """
    zeroedOutFile = None
    try:
        zeroedOutFile = open(prefix + ".zeroed_out", "w")
    except IOError:
        msg = "%(prefix).zeroed_out: can't write file" % locals()
        raise ProgramError(msg)
    print >>zeroedOutFile, "\t".join(["famID", "indID", "snpID"])

    notGoodEnoughFile = None
    try:
        notGoodEnoughFile = open(prefix + ".not_good_enough", "w")
    except IOError:
        msg = "%(prefix)s.not_good_enough: can't write file" % locals()
        raise ProgramError(msg)
    print >>notGoodEnoughFile, "\t".join(["name", "reason"])

    removedFile = None
    try:
        removedFile = open(prefix + ".removed_duplicates", "w")
    except IOError:
        msg = "%(prefix)s.removed_duplicates: can't write file" % locals()
        raise ProgramError(msg)

    notGoodEnoughSnps = set()

    # Split the tped in 'snpInfo' and 'genotypes'
    snpInfo = tped[:, :4]
    genotypes = tped[:, 4:]

    # The sed of index we want to get rid of at the end
    getRidOfIndex = set()

    for snpID, indexes in snps.iteritems():
        if snpID not in snpsToComplete:
            # We don't want to complete this SNP, so we continue to next SNP
            continue

        # Getting the completion
        completionToRemove = set(
            np.where(completion[indexes] < completionT)[0]
        )
        for k in completionToRemove:
            notGoodEnoughSnps.add((snpInfo[indexes][k, 1], "completion"))

        # Getting the concordance
        concordanceToRemove = set(
            np.where(concordance[snpID] < concordanceT)[0]
        )
        for k in concordanceToRemove:
            notGoodEnoughSnps.add((snpInfo[indexes][k, 1], "concordance"))

        # These will be the indexes to remove
        indexesToRemove = set()
        for index in completionToRemove | concordanceToRemove:
            indexesToRemove.add(indexes[index])

        # These are the indexes to keep
        indexesToKeep = []
        for index in indexes:
            if index not in indexesToRemove:
                indexesToKeep.append(index)

        # Getting the chosen SNP
        chosenOne = chosenSNPs[snpID]
        if chosenOne not in set(indexesToKeep):
            # The chosen SNP is not a good SNP, so we go to next SNP
            logger.warning("  - {} chosen but not good enough".format(
                snpInfo[chosenOne, 1],
            ))
            continue

        # Now cycling through the genotypes
        nbSamples = genotypes.shape[1]
        for sampleIndex in xrange(nbSamples):
            # We need to remove the no call and keep the unique genotypes
            curGenotypes = genotypes[indexesToKeep, sampleIndex]
            cleanedCurGenotypes = curGenotypes[
                np.where(curGenotypes != "0 0")
            ]
            uniqueCleanedCurGenotypes = np.unique(cleanedCurGenotypes)

            # Checking the number of unique genotypes
            toComplete = False
            if len(uniqueCleanedCurGenotypes) > 1:
                # There are more than one unique genotype (except 0 0)
                # len = 0 means all were 0 0
                # len = 1 means they are all the same
                # len > 1 means discordance (might need to flip)
                # Just need to check the order of the alleles
                possibleAlleles = [
                    set() for k in xrange(len(uniqueCleanedCurGenotypes))
                ]
                for k, geno in enumerate(uniqueCleanedCurGenotypes):
                    possibleAlleles[k] |= set(geno.split(" "))
                allEqual = True
                for k in xrange(len(possibleAlleles)):
                    for l in xrange(k+1, len(possibleAlleles)):
                        if possibleAlleles[k] != possibleAlleles[l]:
                            allEqual = False

                if not allEqual:
                    # The genotypes are not all equal, we set the chosen
                    # genotype to null (0 0)
                    tped[chosenOne, sampleIndex+4] = "0 0"
                    print >>zeroedOutFile, "\t".join([tfam[sampleIndex, 0],
                                                      tfam[sampleIndex, 1],
                                                      snpInfo[chosenOne, 1]])
                elif genotypes[chosenOne, sampleIndex] == "0 0":
                    toComplete = True
            elif ((len(uniqueCleanedCurGenotypes) == 1) and
                    (genotypes[chosenOne, sampleIndex] == "0 0")):
                toComplete = True

            if toComplete:
                # We complete the current individual
                tped[chosenOne, sampleIndex+4] = uniqueCleanedCurGenotypes[0]

        # We keep only the chose one
        for index in indexes:
            if index != chosenOne:
                getRidOfIndex.add(index)
                print >>removedFile, snpInfo[index, 1]

    # Writing the not good enough file
    for item in notGoodEnoughSnps:
        print >>notGoodEnoughFile, "\t".join(item)

    # Closing the output files
    zeroedOutFile.close()
    notGoodEnoughFile.close()

    # Printing the chosen file
    try:
        shutil.copy(tfamFileName, prefix + ".chosen_snps.tfam")
    except IOError:
        msg = "%(tfamFileName)s: can't copy file to " \
              "%(prefix)s.chosen_snps.tfam" % locals()
        raise ProgramError(msg)
    chosenFile = None
    try:
        chosenFile = open(prefix + ".chosen_snps.tped", "w")
    except IOError:
        msg = "%(prefix)s.chosen_snps.tped: can't write file" % locals()
        raise ProgramError(msg)
    for chosenOne in chosenSNPs.itervalues():
        snpID = (tped[chosenOne, 0], tped[chosenOne, 3])
        if snpID in snpsToComplete:
            print >>chosenFile, "\t".join(tped[chosenOne])
    chosenFile.close()

    return tped, getRidOfIndex


def chooseBestSnps(tped, snps, trueCompletion, trueConcordance, prefix):
    """Choose the best duplicates according to the completion and concordance.

    :param tped: a representation of the ``tped`` of duplicated markers.
    :param snps: the position of the duplicated markers in the ``tped``.
    :param trueCompletion: the completion of each markers.
    :param trueConcordance: the pairwise concordance of each markers.
    :param prefix: the prefix of the output files.

    :type tped: numpy.array
    :type snps: dict
    :type trueCompletion: numpy.array
    :type trueConcordance: dict
    :type prefix: str

    :returns: a tuple containing the chosen indexes (:py:class:`dict`) as the
              first element, the completion (:py:class:`numpy.array`) as the
              second element, and the concordance (:py:class:`dict`) as last
              element.

    It creates two output files: ``prefix.chosen_snps.info`` and
    ``prefix.not_chosen_snps.info``. The first one contains the markers that
    were chosen for completion, and the second one, the markers that weren't.

    It starts by computing the completion of each markers (dividing the number
    of calls divided by the total number of genotypes). Then, for each of the
    duplicated markers, we choose the best one according to completion and
    concordance (see explanation in
    :py:func:`DupSamples.duplicated_samples.chooseBestDuplicates` for more
    details).

    """
    # The output files
    chosenFile = None
    try:
        chosenFile = open(prefix + ".chosen_snps.info", "w")
    except IOError:
        msg = "%(prefix)s.chosen_snps.info: can't write file" % locals()
        raise ProgramError(msg)

    excludedFile = None
    try:
        excludedFile = open(prefix + ".not_chosen_snps.info", "w")
    except IOError:
        msg = "%(prefix)s.not_chosen_snps.info: can't write file" % locals()
        raise ProgramError(msg)

    # Computing the completion
    completion = np.true_divide(trueCompletion[0], trueCompletion[1])

    # For each duplicated SNPs
    chosenIndexes = {}
    snpConcordance = {}
    for snp, indexes in snps.iteritems():
        # Getting the completion for those duplicated SNPs
        currCompletion = completion[indexes]

        # Sorting those completion
        sortedCompletionInsexes = np.argsort(currCompletion)

        # Getting the concordance
        concordance = np.true_divide(trueConcordance[snp][0],
                                     trueConcordance[snp][1])

        currConcordance = [[] for i in xrange(len(indexes))]
        for i in xrange(len(indexes)):
            indexToKeep = list(set(range(len(indexes))) - set([i]))
            currConcordance[i] = np.mean(concordance[i, indexToKeep])
        currConcordance = np.array(currConcordance)
        if snp not in snpConcordance:
            snpConcordance[snp] = currConcordance

        # Sorting the concordance
        sortedConcordanceIndexes = np.argsort(currConcordance)

        # Trying to find the best duplicate to keep
        nbToCheck = 1
        chosenIndex = None
        while nbToCheck <= len(indexes):
            # Getting the `nbToCheck` best value (higher to lower)
            completionValue = currCompletion[
                sortedCompletionInsexes[nbToCheck*-1]
            ]
            concordanceValue = currConcordance[
                sortedConcordanceIndexes[nbToCheck*-1]
            ]

            # Getting the indexes to consider
            completionToConsider = set(
                np.where(currCompletion >= completionValue)[0]
            )
            concordanceToConsider = set(
                np.where(currConcordance >= concordanceValue)[0]
            )

            # Getting the intersection of the indexes
            toConsider = concordanceToConsider & completionToConsider
            if len(toConsider) >= 1:
                chosenIndex = random.choice(list(toConsider))
                break
            nbToCheck += 1

        if chosenIndex is None:
            msg = "Could not choose the best snp ID"
            raise ProgramError(msg)

        # Printing the chosen SNPs
        print >>chosenFile, tped[indexes[chosenIndex], 1]

        # Printing the excluded SNPs
        for i, index in enumerate(indexes):
            if i != chosenIndex:
                print >>excludedFile, tped[index, 1]

        chosenIndexes[snp] = indexes[chosenIndex]

    # Closing the output files
    chosenFile.close()
    excludedFile.close()

    return chosenIndexes, completion, snpConcordance


def computeFrequency(prefix, outPrefix):
    """Computes the frequency of the SNPs using Plink.

    :param prefix: the prefix of the input files.
    :param outPrefix: the prefix of the output files.

    :type prefix: str
    :type outPrefix: str

    :returns: a :py:class:`dict` containing the frequency of each marker.

    Start by computing the frequency of all markers using Plink. Then, it reads
    the output file, and saves the frequency and allele information.

    """
    # The plink command
    plinkCommand = ["plink", "--noweb", "--tfile", prefix, "--freq", "--out",
                    outPrefix + ".duplicated_snps"]
    runCommand(plinkCommand)

    # Reading the frequency file
    snpFreq = {}
    try:
        with open(outPrefix + ".duplicated_snps.frq", "r") as inputFile:
            headerIndex = None
            for i, line in enumerate(inputFile):
                row = createRowFromPlinkSpacedOutput(line)
                if i == 0:
                    # This is the header
                    headerIndex = dict([
                        (row[j], j) for j in xrange(len(row))
                    ])

                    # Checking the column titles
                    for columnTitle in ["SNP", "MAF", "A1", "A2"]:
                        if columnTitle not in headerIndex:
                            msg = "%(outPrefix)s.duplicated_snps.frq: no " \
                                  "column called %(columnTitle)s" % locals()
                            raise ProgramError(msg)
                else:
                    # This is data
                    snpName = row[headerIndex["SNP"]]
                    maf = row[headerIndex["MAF"]]
                    a1 = row[headerIndex["A1"]]
                    a2 = row[headerIndex["A2"]]
                    try:
                        if maf == "NA":
                            maf = 0.0
                        else:
                            maf = float(maf)
                    except ValueError:
                        msg = "%(outPrefix)s.duplicated_snps.frq: %(maf)s: " \
                              "not a valid MAF" % locals()
                        raise ProgramError(msg)
                    snpFreq[snpName] = (maf, (a1, a2))
    except IOError:
        msg = "%(outPrefix)s.duplicated_snps.freq: no such file" % locals()
        raise ProgramError(msg)

    return snpFreq


def runCommand(command):
    """Run the command in Plink.

    :param command: the command to run.

    :type command: list

    Tries to run a command using :py:mod:`subprocess`.

    """
    output = None
    try:
        output = subprocess.check_output(command,
                                         stderr=subprocess.STDOUT, shell=False)
    except subprocess.CalledProcessError:
        msg = "plink: couldn't run plink\n"
        msg += " ".join(command)
        raise ProgramError(msg)


def printDuplicatedTPEDandTFAM(tped, tfamFileName, outPrefix):
    """Print the duplicated SNPs TPED and TFAM.

    :param tped: a representation of the ``tped`` of duplicated markers.
    :param tfamFileName: the name of the original ``tfam`` file.
    :param outPrefix: the output prefix.

    :type tped: numpy.array
    :type tfamFileName: str
    :type outPrefix: str

    First, it copies the original ``tfam`` into
    ``outPrefix.duplicated_snps.tfam``. Then, it prints the ``tped`` of
    duplicated markers in ``outPrefix.duplicated_snps.tped``.

    """
    # Copying the tfam file
    try:
        shutil.copy(tfamFileName, outPrefix + ".duplicated_snps.tfam")
    except IOError:
        msg = "%(tfamFileName)s: can't copy file to " \
              "%(outPrefix)s.duplicated_snps.tfam" % locals()
        raise ProgramError(msg)

    # Writing the tped
    tpedFile = None
    try:
        tpedFile = open(outPrefix + ".duplicated_snps.tped", "w")
    except IOError:
        msg = "%(outPrefix)s.duplicated_snps.tped: can't write " \
              "file" % locals()
        raise ProgramError(msg)
    for row in tped:
        print >>tpedFile, "\t".join(row)
    tpedFile.close()


def printConcordance(concordance, prefix, tped, snps):
    """Print the concordance.

    :param concordance: the concordance.
    :param prefix: the prefix if the output files.
    :param tped: a representation of the ``tped`` of duplicated markers.
    :param snps: the position of the duplicated markers in the ``tped``.

    :type concordance: dict
    :type prefix: str
    :type tped: numpy.array
    :type snps: dict

    Prints the concordance in a file, in the format of a matrix. For each
    duplicated markers, the first line (starting with the `#` signs) contains
    the name of all the markers in the duplicated markers set. Then a :math:`N
    \\times N` matrix is printed to file (where :math:`N` is the number of
    markers in the duplicated marker list), containing the pairwise
    concordance.

    """
    outFile = None
    try:
        outFile = open(prefix + ".concordance", "w")
    except IOError:
        msg = "%s: can't write file" % prefix + ".concordance"
        raise ProgramError(msg)

    for snpID in concordance.iterkeys():
        print >>outFile, "#" + "\t".join(
            list(snpID) + list(tped[snps[snpID], 1])
        )

        # Doing the division
        true_concordance = np.true_divide(concordance[snpID][0],
                                          concordance[snpID][1])

        output = StringIO.StringIO()
        np.savetxt(output, true_concordance, delimiter="\t", fmt="%.8f")
        print >>outFile, output.getvalue().rstrip("\r\n")

    outFile.close()


def printProblems(completion, concordance, tped, snps, frequencies, prefix,
                  diffFreq):
    """Print the statistics.

    :param completion: the completion of each duplicated markers.
    :param concordance: the pairwise concordance between duplicated markers.
    :param tped: a representation of the ``tped`` of duplicated markers.
    :param snps: the positions of the duplicated markers in the ``tped``
    :param frequencies: the frequency of each of the duplicated markers.
    :param prefix: the prefix of the output files.
    :param diffFreq: the frequency difference threshold.

    :type completion: numpy.array
    :type concordance: dict
    :type tped: numpy.array
    :type snps: dict
    :type frequencies: dict
    :type prefix: str
    :type diffFreq: float

    :returns: a :py:class:`set` containing duplicated markers to complete.

    Creates a summary file (``prefix.summary``) containing information about
    duplicated markers: chromosome, position, name, alleles, status, completion
    percentage, completion number and mean concordance.

    The frequency and the minor allele are used to be certain that two
    duplicated markers are exactly the same marker (and not a tri-allelic one,
    for example).

    For each duplicated markers:

    1. Constructs the set of available alleles for the first marker.
    2. Constructs the set of available alleles for the second marker.
    3. If the two sets are different, but the number of alleles is the same, we
       try to flip one of the marker. If the two sets are the same, but the
       number of alleles is 1, we set the status to ``homo_flip``. If the
       markers are heterozygous, we set the status to ``flip``.
    4. If there is a difference in the number of alleles (one is homozygous,
       the other, heterozygous), and that there is on allele in common, we set
       the status to ``homo_hetero``. If there are no allele in common, we try
       to flip one. If the new sets have one allele in common, we set the
       status to ``homo_hetero_flip``.
    5. If the sets of available alleles are the same (without flip), we check
       the frequency and the minor alleles. If the minor allele is different,
       we set the status to ``diff_minor_allele``. If the difference in
       frequencies is higher than a threshold, we set the status to
       ``diff_frequency``.
    6. If all of the above fail, we set the status to ``problem``.

    Problems are written in the ``prefix.problems`` file, and contains the
    following columns: chromosome, position, name and status. This file
    contains all the markers with a status, as explained above.

    """

    completionPercentage = np.true_divide(completion[0], completion[1])

    outSummary = None
    try:
        outSummary = open(prefix + ".summary", "w")
    except IOError:
        msg = "%s: can't write file" % prefix + ".summary"
        raise ProgramError

    # Prints the header of the summary file
    print >>outSummary, "\t".join(["chr", "pos", "name", "alleles", "status",
                                   "% completion", "completion",
                                   "mean concordance"])

    # The data structure containing the problems
    problems = {}
    for snpID, indexes in snps.iteritems():
        for i, index in enumerate(indexes):
            # The SNP information (chromosome and position)
            toPrint = list(snpID)

            # The name of the SNP
            snpName = tped[index, 1]
            toPrint.append(snpName)

            # The frequency of the SNP
            snpFreq, mafAlleles = frequencies[snpName]

            # A list of the other SNP name with problems
            otherSnpNameWithProblem = set()

            # The alleles
            alleles = set()
            otherAlleles = set()
            status = []
            for genotype in np.unique(tped[index, 4:]):
                alleles |= set(genotype.split(" "))
            if "0" in alleles:
                alleles.remove("0")
            for j in xrange(i+1, len(indexes)):
                otherIndex = indexes[j]
                otherSnpName = tped[otherIndex, 1]

                # The frequency of the other SNP
                otherSnpFreq, otherMafAlleles = frequencies[otherSnpName]

                # Checking the alleles
                for genotype in np.unique(tped[otherIndex, 4:]):
                    otherAlleles |= set(genotype.split(" "))
                if "0" in otherAlleles:
                    otherAlleles.remove("0")
                if alleles != otherAlleles:
                    if len(alleles) == len(otherAlleles):
                        # Same number of alleles
                        # Try the flipped ones
                        otherAlleles = flipGenotype(otherAlleles)
                        if alleles == otherAlleles:
                            if len(alleles) == 1:
                                status.append("homo_flip")
                                otherSnpNameWithProblem.add(otherSnpName)
                            else:
                                status.append("flip")
                                otherSnpNameWithProblem.add(otherSnpName)
                        else:
                            status.append("problem")
                            otherSnpNameWithProblem.add(otherSnpName)
                    else:
                        # Different number of alleles
                        if len(alleles & otherAlleles) == 1:
                            status.append("homo_hetero")
                            otherSnpNameWithProblem.add(otherSnpName)
                        else:
                            # Try the flipped one
                            otherAlleles = flipGenotype(otherAlleles)
                            if len(alleles & otherAlleles) == 1:
                                status.append("homo_hetero_flip")
                                otherSnpNameWithProblem.add(otherSnpName)
                            else:
                                status.append("problem")
                                otherSnpNameWithProblem.add(otherSnpName)
                else:
                    # The alleles are the same, so we check the frequency
                    if mafAlleles[0] != otherMafAlleles[0]:
                        # They don't have the same minor allele
                        status.append("diff_minor_allele")
                        otherSnpNameWithProblem.add(otherSnpName)
                    elif math.fabs(snpFreq - otherSnpFreq) > diffFreq:
                        # They don't have same frequency
                        status.append("diff_frequency")
                        otherSnpNameWithProblem.add(otherSnpName)

            alleles = list(alleles)
            alleles.sort()
            if len(alleles) == 1:
                alleles.append(alleles[0])
            toPrint.append(" ".join(alleles))
            toPrint.append(";".join(status))

            # The completion
            toPrint.append("%.8f" % completionPercentage[index])
            toPrint.append("%d/%d" % (completion[0][index],
                                      completion[1][index]))

            # The concordance
            indexToKeep = list(set(range(len(indexes))) - set([i]))
            currConcordance = np.true_divide(
                concordance[snpID][0][i, indexToKeep],
                concordance[snpID][1][i, indexToKeep],
            )
            currConcordance = np.mean(currConcordance)
            toPrint.append("%.8f" % currConcordance)
            print >>outSummary, "\t".join(toPrint)

            # Now updating the problems data structure
            if len(status) != len(otherSnpNameWithProblem):
                msg = "There is a problem with the problematic SNPs"
                raise ProgramError(msg)

            if len(status) > 0:
                if snpID not in problems:
                    tmp = {"snpNames": {snpName}, "problems": set()}
                    problems[snpID] = tmp

                # We have problems
                problems[snpID]["snpNames"] |= otherSnpNameWithProblem
                problems[snpID]["problems"] |= set(status)

    outSummary.close()

    outProblems = None
    try:
        outProblems = open(prefix + ".problems", "w")
    except IOError:
        msg = "%s: can't write file" % prefix + ".problems"
        raise ProgramError

    # Printing the header of the problem file...
    print >>outProblems, "\t".join(["chr", "pos", "name", "status"])
    for snpID in problems.iterkeys():
        toPrint = list(snpID)
        toPrint.append(";".join(list(problems[snpID]["snpNames"])))
        toPrint.append(";".join(list(problems[snpID]["problems"])))
        print >>outProblems, "\t".join(toPrint)

    outProblems.close()

    # Returning the SNPs to complete
    return set(snps.keys()) - set(problems.keys())


def computeStatistics(tped, tfam, snps):
    """Computes the completion and concordance of each SNPs.

    :param tped: a representation of the ``tped``.
    :param tfam: a representation of the ``tfam``
    :param snps: the position of the duplicated markers in the ``tped``.

    :type tped: numpy.array
    :type tfam: list
    :type snps: dict

    :returns: a tuple containing the completion of duplicated markers
              (:py:class:`numpy.array`) as first element, and the concordance
              (:py:class:`dict`) of duplicated markers, as last element.

    A marker's completion is compute using this formula (where :math:`G_i` is
    the set of genotypes for the marker :math:`i`):

    .. math::
        Completion_i = \\frac{||g \\in G_i \\textrm{ where } g \\neq 0||}
                             {||G_i||}

    The pairwise concordance between duplicated markers is compute as follow
    (where :math:`G_i` and :math:`G_j` are the sets of genotypes for markers
    :math:`i` and :math:`j`, respectively):

    .. math::
        Concordance_{i,j} = \\frac{
            ||g \\in G_i \\cup G_j \\textrm{ where } g_i = g_j \\neq 0||
        }{
            ||g \\in G_i \\cup G_j \\textrm{ where } g \\neq 0||
        }

    Hence, we only computes the numerators and denominators of the completion
    and concordance, for future reference.

    .. note::
        When the genotypes are not comparable, the function tries to flip one
        of the genotype to see if it becomes comparable.

    """
    # The completion data type
    completion = np.array([[0 for i in xrange(len(tped))],
                           [0 for i in xrange(len(tped))]])

    # The concordance data type
    concordance = {}
    for snpID in snps.keys():
        nbDup = len(snps[snpID])
        concordance[snpID] = [
            np.asmatrix(np.zeros((nbDup, nbDup), dtype=int)),
            np.asmatrix(np.zeros((nbDup, nbDup), dtype=int))
        ]

    # The women and the no sex
    menIndex = np.where(tfam[:, 4] == "1")
    womenIndex = np.where(tfam[:, 4] == "2")
    noSexIndex = np.where(tfam[:, 4] == "0")

    for snpID, indexes in snps.iteritems():
        nbDup = len(indexes)
        currGenotypes = tped[indexes, 4:]
        chromosome, position = snpID

#         if chromosome == "24":
#             # Remove the heterozygous men
#             menToRemove = getIndexOfHeteroMen(currGenotypes, menIndex)
#             # Remove the women and the no sex
#             currGenotypes = np.delete(currGenotypes,
#                                        np.hstack((womenIndex, noSexIndex,
#                                                    menToRemove)), 1)
#         elif chromosome == "23":
#             # Remove the heterozygous men
#             menToRemove = getIndexOfHeteroMen(currGenotypes, menIndex)
#             # Remove the no sex
#             currGenotypes = np.delete(currGenotypes,
#                                        np.hstack((noSexIndex, menToRemove)),
#                                        1)

        for i in xrange(nbDup):
            # Compute completion here
            completion[0][indexes[i]] = len(
                np.where(currGenotypes[i] != "0 0")[0]
            )
            completion[1][indexes[i]] = len(currGenotypes[i])
            for j in xrange(i+1, nbDup):
                # Compute concordance here
                # Removing samples with at least one null genotype
                nullGenotypeIndexes = np.where(
                    np.any(currGenotypes[[i, j]] == "0 0", 0)
                )
                subGenotypes = np.delete(
                    currGenotypes,
                    nullGenotypeIndexes,
                    1,
                )

                # Finding the errors in the subseted genotypes
                errorIndexes = np.where(subGenotypes[i] != subGenotypes[j])[0]
                nbDiff = len(errorIndexes)

                for k in errorIndexes:
                    # Getting the genotypes
                    genotype1 = set(subGenotypes[i, k].split(" "))
                    genotype2 = set(subGenotypes[j, k].split(" "))

                    # Checking for flips
                    if len(genotype1) == len(genotype2):
                        # Both have the same number of different alleles,
                        # so they might be flipped
                        genotype2 = flipGenotype(genotype2)
                        if genotype1 == genotype2:
                            # The genotypes are equivalent after the flip
                            nbDiff -= 1

                # Updating the concordance
                nbTot = len(subGenotypes[i])
                concordance[snpID][0][i, j] = nbTot - nbDiff
                concordance[snpID][0][j, i] = nbTot - nbDiff
                if nbTot == 0:
                    # We will have a division by 0...
                    nbTot = 1
                concordance[snpID][1][i, j] = nbTot
                concordance[snpID][1][j, i] = nbTot

    for snpID in concordance.iterkeys():
        for i in range(len(concordance[snpID][0])):
            concordance[snpID][0][i, i] = 1
            concordance[snpID][1][i, i] = 1

    return completion, concordance


def getIndexOfHeteroMen(genotypes, menIndex):
    """Get the indexes of heterozygous men.

    :param genotypes: the genotypes of everybody.
    :param menIndex: the indexes of the men (for the genotypes).

    :type genotypes: numpy.array
    :type menIndex: numpy.array

    :returns: a :py:class:`numpy.array` containing the indexes of the genotypes
              to remove.

    Finds the mean that have a heterozygous genotype for this current marker.
    Usually used on sexual chromosomes.

    """
    toRemove = set()
    for i in menIndex[0]:
        for genotype in [set(j.split(" ")) for j in genotypes[:, i]]:
            if len(genotype) != 1:
                # We have an heterozygous
                toRemove.add(i)

    toRemove = list(toRemove)
    toRemove.sort()

    return (np.array(toRemove, dtype=int),)


def flipGenotype(genotype):
    """Flips a genotype.

    :param genotype: the genotype to flip.

    :type genotype: set

    :returns: the new flipped genotype (as a :py:class:`set`)

    .. testsetup::

        from pyGenClean.DupSNPs.duplicated_snps import flipGenotype

    .. doctest::

        >>> flipGenotype({"A", "T"})
        set(['A', 'T'])
        >>> flipGenotype({"C", "T"})
        set(['A', 'G'])
        >>> flipGenotype({"T", "G"})
        set(['A', 'C'])
        >>> flipGenotype({"0", "0"})
        Traceback (most recent call last):
            ...
        ProgramError: 0: unkown allele
        >>> flipGenotype({"A", "N"})
        Traceback (most recent call last):
            ...
        ProgramError: N: unkown allele

    """
    newGenotype = set()
    for allele in genotype:
        if allele == "A":
            newGenotype.add("T")
        elif allele == "C":
            newGenotype.add("G")
        elif allele == "T":
            newGenotype.add("A")
        elif allele == "G":
            newGenotype.add("C")
        else:
            msg = "%(allele)s: unknown allele" % locals()
            raise ProgramError(msg)

    return newGenotype


# def processTPED(uniqueSNPs, duplicatedSNPs, mapF, fileName, tfam, prefix):
def processTPED(uniqueSNPs, mapF, fileName, tfam, prefix):
    """Process the TPED file.

    :param uniqueSNPs: the unique markers.
    :param mapF: a representation of the ``map`` file.
    :param fileName: the name of the ``tped`` file.
    :param tfam: the name of the ``tfam`` file.
    :param prefix: the prefix of all the files.

    :type uniqueSNPs: dict
    :type mapF: list
    :type fileName: str
    :type tfam: str
    :type prefix: str

    :returns: a tuple with the representation of the ``tped`` file
              (:py:class:`numpy.array`) as first element, and the updated
              position of the duplicated markers in the ``tped``
              representation.

    Copies the ``tfam`` file into ``prefix.unique_snps.tfam``. While reading
    the ``tped`` file, creates a new one (``prefix.unique_snps.tped``)
    containing only unique markers.

    """
    # Copying the tfam file
    try:
        shutil.copy(tfam, prefix + ".unique_snps.tfam")
    except IOError:
        msg = "%s: can't write file" % prefix + ".unique_snps.tfam"
        raise ProgramError(msg)

    tped = []
    updatedSNPs = defaultdict(list)
    outputFile = None
    try:
        outputFile = open(prefix + ".unique_snps.tped", "w")
    except IOError:
        msg = "%s: can't write to file" % prefix + ".unique_snps.tped"
        raise ProgramError(msg)
    nbSNP = 0
    with open(fileName, 'r') as inputFile:
        for line in inputFile:
            nbSNP += 1
            row = line.rstrip("\r\n").split("\t")
            snpInfo = row[:4]
            genotype = [i.upper() for i in row[4:]]

            chromosome = snpInfo[0]
            position = snpInfo[3]

            if (chromosome, position) in uniqueSNPs:
                # Printing the new TPED file (unique SNPs only)
                print >>outputFile, "\t".join(snpInfo + genotype)
            else:
                # Saving the TPED file (duplicated samples only)
                currPos = len(tped)
                tped.append(tuple(snpInfo + genotype))
                updatedSNPs[(chromosome, position)].append(currPos)
    outputFile.close()

    if len(mapF) != nbSNP:
        msg = "%(fileName)s: no the same number of SNPs than MAP " \
              "file" % locals()
        raise ProgramError(msg)

    tped = np.array(tped)

    return tped, updatedSNPs


def findUniques(mapF):
    """Finds the unique markers in a MAP.

    :param mapF: representation of a ``map`` file.

    :type mapF: list

    :returns: a :py:class:`dict` containing unique markers (according to their
              genomic localisation).

    """
    uSNPs = {}
    dSNPs = defaultdict(list)
    for i, row in enumerate(mapF):
        chromosome = row[0]
        position = row[3]
        snpID = (chromosome, position)
        if snpID not in uSNPs:
            # This is the first time we see this sample
            uSNPs[snpID] = i
        else:
            # We have seen this sample at least once...
            if snpID not in dSNPs:
                # This is the second time we see this sample...
                dSNPs[snpID].extend([uSNPs[snpID], i])
            else:
                # We have seen this sample multiple times
                dSNPs[snpID].append(i)

    # Removing the duplicates from the unique samples
    for snpID in dSNPs.iterkeys():
        if snpID in uSNPs:
            del uSNPs[snpID]

    return uSNPs


def readTFAM(fileName):
    """Reads the TFAM file.

    :param fileName: the name of the ``tfam`` file.

    :type fileName: str

    :returns: a representation the ``tfam`` file (:py:class:`numpy.array`).

    """
    # Saving the TFAM file
    tfam = None
    with open(fileName, 'r') as inputFile:
        tfam = [
            tuple(i.rstrip("\r\n").split("\t")) for i in inputFile.readlines()
        ]

    tfam = np.array(tfam)

    return tfam


def readMAP(fileName, prefix):
    """Reads the MAP file.

    :param fileName: the name of the ``map`` file.

    :type fileName: str

    :returns: a list of tuples, representing the ``map`` file.

    While reading the ``map`` file, it saves a file
    (``prefix.duplicated_marker_names``) containing the name of the unique
    duplicated markers.

    """
    # Saving the MAP file
    mapF = None
    with open(fileName, 'r') as inputFile:
        mapF = [
            tuple(i.rstrip("\r\n").split("\t")) for i in inputFile.readlines()
        ]

    # Test for uniqueness of names
    marker_names = np.array([i[1] for i in mapF])
    nb_with_same_name = len(marker_names) - len(np.unique(marker_names))
    if nb_with_same_name > 0:
        logger.info("  - {} markers with same name".format(nb_with_same_name))
        u, indices = np.unique(marker_names, return_index=True)
        duplicated_indices = np.setdiff1d(np.arange(len(marker_names)),
                                          indices)
        duplicated_indices = np.in1d(marker_names,
                                     marker_names[duplicated_indices])
        marker_chr_pos = np.array([(i[0], i[3]) for i in mapF],
                                  dtype=[("chr", int),
                                         ("pos", int)])[duplicated_indices]
        marker_names = marker_names[duplicated_indices]
        try:
            duplicated_marker_names = open(prefix + ".duplicated_marker_names",
                                           "w")
        except IOError:
            msg = "{}.duplicated_marker_names: can't write file".format(prefix)
            raise ProgramError(msg)
        for i in xrange(len(marker_names)):
            print >>duplicated_marker_names, "\t".join([marker_names[i]] +
                                                       map(str,
                                                           marker_chr_pos[i]))
        duplicated_marker_names.close()

    return mapF


def checkArgs(args):
    """Checks the arguments and options.

    :param args: an object containing the options of the program.

    :type args: argparse.Namespace

    :returns: ``True`` if everything was OK.

    If there is a problem with an option, an exception is raised using the
    :py:class:`ProgramError` class, a message is printed to the
    :class:`sys.stderr` and the program exists with code 1.

    """
    # Checking the input files
    for suffix in [".tped", ".tfam", ".map"]:
        fileName = args.tfile + suffix
        if not os.path.isfile(fileName):
            msg = "%(fileName)s: no such file" % locals()
            raise ProgramError(msg)

    # Checking the concordance threshold
    if ((args.snp_concordance_threshold < 0) or
            (args.snp_concordance_threshold > 1)):
        msg = "snp-concordance-threshold: must be between 0 and 1 " \
              "(not %f)" % args.snp_concordance_threshold
        raise ProgramError(msg)

    # Checking the completion threshold
    if ((args.snp_completion_threshold < 0) or
            (args.snp_completion_threshold > 1)):
        msg = "snp-completion-threshold: must be between 0 and 1 " \
              "(not %f)" % args.snp_completion_threshold
        raise ProgramError(msg)

    # Checking the difference in frequency
    if (args.frequency_difference < 0) or (args.frequency_difference > 1):
        msg = ("{}: maximal frequency difference: value must be between 0 and "
               "1 (inclusively)".format())
        raise ProgramError(msg)

    return True


def parseArgs(argString=None):  # pragma: no cover
    """Parses the command line options and arguments.

    :param argString: the options.

    :type argString: list

    :returns: A :py:class:`argparse.Namespace` object created by the
              :py:mod:`argparse` module. It contains the values of the
              different options.

    =============================== ====== ====================================
                 Options             Type              Description
    =============================== ====== ====================================
    ``--tfile``                     string The input file prefix (Plink tfile).
    ``--snp-completion-threshold``  float  The completion threshold to consider
                                           a replicate when choosing the best
                                           replicates.
    ``--snp-concordance-threshold`` float  The concordance threshold to
                                           consider a replicate when choosing
                                           the best replicates.
    ``--frequency_difference``      float  The maximum difference in frequency
                                           between duplicated markers.
    ``--out``                       string The prefix of the output files.
    =============================== ====== ====================================

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
pretty_name = "Duplicated markers"
desc = "Extracts and merges duplicated markers."
long_desc = ("The script searches for duplicated markers according to "
             "chromosomal location. It then evaluates concordance, completion "
             "rate, allele calls and minor allele frequency (MAF). The script "
             "keeps markers with different allele calls or with different "
             "MAF. If thresholds are met, the script merges and completes the "
             "genotypes.")
parser = argparse.ArgumentParser(description=desc)
parser.add_argument("-v", "--version", action="version",
                    version="pyGenClean version {}".format(__version__))

# The INPUT files
group = parser.add_argument_group("Input File")
group.add_argument("--tfile", type=str, metavar="FILE", required=True,
                   help=("The input file prefix (will find the tped and tfam "
                         "file by appending the prefix to .tped and .tfam, "
                         "respectively. A .map file is also required."))
# The options
group = parser.add_argument_group("Options")
group.add_argument("--snp-completion-threshold", type=float, metavar="FLOAT",
                   default=0.9, help=("The completion threshold to consider "
                                      "a replicate when choosing the best "
                                      "replicates and for composite creation. "
                                      "[default: %(default).1f]"))
group.add_argument("--snp-concordance-threshold", type=float, metavar="FLOAT",
                   default=0.98, help=("The concordance threshold to consider "
                                       "a replicate when choosing the best "
                                       "replicates and for composite "
                                       "creation. [default: %(default).2f]"))
group.add_argument("--frequency_difference", type=float, metavar="FLOAT",
                   default=0.05, help=("The maximum difference in frequency "
                                       "between duplicated markers [default: "
                                       "%(default).2f]"))
# The OUTPUT files
group = parser.add_argument_group("Output File")
group.add_argument("--out", type=str, metavar="FILE", default="dup_snps",
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
