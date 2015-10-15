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
import shutil
import random
import logging
import argparse
import StringIO
from collections import defaultdict

import numpy as np

from .. import __version__


logger = logging.getLogger("duplicated_samples")


def main(argString=None):
    """Check for duplicated samples in a tfam/tped file.

    :param argString: the options

    :type argString: list

    Here are the steps for the duplicated samples step.

    1.  Prints the options.
    2.  Reads the ``tfam`` file (:py:func:`readTFAM`).
    3.  Separate the duplicated samples from the unique samples
        (:py:func:`findDuplicates`).
    4.  Writes the unique samples into a file named
        ``prefix.unique_samples.tfam`` (:py:func:`printUniqueTFAM`).
    5.  Reads the ``tped`` file and write into ``prefix.unique_samples.tped``
        the pedigree file for the unique samples (:py:func:`processTPED`).
        Saves in memory the pedigree for the duplicated samples. Updates the
        indexes of the duplicated samples.
    6.  If there are no duplicated samples, simply copies the files
        ``prefix.unique_samples`` (``tped`` and ``tfam``) to
        ``prefix.final.tfam`` and ``prefix..final.tped``, respectively.
    7.  Computes the completion (for each of the duplicated samples) and the
        concordance of each sample pairs (:py:func:`computeStatistics`).
    8.  Prints statistics (concordance and completion)
        (:py:func:`printStatistics`).
    9.  We print the concordance matrix for each duplicated samples
        (:py:func:`printConcordance`).
    10. We print the ``tped`` and the ``tfam`` file for the duplicated samples
        (``prefix.duplicated_samples``)
        (:py:func:`printDuplicatedTPEDandTFAM`).
    11. Choose the best of each duplicates (to keep and to complete) according
        to completion and concordance (:py:func:`chooseBestDuplicates`).
    12. Creates a unique ``tped`` and ``tfam`` from the duplicated samples by
        completing the best chosen one with the other samples
        (:py:func:`createAndCleanTPED`).
    13. Merge the two tfiles together (``prefix.unique_samples`` and
        ``prefix.chosen_samples``) to create the final dataset
        (``prefix.final``) (:py:func:`addToTPEDandTFAM`).

    """
    # Getting and checking the options
    args = parseArgs(argString)
    checkArgs(args)

    logger.info("Options used:")
    for key, value in vars(args).iteritems():
        logger.info("  --{} {}".format(key.replace("_", "-"), value))

    # Reading the tfam file
    logger.info("Reading TFAM")
    tfam = readTFAM(args.tfile + ".tfam")

    # Find duplicated samples
    logger.info("Finding duplicated samples")
    uniqueSamples, duplicatedSamples = findDuplicates(tfam)

    # Prints the unique tfam
    logger.info("Creating TFAM for unique samples")
    printUniqueTFAM(tfam, uniqueSamples, args.out)

    # Process the TPED file
    logger.info("Reading TPED (and creating TPED for unique samples)")
    tped, tpedSamples = processTPED(uniqueSamples, duplicatedSamples,
                                    args.tfile + ".tped", args.out)

    if len(duplicatedSamples) == 0:
        logger.info("There are no duplicates in {}.tfam".format(args.tfile))
        # There are no duplicated samples
        try:
            shutil.copy(args.out + ".unique_samples.tfam",
                        args.out + ".final.tfam")
        except IOError:
            msg = "%s.unique_samples.tfam: can't copy to " \
                  "%s.final.tfam" % (args.out, args.out)
            raise ProgramError(msg)
        try:
            shutil.copy(args.out + ".unique_samples.tped",
                        args.out + ".final.tped")
        except IOError:
            msg = "%s.unique_samples.tped: can't copy to " \
                  "%s.final.tped" % (args.out, args.out)
            raise ProgramError(msg)

    else:
        # We continue
        # Compute statistics
        logger.info("Computing the completion and concordance of duplicated "
                    "samples")
        completion, concordance = computeStatistics(tped, tfam, tpedSamples,
                                                    duplicatedSamples,
                                                    args.out)

        # Print the statistics
        logger.info("Printing the statistics")
        completion_percentage = printStatistics(completion, concordance,
                                                tpedSamples, duplicatedSamples,
                                                args.out)

        # Print the concordance file
        logger.info("Printing the concordance file")
        concordance_percentage = printConcordance(concordance, args.out)

        # Print the duplicated TFAM and TPED
        logger.info("Creating TPED and TFAM for duplicated samples")
        printDuplicatedTPEDandTFAM(tped, tfam, tpedSamples, duplicatedSamples,
                                   args.out)

        # Choose the best duplicates
        logger.info("Choosing the best duplicates")
        chosenSamples, comp, conc = chooseBestDuplicates(
            tped,
            tpedSamples,
            duplicatedSamples,
            completion_percentage,
            concordance_percentage,
            args.out,
        )

        # Clean the genotype of the chosen samples
        logger.info("Cleaning and creating unique TPED and TFAM from "
                    "duplicated samples")
        newTPED, newTFAM = createAndCleanTPED(
            tped,
            tfam,
            tpedSamples,
            duplicatedSamples,
            chosenSamples,
            args.out,
            comp,
            args.sample_completion_threshold,
            conc,
            args.sample_concordance_threshold,
        )

        # Add the chosen TPED and TFAM
        logger.info("Creating final TPED and TFAM file")
        addToTPEDandTFAM(newTPED, newTFAM, args.out,
                         args.out + ".unique_samples")


def addToTPEDandTFAM(tped, tfam, prefix, toAddPrefix):
    """Append a tfile to another, creating a new one.

    :param tped: the ``tped`` that will be appended to the other one.
    :param tfam: the ``tfam`` that will be appended to the other one.
    :param prefix: the prefix of all the files.
    :param toAddPrefix: the prefix of the final file.

    :type tped: :py:class:`numpy.array`
    :type tfam: :py:class:`numpy.array`
    :type prefix: str
    :type toAddPrefix: str

    Here are the steps of this function:

    1. Writes the ``tped`` into ``prefix.chosen_samples.tped``.
    2. Writes the ``tfam`` into ``prefix.chosen_samples.tfam``.
    3. Copies the previous ``tfam`` (``toAddPrefix.tfam``) into the final
       ``tfam`` (``prefix.final.tfam``).
    4. Append the ``tfam`` to the final ``tfam`` (``prefix.final.tfam``).
    5. Reads the previous ``tped`` (``toAddPrefix.tped``) and append the new
       ``tped`` to it, writing the final one (``prefix.final.tped``).

    .. warning::
        The ``tped`` and ``tfam`` variables need to contain at least one
        sample.

    """
    # First, writing the chosen tped
    tpedFile = None
    try:
        tpedFile = open(prefix + ".chosen_samples.tped", "w")
    except IOError:
        msg = "%(prefix)s.chosen_samples.tped: can't write file" % locals()
        raise ProgramError(msg)
    for row in tped:
        print >>tpedFile, "\t".join(row)
    tpedFile.close()

    # Second, writing the chosen tfam
    tfamFile = None
    try:
        tfamFile = open(prefix + ".chosen_samples.tfam", "w")
    except IOError:
        msg = "%(prefix)s.chosen_samples.tfam: can't write file" % locals()
        raise ProgramError(msg)
    for row in tfam:
        print >>tfamFile, "\t".join(row)
    tfamFile.close()

    # The final TFAM file (copying and appending)
    try:
        shutil.copyfile(toAddPrefix + ".tfam", prefix + ".final.tfam")
    except IOError:
        msg = "%(toAddPrefix)s.tfam: can't copy to " \
              "%(prefix)s.final.tfam" % locals()
        raise ProgramError(msg)
    newTFAM = None
    try:
        newTFAM = open(prefix + ".final.tfam", "a")
    except IOError:
        msg = "%(prefix)s.final.tfam: can't append" % locals()
        raise ProgramError(msg)
    for row in tfam:
        print >>newTFAM, "\t".join(row)
    newTFAM.close()

    # The final TPED file
    newTPED = None
    try:
        newTPED = open(prefix + ".final.tped", "w")
    except IOError:
        msg = "%(prefix)s.final.tped: can't write file" % locals()
        raise ProgramError(msg)
    try:
        with open(toAddPrefix + ".tped", "r") as inputFile:
            for i, line in enumerate(inputFile):
                row = line.rstrip("\r\n").split("\t")
                originalSNP = row[1]
                toAdd = tped[i]
                if originalSNP != toAdd[1]:
                    msg = "%(toAddPrefix)s.tped: not good sort" % locals()
                    raise ProgramError(msg)
                row.extend(list(toAdd[4:]))
                print >>newTPED, "\t".join(row)
    except IOError:
        msg = "%(toAddPrefix)s.tped: no such file" % locals()
        raise ProgramError(msg)
    newTPED.close()


def createAndCleanTPED(tped, tfam, samples, oldSamples, chosenSamples,
                       prefix, completion, completionT, concordance,
                       concordanceT):
    """Complete a TPED for duplicate samples.

    :param tped: the ``tped`` containing the duplicated samples.
    :param tfam: the ``tfam`` containing the duplicated samples.
    :param samples: the updated position of the samples in the ``tped``
                    containing only duplicated samples.
    :param oldSamples: the original duplicated sample positions.
    :param chosenSamples: the position of the chosen samples.
    :param prefix: the prefix of all the files.
    :param completion: the completion of each of the duplicated samples.
    :param completionT: the completion threshold.
    :param concordance: the pairwise concordance of each of the duplicated
                        samples.
    :param concordanceT: the concordance threshold.

    :type tped: :py:class:`numpy.array`
    :type tfam: :py:class:`numpy.array`
    :type samples: dict
    :type oldSamples: dict
    :type chosenSamples: dict
    :type prefix: str
    :type completion: :py:class:`numpy.array`
    :type completionT: float
    :type concordance: dict
    :type concordanceT: float

    Using a ``tped`` containing duplicated samples, it creates a ``tped``
    containing unique samples by completing a chosen sample with the other
    replicates.

    .. note::
        A chosen sample is not completed using bad replicates (those that
        don't have a concordance or a completion higher than a certain
        threshold). The bad replicates are written in the file
        ``prefix.not_good_enough``.

    """
    zeroedOutFile = None
    try:
        zeroedOutFile = open(prefix + ".zeroed_out", "w")
    except IOError:
        msg = "%(prefix)s.zeroed_out: can't write file" % locals()
        raise ProgramError(msg)
    print >>zeroedOutFile, "\t".join(["famID", "indID", "snpID"])

    notGoodEnoughFile = None
    try:
        notGoodEnoughFile = open(prefix + ".not_good_enough", "w")
    except IOError:
        msg = "%(prefix)s.not_good_enough: can't write file" % locals()
        raise ProgramError(msg)
    print >>notGoodEnoughFile, "\t".join(["origIndex", "dupIndex", "famID",
                                          "indID", "reason"])
    notGoodEnoughSamples = defaultdict(set)

    # Split the tped in 'snpInfo' and 'genotypes'
    snpInfo = tped[:, :4]
    genotypes = tped[:, 4:]

    # Getting the completion and concordance indexes
    completionSampleID = {}
    concordanceSampleID = {}
    for sampleID, indexes in samples.iteritems():
        # Checking the completion
        completionToKeep = np.where(completion[indexes] >= completionT)[0]
        completionToRemove = np.where(completion[indexes] < completionT)[0]
        for k in completionToRemove:
            notGoodEnoughSamples[(str(oldSamples[sampleID][k]+1),
                                  str(indexes[k]+1), sampleID[0],
                                  sampleID[1])].add("completion")

        # Checking the concordance
        concordanceToKeep = np.where(concordance[sampleID] >= concordanceT)[0]
        concordanceToRemove = np.where(
            concordance[sampleID] < concordanceT
        )[0]
        for k in concordanceToRemove:
            notGoodEnoughSamples[(str(oldSamples[sampleID][k]+1),
                                  str(indexes[k]+1), sampleID[0],
                                  sampleID[1])].add("concordance")

        completionSampleID[sampleID] = completionToKeep
        concordanceSampleID[sampleID] = concordanceToKeep

    for i, curSNPgenotypes in enumerate(genotypes):
        for sampleID, indexes in samples.iteritems():
            # Getting the concordance and the completion
            concordanceToKeep = concordanceSampleID[sampleID]
            completionToKeep = completionSampleID[sampleID]

            # Getting the chosen sample index
            chosenOne = chosenSamples[sampleID]

            # Getting the duplicates' genotypes
            curGenotypes = curSNPgenotypes[indexes]

            # Subsetting for the good completion and concordance
            curGenotypes = curGenotypes[np.intersect1d(concordanceToKeep,
                                                       completionToKeep)]

            # Removing the no call from the genotypes
            cleanedCurGenotype = curGenotypes[np.where(curGenotypes != "0 0")]

            # Checking the number of unique genotypes
            uniqueGenotypes = np.unique(cleanedCurGenotype)

            # Checking the number of unique genotypes (except 0 0)
            if len(uniqueGenotypes) > 1:
                # There are more than one unique genotypes (except 0 0)
                # len = 0 means all were 0 0
                # len = 1 means all the same
                # len > 1 means more than one unique genotypes
                possibleAlleles = [
                    set() for k in xrange(len(uniqueGenotypes))
                ]
                for k, geno in enumerate(uniqueGenotypes):
                    possibleAlleles[k] |= set(geno.split(" "))
                allEqual = True
                for k in xrange(len(possibleAlleles)):
                    for l in xrange(k+1, len(possibleAlleles)):
                        if possibleAlleles[k] != possibleAlleles[l]:
                            allEqual = False
                            break
                    if not allEqual:
                        break
                if not allEqual:
                    # Changing the chosen sample's genotype to no call
                    tped[i][chosenOne+4] = "0 0"
                    print >>zeroedOutFile, "\t".join([sampleID[0], sampleID[1],
                                                      snpInfo[i][1]])
            else:
                # Do we want to complete the chosen sample's genotype???
                pass

    # Writing the not good enough file
    for key, value in notGoodEnoughSamples.iteritems():
        print >>notGoodEnoughFile, "\t".join(list(key) + [";".join(value)])

    # Closing the output files
    zeroedOutFile.close()
    notGoodEnoughFile.close()

    # Getting the indexes to keep
    sampleToKeep = np.array(chosenSamples.values()) + 4
    sampleToKeep.sort()
    toKeep = np.append(np.array(range(4)), sampleToKeep)

    # Subsetting the tped
    newTPED = tped[:, toKeep]

    # Getting the indexes of the tfam
    toKeep = []
    for sampleID, indexes in oldSamples.iteritems():
        i = samples[sampleID].index(chosenSamples[sampleID])
        toKeep.append(indexes[i])
    toKeep = np.array(toKeep)
    toKeep.sort()

    # Subsetting the tfam
    newTFAM = tfam[toKeep]

    return newTPED, newTFAM


def chooseBestDuplicates(tped, samples, oldSamples, completion,
                         concordance_all, prefix):
    """Choose the best duplicates according to the completion rate.

    :param tped: the ``tped`` containing the duplicated samples.
    :param samples: the updated position of the samples in the tped containing
                    only duplicated samples.
    :param oldSamples: the original duplicated sample positions.
    :param completion: the completion of each of the duplicated samples.
    :param concordance_all: the concordance of every duplicated samples.
    :param prefix: the prefix of all the files.

    :type tped: :py:class:`numpy.array`
    :type samples: dict
    :type oldSamples: dict
    :type completion: :py:class:`numpy.array`
    :type concordance_all: dict
    :type prefix: str

    :returns: a tuple where the first element is a list of the chosen samples'
              indexes, the second on is the completion and the last one is the
              concordance (a map).

    These are the steps to find the best duplicated sample:

    1. Sort the list of concordances.
    2. Sort the list of completions.
    3. Choose the best of the concordance and put in a set.
    4. Choose the best of the completion and put it in a set.
    5. Compute the intersection of the two sets. If there is one sample or
       more, then randomly choose one sample.
    6. If the intersection doesn't contain at least one sample, redo steps 3
       and 4, but increase the number of chosen best by one. Redo step 5 and 6
       (if required).

    The chosen samples are written in ``prefix.chosen_samples.info``. The rest
    are written in ``prefix.excluded_samples.info``.

    """
    # The output files
    chosenFile = None
    try:
        chosenFile = open(prefix + ".chosen_samples.info", "w")
    except IOError:
        msg = "%(prefix)s.chosen_samples.info: can't write file" % locals()
        raise ProgramError(msg)
    print >>chosenFile, "\t".join(["origIndex", "dupIndex", "famID", "indID"])

    excludedFile = None
    try:
        excludedFile = open(prefix + ".excluded_samples.info", "w")
    except IOError:
        msg = "%(prefix)s.excluded_samples.info: can't write file" % locals()
        raise ProgramError(msg)
    print >>excludedFile, "\t".join(["origIndex", "dupIndex", "famID",
                                     "indID"])

    # For each duplicated sample
    chosenIndexes = {}
    sampleConcordance = {}
    for sample, indexes in samples.iteritems():
        # Getting the completion for those duplicated samples
        currCompletion = completion[indexes]

        # Sorting those completion
        sortedCompletionIndexes = np.argsort(currCompletion)

        # Getting the concordance
        concordance = concordance_all[sample]
        currConcordance = [[] for i in xrange(len(indexes))]
        for i in xrange(len(indexes)):
            indexToKeep = list(set(range(len(indexes))) - set([i]))
            currConcordance[i] = np.mean(concordance[i, indexToKeep])
        currConcordance = np.array(currConcordance)
        if sample not in sampleConcordance:
            sampleConcordance[sample] = currConcordance

        # Sorting the concordance
        sortedConcordanceIndexes = np.argsort(currConcordance)

        # Trying to find the best duplicate to keep
        nbToCheck = 1
        chosenIndex = None
        while nbToCheck <= len(indexes):
            # Getting the `nbToCheck` best value (higher to lower)
            completionValue = currCompletion[
                sortedCompletionIndexes[nbToCheck*-1]
            ]
            concordanceValue = currConcordance[
                sortedConcordanceIndexes[nbToCheck*-1]
            ]

            # Getting the indexes to consider
            completionToConsider = set(
                np.where(currCompletion >= completionValue)[0]
            )
            concordanceToConsider = set(
                np.where(currConcordance >= concordanceValue)[0])

            # Getting the intersection of the indexes
            toConsider = concordanceToConsider & completionToConsider
            if len(toConsider) >= 1:
                chosenIndex = random.choice(list(toConsider))
                break
            nbToCheck += 1

        if chosenIndex is None:
            msg = "Could not choose the best sample ID for {}".format(sample)
            raise ProgramError(msg)

        # Printing the chosen samples
        print >>chosenFile, "\t".join([str(oldSamples[sample][chosenIndex]+1),
                                       str(indexes[chosenIndex]+1), sample[0],
                                       sample[1]])

        # Printing the excluded samples
        for i, index in enumerate(indexes):
            if i != chosenIndex:
                print >>excludedFile, "\t".join([str(oldSamples[sample][i]+1),
                                                 str(index+1), sample[0],
                                                 sample[1]])

        chosenIndexes[sample] = indexes[chosenIndex]

    # Closing the output files
    chosenFile.close()
    excludedFile.close()

    return chosenIndexes, completion, sampleConcordance


def printDuplicatedTPEDandTFAM(tped, tfam, samples, oldSamples, prefix):
    """Print the TPED and TFAM of the duplicated samples.

    :param tped: the ``tped`` containing duplicated samples.
    :param tfam: the ``tfam`` containing duplicated samples.
    :param samples: the updated position of the samples in the tped containing
                    only duplicated samples.
    :param oldSamples: the original duplicated sample positions.
    :param prefix: the prefix of all the files.

    :type tped: :py:class:`numpy.array`
    :type tfam: :py:class:`numpy.array`
    :type samples: dict
    :type oldSamples: dict
    :type prefix: str

    The ``tped`` and ``tfam`` files are written in
    ``prefix.duplicated_samples.tped`` and ``prefix.duplicated_samples.tfam``,
    respectively.

    """
    # Print the TPED
    outputTPED = None
    try:
        outputTPED = open(prefix + ".duplicated_samples.tped", "w")
    except IOError:
        msg = "%(prefix)s.duplicated_samples.tped: can't write " \
              "file" % locals()
        raise ProgramError(msg)
    for row in tped:
        print >>outputTPED, "\t".join(row)
    outputTPED.close()

    # Print the TFAM
    nbSamples = len(tped[0][4:])
    newTFAM = [0 for i in xrange(nbSamples)]
    for samples, indexes in samples.iteritems():
        oldIndexes = oldSamples[samples]
        for i, index in enumerate(indexes):
            oldIndex = oldIndexes[i]
            newTFAM[index] = tfam[oldIndex]
    outputTFAM = None
    try:
        outputTFAM = open(prefix + ".duplicated_samples.tfam", "w")
    except IOError:
        msg = "%(prefix)s.duplicated_samples.tfam: can't write " \
              "file" % locals()
        raise ProgramError(msg)
    for row in newTFAM:
        print >>outputTFAM, "\t".join(row)
    outputTFAM.close()


def printConcordance(concordance, prefix):
    """Print the concordance.

    :param concordance: the concordance of each sample.
    :param prefix: the prefix of all the files.

    :type concordance: dict
    :type prefix: str

    :returns: the concordance percentage (dict)

    The concordance is the number of genotypes that are equal when comparing a
    duplicated samples with another one, divided by the total number of
    genotypes (excluding genotypes that are no call [*i.e.* ``0``]). If a
    duplicated sample has 100% of no calls, the concordance will be zero.

    The file ``prefix.concordance`` will contain :math:`N \\times N` matrices
    for each set of duplicated samples.

    """
    outFile = None
    try:
        outFile = open(prefix + ".concordance", "w")
    except IOError:
        msg = "%s: can't write file" % prefix + ".concordance"
        raise ProgramError(msg)

    concordance_percentage = {}
    for key in concordance.iterkeys():
        print >>outFile, "#%s\t%s" % key

        # Doing the division
        none_zero = concordance[key][1] != 0
        true_concordance = np.zeros(np.multiply(*concordance[key][1].shape))
        true_concordance[np.ravel(none_zero)] = np.true_divide(
            concordance[key][0][none_zero],
            concordance[key][1][none_zero],
        )
        true_concordance.shape = concordance[key][1].shape
        true_concordance = np.asmatrix(true_concordance)
        concordance_percentage[key] = true_concordance

        output = StringIO.StringIO()
        np.savetxt(output, true_concordance, delimiter="\t", fmt="%.8f")
        print >>outFile, output.getvalue().rstrip("\r\n")

    outFile.close()

    return concordance_percentage


def printStatistics(completion, concordance, tpedSamples, oldSamples, prefix):
    """Print the statistics in a file.

    :param completion: the completion of each duplicated samples.
    :param concordance: the concordance of each duplicated samples.
    :param tpedSamples: the updated position of the samples in the tped
                        containing only duplicated samples.
    :param oldSamples: the original duplicated sample positions.
    :param prefix: the prefix of all the files.

    :type completion: :py:class:`numpy.array`
    :type concordance: dict
    :type tpedSamples: dict
    :type oldSamples: dict
    :type prefix: str

    :returns: the completion for each duplicated samples, as a
              :py:class:`numpy.array`.

    Prints the statistics (completion of each samples and pairwise concordance
    between duplicated samples) in a file (``prefix.summary``).

    """

    # Compute the completion percentage on none zero values
    none_zero_indexes = np.where(completion[1] != 0)
    completionPercentage = np.zeros(len(completion[0]), dtype=float)
    completionPercentage[none_zero_indexes] = np.true_divide(
        completion[0, none_zero_indexes],
        completion[1, none_zero_indexes],
    )

    # The output file containing the summary statistics (for each of the
    # duplicated samples, print the mean concordance and the completion).
    outputFile = None
    try:
        outputFile = open(prefix + ".summary", "w")
    except IOError:
        msg = "%(prefix)s.summary: can't write file" % locals()
        raise ProgramError(msg)

    print >>outputFile, "\t".join(["origIndex", "dupIndex", "famID", "indID",
                                   "% completion", "completion",
                                   "mean concordance"])
    for sampleID, indexes in tpedSamples.iteritems():
        for i, index in enumerate(indexes):
            # The indexes
            toPrint = [str(oldSamples[sampleID][i]+1), str(index+1)]
            # The samples
            toPrint.extend(list(sampleID))

            # The completion
            toPrint.append("%.8f" % completionPercentage[index])
            toPrint.append("%d/%d" % (completion[0][index],
                                      completion[1][index]))

            # The concordance (not on total values = 0)
            indexToKeep = list(set(range(len(indexes))) - set([i]))
            values = np.ravel(
                np.asarray(concordance[sampleID][0][i, indexToKeep])
            )
            total_values = np.ravel(
                np.asarray(concordance[sampleID][1][i, indexToKeep])
            )
            currConcordance = np.zeros(len(indexToKeep), dtype=float)
            none_zero_indexes = np.where(total_values != 0)
            currConcordance[none_zero_indexes] = np.true_divide(
                values[none_zero_indexes],
                total_values[none_zero_indexes],
            )
            currConcordance = np.mean(currConcordance)
            toPrint.append("%.8f" % currConcordance)
            print >>outputFile, "\t".join(toPrint)

    # Closing the output file
    outputFile.close()

    return completionPercentage


def computeStatistics(tped, tfam, samples, oldSamples, prefix):
    """Computes the completion and concordance of each samples.

    :param tped: the ``tped`` containing duplicated samples.
    :param tfam: the ``tfam`` containing duplicated samples.
    :param samples: the updated position of the samples in the tped containing
                    only duplicated samples.
    :param oldSamples: the original duplicated sample positions.
    :param prefix: the prefix of all the files.

    :type tped: :py:class:`numpy.array`
    :type tfam: :py:class:`numpy.array`
    :type samples: dict
    :type oldSamples: dict
    :type prefix: str

    :returns: a tuple containing the completion (:py:class:`numpy.array`) as
              first element, and the concordance (:py:class:`dict`) as last
              element.

    Reads the ``tped`` file and compute the completion for each duplicated
    samples and the pairwise concordance between duplicated samples.

    .. note::
        The completion and concordance computation excludes a markers if it's
        on chromosome 24 and if the sample is a female.

    .. note::
        A missing genotype is encoded by ``0``.

    .. note::
        No percentage is computed here, only the numbers. Percentages are
        computing in other functions: :py:func:`printStatistics`, for
        completion, and :py:func:`printConcordance`, for concordance.


    **Completion**

    Computes the completion of none zero values (where all genotypes of at
    least one duplicated sample are no call [*i.e.* ``0``]). The completion of
    sample :math:`i` (*i.e.* :math:`Comp_i`) is the number of genotypes
    that have a call divided by the total number of genotypes (the set
    :math:`G_i`):

    .. math::
        Comp_i = \\frac{||g \\in G_i\\textrm{ where }g \\neq 0||}{||G_i||}

    .. note::
        We consider a genotype as being missing if the sample is a male and if
        a marker on chromosome 23 or 24 is heterozygous.


    **Concordance**

    Computes the pairwise concordance between duplicated samples. For each
    marker, if both genotypes are not missing, we add one to the total number
    of compared markers. If both genotypes are the same, we add one to the
    number of concordant calls. We write the observed genotype difference in
    the file ``prefix.diff``. The concordance between sample :math:`i` and
    :math:`j` (*i.e.* :math:`Concordance_{i,j}`) is the number of genotypes
    that are equal divided by the total number of genotypes (excluding the no
    calls):

    .. math::
        Concordance_{i,j} = \\frac{
            ||g \\in G_i \\cup G_j \\textrm{ where } g_i = g_j \\neq 0||
        }{
            ||g \\in G_i \\cup G_j \\textrm{ where } g \\neq 0||
        }

    """
    # The diff file containing the genotype differences between a pair of
    # duplicated samples
    diffOutput = None
    try:
        diffOutput = open(prefix + ".diff", "w")
    except IOError:
        msg = "%(prefix)s.diff: can't write file" % locals()
        raise ProgramError(msg)
    print >>diffOutput, "\t".join(["name", "famID", "indID", "dupIndex_1",
                                   "dupIndex_2", "genotype_1", "genotype_2"])

    # The completion data type
    completion = np.array([[0 for i in xrange(len(tped[0][4:]))],
                           [0 for i in xrange(len(tped[0][4:]))]])

    # The concordance data type
    concordance = {}
    for sampleID in samples.keys():
        nbDup = len(samples[sampleID])
        concordance[sampleID] = [
            np.asmatrix(np.zeros((nbDup, nbDup), dtype=int)),
            np.asmatrix(np.zeros((nbDup, nbDup), dtype=int)),
        ]

    #####################################################################
    # Add options for concordance only for hetero SNPs... (all samples) #
    #####################################################################

    for row in tped:
        genotype = row[4:]
        chromosome = row[0]
        snpName = row[1]
        for key, indexes in samples.iteritems():
            for i, index in enumerate(indexes):
                sex = tfam[oldSamples[key][i]][4]
                genotype1 = set(genotype[index].split(" "))

                # Updating the completion
                if not ((chromosome == "24") and (sex == "2")):
                    if (sex == "1") and (chromosome in ["23", "24"]):
                        if ("0" not in genotype1) and (len(genotype1) != 2):
                            completion[0][index] += 1
                    elif "0" not in genotype1:
                        completion[0][index] += 1
                    completion[1][index] += 1

                    # Updating the concordance
                    for j in xrange(i + 1, len(indexes)):
                        genotype2 = set(genotype[indexes[j]].split(" "))

# This is for man on chr 23 and 24, if heterozygous!!! #
# if (sex == "1") and (chromosome in ["23", "24"]):
#     if ("0" not in genotype1) and ("0" not in genotype2):
#         if (len(genotype1) != 2) and (len(genotype2) != 2):
#             concordance[key][1][i,j] += 1
#             concordance[key][1][j,i] += 1
#             if genotype1 == genotype2:
#                 concordance[key][0][i,j] += 1
#                 concordance[key][0][j,i] += 1
#             else:
#                 # Print to diff file
#                 toPrint = [snpName, key[0], key[1],
#                            str(index+1),
#                            str(indexes[j]+1)]
#                 toPrint.append(" ".join(genotype1))
#                 if len(genotype1) == 1:
#                     toPrint[-1] += " %s" % toPrint[-1]
#                 toPrint.append(" ".join(genotype2))
#                 if len(genotype2) == 1:
#                     toPrint[-1] += " %s" % toPrint[-1]
#                 print >>diffOutput, "\t".join(toPrint)
#
# elif ("0" not in genotype1) and ("0" not in genotype2):

                        if ("0" not in genotype1) and ("0" not in genotype2):
                            # Both have calls
                            concordance[key][1][i, j] += 1
                            concordance[key][1][j, i] += 1
                            if genotype1 == genotype2:
                                concordance[key][0][i, j] += 1
                                concordance[key][0][j, i] += 1
                            else:
                                # Print to diff file
                                toPrint = [snpName, key[0], key[1],
                                           str(index+1), str(indexes[j]+1)]
                                toPrint.append(" ".join(genotype1))
                                if len(genotype1) == 1:
                                    toPrint[-1] += " %s" % toPrint[-1]
                                toPrint.append(" ".join(genotype2))
                                if len(genotype2) == 1:
                                    toPrint[-1] += " %s" % toPrint[-1]
                                print >>diffOutput, "\t".join(toPrint)

    diffOutput.close()

    for key in concordance.iterkeys():
        for i in range(len(concordance[key][0])):
            concordance[key][0][i, i] = 1
            concordance[key][1][i, i] = 1

    return completion, concordance


def processTPED(uniqueSamples, duplicatedSamples, fileName, prefix):
    """Process the TPED file.

    :param uniqueSamples: the position of unique samples.
    :param duplicatedSamples: the position of duplicated samples.
    :param fileName: the name of the file.
    :param prefix: the prefix of all the files.

    :type uniqueSamples: dict
    :type duplicatedSamples: collections.defaultdict
    :type fileName: str
    :type prefix: str

    :returns: a tuple containing the ``tped`` (:py:class:`numpy.array`) as
              first element, and the updated positions of the duplicated
              samples (:py:class:`dict`)

    Reads the entire ``tped`` and prints another one containing only unique
    samples (``prefix.unique_samples.tped``). It then creates a
    :py:class:`numpy.array` containing the duplicated samples.

    """
    # Getting the indexes
    uI = sorted(uniqueSamples.values())
    dI = []
    for item in duplicatedSamples.itervalues():
        dI.extend(item)
    dI.sort()

    tped = []
    outputFile = None
    try:
        outputFile = open(prefix + ".unique_samples.tped", "w")
    except IOError:
        msg = "%s: can't write to file" % prefix + ".unique_samples.tped"
        raise ProgramError
    with open(fileName, 'r') as inputFile:
        for line in inputFile:
            row = line.rstrip("\r\n").split("\t")
            snpInfo = row[:4]
            genotype = row[4:]

            # Printing the new TPED file (unique samples only)
            print >>outputFile, "\t".join(snpInfo + [genotype[i] for i in uI])

            # Saving the TPED file (duplicated samples only)
            tped.append(tuple(snpInfo + [genotype[i] for i in dI]))

    # Closing the output file
    outputFile.close()

    # Updating the positions
    updatedSamples = {}
    for sampleID in duplicatedSamples:
        positions = duplicatedSamples[sampleID]
        newPositions = range(len(positions))
        for i, pos in enumerate(positions):
            newPositions[i] = dI.index(pos)
        updatedSamples[sampleID] = newPositions

    tped = np.array(tped)

    return tped, updatedSamples


def printUniqueTFAM(tfam, samples, prefix):
    """Prints a new TFAM with only unique samples.

    :param tfam: a representation of a TFAM file.
    :param samples: the position of the samples
    :param prefix: the prefix of the output file name

    :type tfam: list
    :type samples: dict
    :type prefix: str

    """
    fileName = prefix + ".unique_samples.tfam"
    try:
        with open(fileName, "w") as outputFile:
            for i in sorted(samples.values()):
                print >>outputFile, "\t".join(tfam[i])
    except IOError:
        msg = "%(fileName)s: no such file"
        raise ProgramError(msg)


def findDuplicates(tfam):
    """Finds the duplicates in a TFAM.

    :param tfam: representation of a ``tfam`` file.
    :type tfam: list

    :returns: two :py:class:`dict`, containing unique and duplicated samples
              position.

    """
    uSamples = {}
    dSamples = defaultdict(list)
    for i, row in enumerate(tfam):
        sampleID = tuple(row[:2])
        if sampleID not in uSamples:
            # This is the first time we see this sample
            uSamples[sampleID] = i
        else:
            # We have seen this sample at least once...
            if sampleID not in dSamples:
                # This is the second time we see this sample...
                dSamples[sampleID].extend([uSamples[sampleID], i])
            else:
                # We have seen this sample multiple times
                dSamples[sampleID].append(i)

    # Removing the duplicates from the unique samples
    for sampleID in dSamples.iterkeys():
        if sampleID in uSamples:
            del uSamples[sampleID]

    return uSamples, dSamples


def readTFAM(fileName):
    """Reads the TFAM file.

    :param fileName: the name of the ``tfam`` file.

    :type fileName: str

    :returns: a list of tuples, representing the ``tfam`` file.

    """
    # Saving the TFAM file
    tfam = None
    with open(fileName, 'r') as inputFile:
        tfam = [
            tuple(i.rstrip("\r\n").split("\t")) for i in inputFile.readlines()
        ]

    tfam = np.array(tfam)

    return tfam


def checkArgs(args):
    """Checks the arguments and options.

    :param args: a :py:class:`argparse.Namespace` object containing the options
                 of the program.

    :type args: :py:class:`argparse.Namespace`

    :returns: ``True`` if everything was OK.

    If there is a problem with an option, an exception is raised using the
    :py:class:`ProgramError` class, a message is printed to the
    :class:`sys.stderr` and the program exists with code 1.

    """
    # Checking the input files
    for suffix in [".tped", ".tfam"]:
        fileName = args.tfile + suffix
        if not os.path.isfile(fileName):
            msg = "%(fileName)s: no such file" % locals()
            raise ProgramError(msg)

    # Checking the concordance threshold
    if ((args.sample_concordance_threshold < 0) or
            (args.sample_concordance_threshold > 1)):
        msg = "sample-concordance-threshold: must be between 0 and 1 " \
              "(not %f)" % args.sample_concordance_threshold
        raise ProgramError(msg)

    # Checking the completion threshold
    if ((args.sample_completion_threshold < 0) or
            (args.sample_completion_threshold > 1)):
        msg = "sample-completion-threshold: must be between 0 and 1 " \
              "(not %f)" % args.sample_completion_threshold
        raise ProgramError(msg)

    return True


def parseArgs(argString=None):  # pragma: no cover
    """Parses the command line options and arguments.

    :param argString: the options

    :type argString: list

    :returns: A :py:class:`argparse.Namespace` object created by
              the :py:mod:`argparse` module. It contains the values of the
              different options.

    ================================== ====== =================================
                  Options               Type             Description
    ================================== ====== =================================
    ``--tfile``                        string The input file prefix (of type
                                              ``tfile``).
    ``--sample-completion-threshold``  float  The completion threshold.
    ``--sample-concordance-threshold`` float  The concordance threshold.
    ``--out``                          string The prefix of the output files.
    ================================== ====== =================================

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
pretty_name = "Duplicated samples"
desc = "Extracts and merges duplicated samples."
long_desc = ("The script evaluates concordance and completion rate. If the "
             "thresholds are met, the script merges and completes the "
             "genotypes.")
parser = argparse.ArgumentParser(description=desc)
parser.add_argument("-v", "--version", action="version",
                    version="pyGenClean version {}".format(__version__))

# The INPUT files
group = parser.add_argument_group("Input File")
group.add_argument("--tfile", type=str, metavar="FILE", required=True,
                   help=("The input file prefix (will find the tped and tfam "
                         "file by appending the prefix to .tped and .tfam, "
                         "respectively.) The duplicated samples should have "
                         "the same identification numbers (both family and "
                         "individual ids.)"))
# The options
group = parser.add_argument_group("Options")
group.add_argument("--sample-completion-threshold", type=float,
                   metavar="FLOAT", default=0.9,
                   help=("The completion threshold to consider a replicate "
                         "when choosing the best replicates and for creating "
                         "the composite samples. [default: %(default).1f]"))
group.add_argument("--sample-concordance-threshold", type=float,
                   metavar="FLOAT", default=0.97,
                   help=("The concordance threshold to consider a replicate "
                         "when choosing the best replicates and for creating "
                         "the composite samples. [default: %(default).2f]"))
# The OUTPUT files
group = parser.add_argument_group("Output File")
group.add_argument("--out", type=str, metavar="FILE",
                   default="dup_samples",
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
