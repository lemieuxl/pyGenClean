#!/usr/bin/env python2.7

import os
import sys
import shutil
import random
import argparse
import StringIO
from collections import defaultdict

import numpy as npy

def main(argString=None):
    """Check for duplicated samples in a tfam/tped file.

    Here are the steps for duplicated samples.

    1.  Reads the ``tfam`` file.
    2.  Seperate the duplicated samples from the unique samples.
    3.  Writes the unique samples into a file named
        ``prefix.unique_samples.tfam``.
    4.  Reads the ``tped`` file and write into ``prefix.unique_samples.tped`` the
        pedigree file for the unique samples. Saves in memory the pedigree for
        the duplicated samples. Updates the indexes of the duplicated samples.
    5.  If there are no duplicated samples, simply copies the files
        ``prefix.unique_samples`` (``tped`` and ``tfam``) to
        ``prefix.final.tfam`` and ``prefix..final.tped``, respectively.
    6.  Computes the completion (for each of the duplicated samples) and the
        concordance of each sample pairs.
    7.  We print statistics (concordance and completion).
    8.  We print the concordance matrix for each duplicated samples.
    9.  We print the the ``tped`` and the ``tfam`` file for the duplicated
        samples (``prefix.duplicated_samples``)
    10. Choose the best of each duplicates (to keep and to complete) according
        to completion and concordance.

    """
    # Getting and checking the options
    args = parseArgs(argString)
    checkArgs(args)

    print "   - Options used:"
    for key, value in vars(args).iteritems():
        print "      --{} {}".format(key, value)

    # Reading the tfam file
    print "   - Reading TFAM"
    tfam = readTFAM(args.tfile + ".tfam")

    # Find duplicated samples
    print "   - Finding duplicated samples"
    uniqueSamples, duplicatedSamples = findDuplicates(tfam)

    # Prints the unique tfam
    print "   - Creating TFAM for unique samples"
    printUniqueTFAM(tfam, uniqueSamples, args.out)

    # Process the TPED file
    print "   - Reading TPED (and creating TPED for unique samples)"
    tped, tpedSamples = processTPED(uniqueSamples, duplicatedSamples,
                                    args.tfile + ".tped", args.out)

    if len(duplicatedSamples) == 0:
        print "   - There are no duplicates in {}.tfam".format(args.tfile)
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
        print ("   - Computing the completion and concordance of duplicated\n"
               "     samples")
        completion, concordance = computeStatistics(tped, tfam, tpedSamples,
                                                    duplicatedSamples, args.out)

        # Print the statistics
        print "   - Printing the statistics"
        completion_percentage = printStatistics(completion, concordance,
                                                tpedSamples, duplicatedSamples,
                                                args.out)

        # Print the concordance file
        print "   - Printing the concordance file"
        concordance_percentage = printConcordance(concordance, args.out)

        # Print the duplicated TFAM and TPED
        print "   - Creating TPED and TFAM for duplicated samples"
        printDuplicatedTPEDandTFAM(tped, tfam, tpedSamples, duplicatedSamples,
                                   args.out)

        # Choose the best duplicates
        print "   - Choosing the best duplicates"
        chosenSamples, comp, conc = chooseBestDuplicates(tped, tpedSamples,
                                                         duplicatedSamples,
                                                         completion_percentage,
                                                         concordance_percentage,
                                                         args.out)

        # Clean the genotype of the chosen samples
        print ("   - Cleaning and creating unique TPED and TFAM from\n"
               "     duplicated samples")
        newTPED, newTFAM = createAndCleanTPED(tped, tfam, tpedSamples,
                                              duplicatedSamples,
                                              chosenSamples, args.out, comp,
                                              args.sample_completion_threshold,
                                              conc,
                                              args.sample_concordance_threshold)

        # Add the chosen TPED and TFAM
        print "   - Creating final TPED and TFAM file"
        addToTPEDandTFAM(newTPED, newTFAM, args.out,
                         args.out + ".unique_samples")


def addToTPEDandTFAM(tped, tfam, prefix, toAddPrefix):
    """Add the tped and the tfam to the prefix.tped and prefix.tfam."""
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
    """Complete a TPED for duplicate samples."""
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
    snpInfo = tped[:,:4]
    genotypes = tped[:,4:]

    # Getting the completion and concordance indexes
    completionSampleID = {}
    concordanceSampleID = {}
    for sampleID, indexes in samples.iteritems():
        # Checking the completion
        completionToKeep = npy.where(completion[indexes] >= completionT)[0]
        completionToRemove = npy.where(completion[indexes] < completionT)[0]
        for k in completionToRemove:
            notGoodEnoughSamples[(str(oldSamples[sampleID][k]+1),
                                  str(indexes[k]+1), sampleID[0],
                                  sampleID[1])].add("completion")

        # Checking the concordance
        concordanceToKeep = npy.where(concordance[sampleID] >= concordanceT)[0]
        concordanceToRemove = npy.where(concordance[sampleID] < concordanceT)[0]
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
            curGenotypes = curGenotypes[npy.intersect1d(concordanceToKeep,
                                                        completionToKeep)]

            # Removing the no call from the genotypes
            cleanedCurGenotype = curGenotypes[npy.where(curGenotypes != "0 0")]

            # Checking the number of unique genotypes
            uniqueGenotypes = npy.unique(cleanedCurGenotype)

            # Checking the number of unique genotypes (except 0 0)
            if len(uniqueGenotypes) > 1:
                # There are more than one unique genotypes (except 0 0)
                # len = 0 means all were 0 0
                # len = 1 means all the same
                # len > 1 means more than one unique genotypes
                possibleAlleles = [set() \
                                    for k in xrange(len(uniqueGenotypes))]
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
    sampleToKeep = npy.array(chosenSamples.values()) + 4
    sampleToKeep.sort()
    toKeep = npy.append(npy.array(range(4)), sampleToKeep)

    # Subsetting the tped
    newTPED = tped[:,toKeep]

    # Getting the indexes of the tfam
    toKeep = []
    for sampleID, indexes in oldSamples.iteritems():
        i = samples[sampleID].index(chosenSamples[sampleID])
        toKeep.append(indexes[i])
    toKeep = npy.array(toKeep)
    toKeep.sort()

    # Subsetting the tfam
    newTFAM = tfam[toKeep]

    return newTPED, newTFAM


def chooseBestDuplicates(tped, samples, oldSamples, completion,
                        concordance_all, prefix):
    """Choose the best duplicates according to the completion rate.
    
    Describe what the f*** is going on here...
    
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
        sortedCompletionIndexes = npy.argsort(currCompletion)

        # Getting the concordance
        concordance = concordance_all[sample]
        currConcordance = [[] for i in xrange(len(indexes))]
        for i in xrange(len(indexes)):
            indexToKeep = list(set(range(len(indexes))) - set([i]))
            currConcordance[i] = npy.mean(concordance[i, indexToKeep])
        currConcordance = npy.array(currConcordance)
        if sample not in sampleConcordance:
            sampleConcordance[sample] = currConcordance

        # Sorting the concordance
        sortedConcordanceIndexes = npy.argsort(currConcordance)

        # Trying to find the best duplicate to keep
        nbToCheck = 1
        chosenIndex = None
        while nbToCheck <= len(indexes):
            # Getting the `nbToCheck` best value (higher to lower)
            completionValue = currCompletion[sortedCompletionIndexes[nbToCheck*-1]]
            concordanceValue = currConcordance[sortedConcordanceIndexes[nbToCheck*-1]]

            # Getting the indexes to consider
            completionToConsider = set(npy.where(currCompletion >= \
                                                    completionValue)[0])
            concordanceToConsider = set(npy.where(currConcordance >= \
                                                    concordanceValue)[0])

            # Getting the intesection of the indexes
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

        # Printing the excuded samples
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
    """Print the TPED and TFAM of the duplicated samples."""
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
    """Print the concordance."""
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
        true_concordance = npy.zeros(npy.multiply(*concordance[key][1].shape))
        true_concordance[npy.ravel(none_zero)] = npy.true_divide(concordance[key][0][none_zero],
                                                                 concordance[key][1][none_zero])
        true_concordance.shape = concordance[key][1].shape
        true_concordance = npy.asmatrix(true_concordance)
        concordance_percentage[key] = true_concordance

        output = StringIO.StringIO()
        npy.savetxt(output, true_concordance, delimiter="\t", fmt="%.8f")
        print >>outFile, output.getvalue().rstrip("\r\n")

    outFile.close()

    return concordance_percentage


def printStatistics(completion, concordance, tpedSamples, oldSamples, prefix):
    """Print the statistics."""

    # Compute the completion percentage on none zero values
    none_zero_indexes = npy.where(completion[1] != 0)
    completionPercentage = npy.zeros(len(completion[0]), dtype=float)
    completionPercentage[none_zero_indexes] = npy.true_divide(completion[0,none_zero_indexes],
                                                              completion[1,none_zero_indexes])

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
            values = npy.ravel(npy.asarray(concordance[sampleID][0][i,indexToKeep]))
            total_values = npy.ravel(npy.asarray(concordance[sampleID][1][i,indexToKeep]))
            currConcordance = npy.zeros(len(indexToKeep), dtype=float)
            none_zero_indexes = npy.where(total_values != 0)
            currConcordance[none_zero_indexes] = npy.true_divide(values[none_zero_indexes],
                                                                 total_values[none_zero_indexes])
            currConcordance = npy.mean(currConcordance)
            toPrint.append("%.8f" % currConcordance)
            print >>outputFile, "\t".join(toPrint)

    # Closing the output file
    outputFile.close()

    return completionPercentage


def computeStatistics(tped, tfam, samples, oldSamples, prefix):
    """Computes the completion and concordance of each samples.

    We don't compute completion and concordance for markers on chromosome 24 is
    sample is female.


    **Completion**

    We consider a genotype as being missing if sample is male and marker is
    heterozygous on chromosome 23 and 24. Also, code for missing genotype is 0.


    **Concordance**

    For each samples, we find the duplicated ones, and compute the mean
    concordance by pairs. For each marker, if both are genotypes are not
    missing, we add one to the total number of compared markers. If both
    genotypes are the same, we add one to the number of concordant calls. We
    write the observed genotype difference in the file ``prefix.diff``.

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
    completion = npy.array([[0 for i in xrange(len(tped[0][4:]))],
                            [0 for i in xrange(len(tped[0][4:]))]])

    # The concordance data type
    concordance = {}
    for sampleID in samples.keys():
        nbDup = len(samples[sampleID])
        concordance[sampleID] = [npy.asmatrix(npy.zeros((nbDup, nbDup),
                                                         dtype=int)),
                                 npy.asmatrix(npy.zeros((nbDup, nbDup),
                                                         dtype=int))]

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
#########################################################
## This is for man on chr 23 and 24, if heterozygous!!! #
#########################################################
##                         if (sex == "1") and (chromosome in ["23", "24"]):
##                             if ("0" not in genotype1) and ("0" not in genotype2):
##                                 if (len(genotype1) != 2) and (len(genotype2) != 2):
##                                     concordance[key][1][i,j] += 1
##                                     concordance[key][1][j,i] += 1
##                                     if genotype1 == genotype2:
##                                         concordance[key][0][i,j] += 1
##                                         concordance[key][0][j,i] += 1
##                                     else:
##                                         # Print to diff file
##                                         toPrint = [snpName, key[0], key[1],
##                                                    str(index+1),
##                                                    str(indexes[j]+1)]
##                                         toPrint.append(" ".join(genotype1))
##                                         if len(genotype1) == 1:
##                                             toPrint[-1] += " %s" % toPrint[-1]
##                                         toPrint.append(" ".join(genotype2))
##                                         if len(genotype2) == 1:
##                                             toPrint[-1] += " %s" % toPrint[-1]
##                                         print >>diffOutput, "\t".join(toPrint)
## 
##                         elif ("0" not in genotype1) and ("0" not in genotype2):
                        if ("0" not in genotype1) and ("0" not in genotype2):
                            # Both have calls
                            concordance[key][1][i,j] += 1
                            concordance[key][1][j,i] += 1
                            if genotype1 == genotype2:
                                concordance[key][0][i,j] += 1
                                concordance[key][0][j,i] += 1
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
            concordance[key][0][i,i] = 1
            concordance[key][1][i,i] = 1

    return completion, concordance


def processTPED(uniqueSamples, duplicatedSamples, fileName, prefix):
    """Process the TPED file."""
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

    tped = npy.array(tped)

    return tped, updatedSamples


def printUniqueTFAM(tfam, samples, prefix):
    """Prints a new TFAM with only unique samples.

    :param tfam: a representation of a TFAM file.
    :type tfam: list of tuple of string

    :param samples: a map of sample position
    :type samples: map of int

    :param prefix: the prefix of the output file name
    :type prefix: string

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

    :param tfam: representation of a TFAM file.
    :type tfam: list of tuple of string

    :returns: two map, containing unique and duplicated samples position.

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

    :param fileName: the name of the TFAM file.
    :type fileName: string

    :returns: a list of tuples, representing the TFAM file.

    """
    # Saving the TFAM file
    tfam = None
    with open(fileName, 'r') as inputFile:
        tfam = [tuple(i.rstrip("\r\n").split("\t")) for i in inputFile.readlines()]

    tfam = npy.array(tfam)

    return tfam


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
    # Checking the input files
    for suffix in [".tped", ".tfam"]:
        fileName = args.tfile + suffix
        if not os.path.isfile(fileName):
            msg = "%(fileName)s: no such file" % locals()
            raise ProgramError(msg)

    # Checking the concordance threshold
    if (args.sample_concordance_threshold < 0) or (args.sample_concordance_threshold > 1):
        msg = "sample-concordance-threshold: must be between 0 and 1 " \
              "(not %f)" % args.sample_concordance_threshold
        raise ProgramError(msg)

    # Checking the completion threshold
    if (args.sample_completion_threshold < 0) or (args.sample_completion_threshold > 1):
        msg = "sample-completion-threshold: must be between 0 and 1 " \
              "(not %f)" % args.sample_completion_threshold
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
prog = "duplicated_samples"
desc = """Extract and work with duplicated samples."""
parser = argparse.ArgumentParser(description=desc, prog=prog)

# The INPUT files
group = parser.add_argument_group("Input File")
group.add_argument("--tfile", type=str, metavar="FILE", required=True,
                   help=("The input file prefix (will find the tped and tfam "
                         "file by appending the prefix to .tped and .tfam, "
                         "respectively."))
# The options
group = parser.add_argument_group("Options")
group.add_argument("--sample-completion-threshold", type=float, metavar="FLOAT",
                   default=0.9, help=("The completion threshold to consider a "
                                      "replicate when choosing the best "
                                      "replicates. [default: %(default).1f]"))
group.add_argument("--sample-concordance-threshold", type=float,
                   metavar="FLOAT", default=0.97,
                   help=("The concordance threshold to consider a replicate "
                         "when choosing the best replicates. [default: "
                         "%(default).2f]"))
# The OUTPUT files
group = parser.add_argument_group("Output File")
group.add_argument("--out", type=str, metavar="FILE",
                   default="dup_samples",
                   help=("The prefix of the output files. [default: "
                         "%(default)s]"))

if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print >>sys.stderr, "Cancelled by user"
        sys.exit(0)
    except ProgramError as e:
        parser.error(e.message)
