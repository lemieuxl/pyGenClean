"""Extracts and merges duplicated markers."""


import argparse
import itertools
import logging
import shutil
from collections import defaultdict
from typing import Dict, List, Optional, Tuple

import numpy as np
from pyplink import PyPlink

from ...error import ProgramError
from ...utils import check_allele_status, compatible_alleles
from ...utils import plink as plink_utils
from ...utils import timer
from ...version import pygenclean_version as __version__


SCRIPT_NAME = "duplicated-markers"
DESCRIPTION = "Extracts and merges duplicated markers."


logger = logging.getLogger(__name__)


@timer(logger)
def main(args: Optional[argparse.Namespace] = None,
         argv: Optional[List[str]] = None) -> None:
    """Checks for duplicated markers.

    Args:
        args (argparse.Namespace): the arguments and options.
        argv (list): the arguments as a list.

    """
    if args is None:
        args = parse_args(argv)
    check_args(args)

    # Reading the map file
    logger.info("Reading BIM file")
    duplicated_markers = parse_bim(args.bfile + ".bim", args.out)

    # If there are no duplicated genomic locations, then the original Plink
    # files are copied over
    if not duplicated_markers:
        logger.info("There are no duplicated genomic locations")
        logger.info("  - Creating final Plink files")
        for extension in (".bed", ".bim", ".fam"):
            shutil.copy(
                args.bfile + extension, args.out + ".final" + extension,
            )
        return

    # Creating binary files with only unique genomic locations
    logger.info("Creating the Plink files containing unique genomic locations")
    plink_utils.subset(
        bfile=args.bfile,
        out=args.out + ".unique_markers",
        markers=args.out + ".unique_markers",
        marker_subset_type="extract",
        use_original_plink=args.plink_107,
    )

    # Creating binary files with only the duplicated genomic locations
    logger.info("Creating the Plink files containing duplicated genomic "
                "locations")
    plink_utils.subset(
        bfile=args.bfile,
        out=args.out + ".duplicated_markers",
        markers=args.out + ".duplicated_markers",
        marker_subset_type="extract",
        use_original_plink=args.plink_107,
    )

    # Computing the frequency of the duplicated SNPs
    logger.info("Computing duplicated markers' frequency")
    frequencies = get_frequencies(
        args.out + ".duplicated_markers",
        args.out + ".duplicated_markers",
        use_plink_107=args.plink_107,
    )

    # Compute statistics
    logger.info("Computing concordance and completion of duplicated SNPs")
    completion, concordance = compute_statistics(
        bfile=args.out + ".duplicated_markers",
        markers=duplicated_markers,
        frequencies=frequencies,
        diff_freq=args.frequency_difference,
        out=args.out,
    )

#         # Print the statistics
#         logger.info("Printing duplicated SNPs summary file and finding errors "
#                     "within duplicates")
#         snpsToComplete = printProblems(completion, concordance, tped,
#                                        duplicatedSNPs, dupSNPsFreq, args.out,
#                                        args.frequency_difference)

    # # Print the concordance file
    # logger.info("Printing concordance file")
    # print_concordance(concordance, args.out, tped, duplicatedSNPs)

#         # Choose the best SNP
#         logger.info("Choosing best SNP for each duplicates")
#         chosenSNPs, comp, conc = chooseBestSnps(tped, duplicatedSNPs,
#                                                 completion, concordance,
#                                                 args.out)

#         # Complete the SNPs
#         logger.info("Completing chosen duplicates (removing discordant "
#                     "genotypes)")
#         newTPED, snpToRemove = createAndCleanTPED(
#             tped,
#             tfam,
#             duplicatedSNPs,
#             args.out,
#             chosenSNPs,
#             comp,
#             conc,
#             snpsToComplete,
#             args.tfile + ".tfam",
#             args.snp_completion_threshold,
#             args.snp_concordance_threshold,
#         )

#         # Creates the final tped
#         logger.info("Writing final TPED and TFAM file")
#         createFinalTPEDandTFAM(newTPED, args.out + ".unique_snps",
#                                args.out, snpToRemove)


# def createFinalTPEDandTFAM(tped, toReadPrefix, prefix, snpToRemove):
#     """Creates the final TPED and TFAM.

#     :param tped: a representation of the ``tped`` of duplicated markers.
#     :param toReadPrefix: the prefix of the unique files.
#     :param prefix: the prefix of the output files.
#     :param snpToRemove: the markers to remove.

#     :type tped: numpy.array
#     :type toReadPrefix: str
#     :type prefix: str
#     :type snpToRemove: set

#     Starts by copying the unique markers' ``tfam`` file to
#     ``prefix.final.tfam``. Then, it copies the unique markers' ``tped`` file,
#     in which the chosen markers will be appended.

#     The final data set will include the unique markers, the chosen markers
#     which were completed, and the problematic duplicated markers (for further
#     analysis). The markers that were used to complete the chosen ones are not
#     present in the final data set.

#     """
#     # First, copying the tfam
#     try:
#         shutil.copy(toReadPrefix + ".tfam", prefix + ".final.tfam")
#     except IOError:
#         msg = "%(toReadPrefix)s.tfam: can't copy file to " \
#               "%(prefix)s.final.tfam" % locals()
#         raise ProgramError(msg)

#     # Next, copy the tped, and append at the end
#     try:
#         shutil.copy(toReadPrefix + ".tped", prefix + ".final.tped")
#     except IOError:
#         msg = "%(toReadPrefix)s.tped: can't copy fil to " \
#               "%(prefix)s.final.tped" % locals()
#         raise ProgramError(msg)
#     tpedFile = None
#     try:
#         tpedFile = open(prefix + ".final.tped", "a")
#     except IOError:
#         msg = "%(prefix)s.final.tped: can't append to file" % locals()
#         raise ProgramError(msg)
#     for i, row in enumerate(tped):
#         if i not in snpToRemove:
#             print >>tpedFile, "\t".join(row)
#     tpedFile.close()


# def createAndCleanTPED(tped, tfam, snps, prefix, chosenSNPs, completion,
#                        concordance, snpsToComplete, tfamFileName, completionT,
#                        concordanceT):
#     """Complete a TPED for duplicated SNPs.

#     :param tped: a representation of the ``tped`` of duplicated markers.
#     :param tfam: a representation of the ``tfam``.
#     :param snps: the position of duplicated markers in the ``tped``.
#     :param prefix: the prefix of the output files.
#     :param chosenSNPs: the markers that were chosen for completion (including
#                        problems).
#     :param completion: the completion of each of the duplicated markers.
#     :param concordance: the pairwise concordance of the duplicated markers.
#     :param snpsToComplete: the markers that will be completed (excluding
#                            problems).
#     :param tfamFileName: the name of the original ``tfam`` file.
#     :param completionT: the completion threshold.
#     :param concordanceT: the concordance threshold.

#     :type tped: numpy.array
#     :type tfam: list
#     :type snps: dict
#     :type prefix: str
#     :type chosenSNPs: dict
#     :type completion: numpy.array
#     :type concordance: dict
#     :type snpsToComplete: set
#     :type tfamFileName: str
#     :type completionT: float
#     :type concordanceT: float

#     :returns: a tuple containing the new ``tped`` after completion
#               (:py:class:`numpy.array` as the first element, and the index of
#               the markers that will need to be rid of (:py:class:`set`) as the
#               last element.

#     It creates three different files:

#     * ``prefix.zeroed_out``: contains information about markers and samples
#                              where the genotyped was zeroed out.
#     * ``prefix.not_good_enough``: contains information about markers that were
#                                   not good enough to help in completing the
#                                   chosen markers (because of concordance or
#                                   completion).
#     * ``prefix.removed_duplicates``: the list of markers that where used for
#                                      completing the chosen one, hence they will
#                                      be removed from the final data
#                                      set.

#     Cycling through every genotypes of every samples of every duplicated
#     markers, checks if the genotypes are all the same. If the chosen one was
#     not called, but the other ones were, then we complete the chosen one with
#     the genotypes for the others (assuming that they are all the same). If
#     there is a difference between the genotypes, it is zeroed out for the
#     chosen marker.

#     """
#     zeroedOutFile = None
#     try:
#         zeroedOutFile = open(prefix + ".zeroed_out", "w")
#     except IOError:
#         msg = "%(prefix).zeroed_out: can't write file" % locals()
#         raise ProgramError(msg)
#     print >>zeroedOutFile, "\t".join(["famID", "indID", "snpID"])

#     notGoodEnoughFile = None
#     try:
#         notGoodEnoughFile = open(prefix + ".not_good_enough", "w")
#     except IOError:
#         msg = "%(prefix)s.not_good_enough: can't write file" % locals()
#         raise ProgramError(msg)
#     print >>notGoodEnoughFile, "\t".join(["name", "reason"])

#     removedFile = None
#     try:
#         removedFile = open(prefix + ".removed_duplicates", "w")
#     except IOError:
#         msg = "%(prefix)s.removed_duplicates: can't write file" % locals()
#         raise ProgramError(msg)

#     notGoodEnoughSnps = set()

#     # Split the tped in 'snpInfo' and 'genotypes'
#     snpInfo = tped[:, :4]
#     genotypes = tped[:, 4:]

#     # The sed of index we want to get rid of at the end
#     getRidOfIndex = set()

#     for snpID, indexes in snps.iteritems():
#         if snpID not in snpsToComplete:
#             # We don't want to complete this SNP, so we continue to next SNP
#             continue

#         # Getting the completion
#         completionToRemove = set(
#             np.where(completion[indexes] < completionT)[0]
#         )
#         for k in completionToRemove:
#             notGoodEnoughSnps.add((snpInfo[indexes][k, 1], "completion"))

#         # Getting the concordance
#         concordanceToRemove = set(
#             np.where(concordance[snpID] < concordanceT)[0]
#         )
#         for k in concordanceToRemove:
#             notGoodEnoughSnps.add((snpInfo[indexes][k, 1], "concordance"))

#         # These will be the indexes to remove
#         indexesToRemove = set()
#         for index in completionToRemove | concordanceToRemove:
#             indexesToRemove.add(indexes[index])

#         # These are the indexes to keep
#         indexesToKeep = []
#         for index in indexes:
#             if index not in indexesToRemove:
#                 indexesToKeep.append(index)

#         # Getting the chosen SNP
#         chosenOne = chosenSNPs[snpID]
#         if chosenOne not in set(indexesToKeep):
#             # The chosen SNP is not a good SNP, so we go to next SNP
#             logger.warning("  - {} chosen but not good enough".format(
#                 snpInfo[chosenOne, 1],
#             ))
#             continue

#         # Now cycling through the genotypes
#         nbSamples = genotypes.shape[1]
#         for sampleIndex in xrange(nbSamples):
#             # We need to remove the no call and keep the unique genotypes
#             curGenotypes = genotypes[indexesToKeep, sampleIndex]
#             cleanedCurGenotypes = curGenotypes[
#                 np.where(curGenotypes != "0 0")
#             ]
#             uniqueCleanedCurGenotypes = np.unique(cleanedCurGenotypes)

#             # Checking the number of unique genotypes
#             toComplete = False
#             if len(uniqueCleanedCurGenotypes) > 1:
#                 # There are more than one unique genotype (except 0 0)
#                 # len = 0 means all were 0 0
#                 # len = 1 means they are all the same
#                 # len > 1 means discordance (might need to flip)
#                 # Just need to check the order of the alleles
#                 possibleAlleles = [
#                     set() for k in xrange(len(uniqueCleanedCurGenotypes))
#                 ]
#                 for k, geno in enumerate(uniqueCleanedCurGenotypes):
#                     possibleAlleles[k] |= set(geno.split(" "))
#                 allEqual = True
#                 for k in xrange(len(possibleAlleles)):
#                     for l in xrange(k+1, len(possibleAlleles)):
#                         if possibleAlleles[k] != possibleAlleles[l]:
#                             allEqual = False

#                 if not allEqual:
#                     # The genotypes are not all equal, we set the chosen
#                     # genotype to null (0 0)
#                     tped[chosenOne, sampleIndex+4] = "0 0"
#                     print >>zeroedOutFile, "\t".join([tfam[sampleIndex, 0],
#                                                       tfam[sampleIndex, 1],
#                                                       snpInfo[chosenOne, 1]])
#                 elif genotypes[chosenOne, sampleIndex] == "0 0":
#                     toComplete = True
#             elif ((len(uniqueCleanedCurGenotypes) == 1) and
#                     (genotypes[chosenOne, sampleIndex] == "0 0")):
#                 toComplete = True

#             if toComplete:
#                 # We complete the current individual
#                 tped[chosenOne, sampleIndex+4] = uniqueCleanedCurGenotypes[0]

#         # We keep only the chose one
#         for index in indexes:
#             if index != chosenOne:
#                 getRidOfIndex.add(index)
#                 print >>removedFile, snpInfo[index, 1]

#     # Writing the not good enough file
#     for item in notGoodEnoughSnps:
#         print >>notGoodEnoughFile, "\t".join(item)

#     # Closing the output files
#     zeroedOutFile.close()
#     notGoodEnoughFile.close()

#     # Printing the chosen file
#     try:
#         shutil.copy(tfamFileName, prefix + ".chosen_snps.tfam")
#     except IOError:
#         msg = "%(tfamFileName)s: can't copy file to " \
#               "%(prefix)s.chosen_snps.tfam" % locals()
#         raise ProgramError(msg)
#     chosenFile = None
#     try:
#         chosenFile = open(prefix + ".chosen_snps.tped", "w")
#     except IOError:
#         msg = "%(prefix)s.chosen_snps.tped: can't write file" % locals()
#         raise ProgramError(msg)
#     for chosenOne in chosenSNPs.itervalues():
#         snpID = (tped[chosenOne, 0], tped[chosenOne, 3])
#         if snpID in snpsToComplete:
#             print >>chosenFile, "\t".join(tped[chosenOne])
#     chosenFile.close()

#     return tped, getRidOfIndex


# def chooseBestSnps(tped, snps, trueCompletion, trueConcordance, prefix):
#     """Choose the best duplicates according to the completion and concordance.

#     :param tped: a representation of the ``tped`` of duplicated markers.
#     :param snps: the position of the duplicated markers in the ``tped``.
#     :param trueCompletion: the completion of each markers.
#     :param trueConcordance: the pairwise concordance of each markers.
#     :param prefix: the prefix of the output files.

#     :type tped: numpy.array
#     :type snps: dict
#     :type trueCompletion: numpy.array
#     :type trueConcordance: dict
#     :type prefix: str

#     :returns: a tuple containing the chosen indexes (:py:class:`dict`) as the
#               first element, the completion (:py:class:`numpy.array`) as the
#               second element, and the concordance (:py:class:`dict`) as last
#               element.

#     It creates two output files: ``prefix.chosen_snps.info`` and
#     ``prefix.not_chosen_snps.info``. The first one contains the markers that
#     were chosen for completion, and the second one, the markers that weren't.

#     It starts by computing the completion of each markers (dividing the number
#     of calls divided by the total number of genotypes). Then, for each of the
#     duplicated markers, we choose the best one according to completion and
#     concordance (see explanation in
#     :py:func:`DupSamples.duplicated_samples.chooseBestDuplicates` for more
#     details).

#     """
#     # The output files
#     chosenFile = None
#     try:
#         chosenFile = open(prefix + ".chosen_snps.info", "w")
#     except IOError:
#         msg = "%(prefix)s.chosen_snps.info: can't write file" % locals()
#         raise ProgramError(msg)

#     excludedFile = None
#     try:
#         excludedFile = open(prefix + ".not_chosen_snps.info", "w")
#     except IOError:
#         msg = "%(prefix)s.not_chosen_snps.info: can't write file" % locals()
#         raise ProgramError(msg)

#     # Computing the completion
#     completion = np.true_divide(trueCompletion[0], trueCompletion[1])

#     # For each duplicated SNPs
#     chosenIndexes = {}
#     snpConcordance = {}
#     for snp, indexes in snps.iteritems():
#         # Getting the completion for those duplicated SNPs
#         currCompletion = completion[indexes]

#         # Sorting those completion
#         sortedCompletionInsexes = np.argsort(currCompletion)

#         # Getting the concordance
#         concordance = np.true_divide(trueConcordance[snp][0],
#                                      trueConcordance[snp][1])

#         currConcordance = [[] for i in xrange(len(indexes))]
#         for i in xrange(len(indexes)):
#             indexToKeep = list(set(range(len(indexes))) - set([i]))
#             currConcordance[i] = np.mean(concordance[i, indexToKeep])
#         currConcordance = np.array(currConcordance)
#         if snp not in snpConcordance:
#             snpConcordance[snp] = currConcordance

#         # Sorting the concordance
#         sortedConcordanceIndexes = np.argsort(currConcordance)

#         # Trying to find the best duplicate to keep
#         nbToCheck = 1
#         chosenIndex = None
#         while nbToCheck <= len(indexes):
#             # Getting the `nbToCheck` best value (higher to lower)
#             completionValue = currCompletion[
#                 sortedCompletionInsexes[nbToCheck*-1]
#             ]
#             concordanceValue = currConcordance[
#                 sortedConcordanceIndexes[nbToCheck*-1]
#             ]

#             # Getting the indexes to consider
#             completionToConsider = set(
#                 np.where(currCompletion >= completionValue)[0]
#             )
#             concordanceToConsider = set(
#                 np.where(currConcordance >= concordanceValue)[0]
#             )

#             # Getting the intersection of the indexes
#             toConsider = concordanceToConsider & completionToConsider
#             if len(toConsider) >= 1:
#                 chosenIndex = random.choice(list(toConsider))
#                 break
#             nbToCheck += 1

#         if chosenIndex is None:
#             msg = "Could not choose the best snp ID"
#             raise ProgramError(msg)

#         # Printing the chosen SNPs
#         print >>chosenFile, tped[indexes[chosenIndex], 1]

#         # Printing the excluded SNPs
#         for i, index in enumerate(indexes):
#             if i != chosenIndex:
#                 print >>excludedFile, tped[index, 1]

#         chosenIndexes[snp] = indexes[chosenIndex]

#     # Closing the output files
#     chosenFile.close()
#     excludedFile.close()

#     return chosenIndexes, completion, snpConcordance


def get_frequencies(
    bfile: str,
    out: str,
    use_plink_107: bool,
) -> Dict[str, Tuple[float, Tuple[str, str]]]:
    """Computes the frequency of the markers using Plink.

    Args:
        bfile (str): the prefix of the input files.
        out (str): the prefix of the output files.

    Returns:
        dict: the frequencies of each marker along with the two alleles.

    Start by computing the frequency of all markers using Plink. Then, it reads
    the output file, and saves the frequency and allele information.

    """
    plink_utils.compute_freq(
        bfile=bfile, out=out, use_original_plink=use_plink_107,
    )

    # Reading the frequency file
    frequencies: Dict[str, Tuple[float, Tuple[str, str]]] = {}
    with open(out + ".frq", "r") as f:
        header = None
        for line in f:
            row = plink_utils.split_line(line)

            if header is None:
                header = {name: i for i, name in enumerate(row)}
                continue

            # Getting the MAF (if NA, we set as 0)
            maf = row[header["MAF"]]
            if maf == "NA":
                maf = 0.0
            else:
                maf = float(maf)

            frequencies[row[header["SNP"]]] = (
                maf, (row[header["A1"]], row[header["A2"]]),
            )

    return frequencies


# def print_concordance(concordance, prefix, tped, snps):
#     """Print the concordance.

#     :param concordance: the concordance.
#     :param prefix: the prefix if the output files.
#     :param tped: a representation of the ``tped`` of duplicated markers.
#     :param snps: the position of the duplicated markers in the ``tped``.

#     :type concordance: dict
#     :type prefix: str
#     :type tped: numpy.array
#     :type snps: dict

#     Prints the concordance in a file, in the format of a matrix. For each
#     duplicated markers, the first line (starting with the `#` signs) contains
#     the name of all the markers in the duplicated markers set. Then a :math:`N
#     \\times N` matrix is printed to file (where :math:`N` is the number of
#     markers in the duplicated marker list), containing the pairwise
#     concordance.

#     """
#     outFile = None
#     try:
#         outFile = open(prefix + ".concordance", "w")
#     except IOError:
#         msg = "%s: can't write file" % prefix + ".concordance"
#         raise ProgramError(msg)

#     for snpID in concordance.iterkeys():
#         print >>outFile, "#" + "\t".join(
#             list(snpID) + list(tped[snps[snpID], 1])
#         )

#         # Doing the division
#         true_concordance = np.true_divide(concordance[snpID][0],
#                                           concordance[snpID][1])

#         output = StringIO.StringIO()
#         np.savetxt(output, true_concordance, delimiter="\t", fmt="%.8f")
#         print >>outFile, output.getvalue().rstrip("\r\n")

#     outFile.close()


# def printProblems(completion, concordance, tped, snps, frequencies, prefix,
#                   diffFreq):
#     """Print the statistics.

#     :param completion: the completion of each duplicated markers.
#     :param concordance: the pairwise concordance between duplicated markers.
#     :param tped: a representation of the ``tped`` of duplicated markers.
#     :param snps: the positions of the duplicated markers in the ``tped``
#     :param frequencies: the frequency of each of the duplicated markers.
#     :param prefix: the prefix of the output files.
#     :param diffFreq: the frequency difference threshold.

#     :type completion: numpy.array
#     :type concordance: dict
#     :type tped: numpy.array
#     :type snps: dict
#     :type frequencies: dict
#     :type prefix: str
#     :type diffFreq: float

#     :returns: a :py:class:`set` containing duplicated markers to complete.

#     Creates a summary file (``prefix.summary``) containing information about
#     duplicated markers: chromosome, position, name, alleles, status, completion
#     percentage, completion number and mean concordance.

#     The frequency and the minor allele are used to be certain that two
#     duplicated markers are exactly the same marker (and not a tri-allelic one,
#     for example).

#     For each duplicated markers:

#     1. Constructs the set of available alleles for the first marker.
#     2. Constructs the set of available alleles for the second marker.
#     3. If the two sets are different, but the number of alleles is the same, we
#        try to flip one of the marker. If the two sets are the same, but the
#        number of alleles is 1, we set the status to ``homo_flip``. If the
#        markers are heterozygous, we set the status to ``flip``.
#     4. If there is a difference in the number of alleles (one is homozygous,
#        the other, heterozygous), and that there is on allele in common, we set
#        the status to ``homo_hetero``. If there are no allele in common, we try
#        to flip one. If the new sets have one allele in common, we set the
#        status to ``homo_hetero_flip``.
#     5. If the sets of available alleles are the same (without flip), we check
#        the frequency and the minor alleles. If the minor allele is different,
#        we set the status to ``diff_minor_allele``. If the difference in
#        frequencies is higher than a threshold, we set the status to
#        ``diff_frequency``.
#     6. If all of the above fail, we set the status to ``problem``.

#     Problems are written in the ``prefix.problems`` file, and contains the
#     following columns: chromosome, position, name and status. This file
#     contains all the markers with a status, as explained above.

#     """

#     completionPercentage = np.true_divide(completion[0], completion[1])

#     outSummary = None
#     try:
#         outSummary = open(prefix + ".summary", "w")
#     except IOError:
#         msg = "%s: can't write file" % prefix + ".summary"
#         raise ProgramError

#     # Prints the header of the summary file
#     print >>outSummary, "\t".join(["chr", "pos", "name", "alleles", "status",
#                                    "% completion", "completion",
#                                    "mean concordance"])

#     # The data structure containing the problems
#     problems = {}
#     for snpID, indexes in snps.iteritems():
#         for i, index in enumerate(indexes):
#             # The SNP information (chromosome and position)
#             toPrint = list(snpID)

#             # The name of the SNP
#             snpName = tped[index, 1]
#             toPrint.append(snpName)

#             # The frequency of the SNP
#             snpFreq, mafAlleles = frequencies[snpName]

#             # A list of the other SNP name with problems
#             otherSnpNameWithProblem = set()

#             # The alleles
#             alleles = set()
#             otherAlleles = set()
#             status = []
#             for genotype in np.unique(tped[index, 4:]):
#                 alleles |= set(genotype.split(" "))
#             if "0" in alleles:
#                 alleles.remove("0")
#             for j in xrange(i+1, len(indexes)):
#                 otherIndex = indexes[j]
#                 otherSnpName = tped[otherIndex, 1]

#                 # The frequency of the other SNP
#                 otherSnpFreq, otherMafAlleles = frequencies[otherSnpName]

#                 # Checking the alleles
#                 for genotype in np.unique(tped[otherIndex, 4:]):
#                     otherAlleles |= set(genotype.split(" "))
#                 if "0" in otherAlleles:
#                     otherAlleles.remove("0")
#                 if alleles != otherAlleles:
#                     if len(alleles) == len(otherAlleles):
#                         # Same number of alleles
#                         # Try the flipped ones
#                         otherAlleles = flipGenotype(otherAlleles)
#                         if alleles == otherAlleles:
#                             if len(alleles) == 1:
#                                 status.append("homo_flip")
#                                 otherSnpNameWithProblem.add(otherSnpName)
#                             else:
#                                 status.append("flip")
#                                 otherSnpNameWithProblem.add(otherSnpName)
#                         else:
#                             status.append("problem")
#                             otherSnpNameWithProblem.add(otherSnpName)
#                     else:
#                         # Different number of alleles
#                         if len(alleles & otherAlleles) == 1:
#                             status.append("homo_hetero")
#                             otherSnpNameWithProblem.add(otherSnpName)
#                         else:
#                             # Try the flipped one
#                             otherAlleles = flipGenotype(otherAlleles)
#                             if len(alleles & otherAlleles) == 1:
#                                 status.append("homo_hetero_flip")
#                                 otherSnpNameWithProblem.add(otherSnpName)
#                             else:
#                                 status.append("problem")
#                                 otherSnpNameWithProblem.add(otherSnpName)
#                 else:
#                     # The alleles are the same, so we check the frequency
#                     if mafAlleles[0] != otherMafAlleles[0]:
#                         # They don't have the same minor allele
#                         status.append("diff_minor_allele")
#                         otherSnpNameWithProblem.add(otherSnpName)
#                     elif math.fabs(snpFreq - otherSnpFreq) > diffFreq:
#                         # They don't have same frequency
#                         status.append("diff_frequency")
#                         otherSnpNameWithProblem.add(otherSnpName)

#             alleles = list(alleles)
#             alleles.sort()
#             if len(alleles) == 1:
#                 alleles.append(alleles[0])
#             toPrint.append(" ".join(alleles))
#             toPrint.append(";".join(status))

#             # The completion
#             toPrint.append("%.8f" % completionPercentage[index])
#             toPrint.append("%d/%d" % (completion[0][index],
#                                       completion[1][index]))

#             # The concordance
#             indexToKeep = list(set(range(len(indexes))) - set([i]))
#             currConcordance = np.true_divide(
#                 concordance[snpID][0][i, indexToKeep],
#                 concordance[snpID][1][i, indexToKeep],
#             )
#             currConcordance = np.mean(currConcordance)
#             toPrint.append("%.8f" % currConcordance)
#             print >>outSummary, "\t".join(toPrint)

#             # Now updating the problems data structure
#             if len(status) != len(otherSnpNameWithProblem):
#                 msg = "There is a problem with the problematic SNPs"
#                 raise ProgramError(msg)

#             if len(status) > 0:
#                 if snpID not in problems:
#                     tmp = {"snpNames": {snpName}, "problems": set()}
#                     problems[snpID] = tmp

#                 # We have problems
#                 problems[snpID]["snpNames"] |= otherSnpNameWithProblem
#                 problems[snpID]["problems"] |= set(status)

#     outSummary.close()

#     outProblems = None
#     try:
#         outProblems = open(prefix + ".problems", "w")
#     except IOError:
#         msg = "%s: can't write file" % prefix + ".problems"
#         raise ProgramError

#     # Printing the header of the problem file...
#     print >>outProblems, "\t".join(["chr", "pos", "name", "status"])
#     for snpID in problems.iterkeys():
#         toPrint = list(snpID)
#         toPrint.append(";".join(list(problems[snpID]["snpNames"])))
#         toPrint.append(";".join(list(problems[snpID]["problems"])))
#         print >>outProblems, "\t".join(toPrint)

#     outProblems.close()

#     # Returning the SNPs to complete
#     return set(snps.keys()) - set(problems.keys())


def compute_statistics(
    bfile: str,
    markers: Dict[Tuple[int, int], List[str]],
    frequencies: Dict[str, Tuple[float, Tuple[str, str]]],
    diff_freq: float,
    out: str,
) -> Tuple[
    Dict[Tuple[int, int], np.ndarray],
    Dict[Tuple[int, int], np.ndarray]
]:
    """Computes the completion and concordance of each SNPs.

    Args:
        bfile (str): the input file prefix.
        markers (dict): the duplicated markers.
        out (str): the output prefix.

    Returns:
        tuple: a tuple containing the completion of duplicated markers
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
    # The completion
    completion = {
        location: np.zeros(len(duplicates), dtype=float)
        for location, duplicates in markers.items()
    }

    # The concordance
    concordance = {
        location: np.zeros((len(duplicates), len(duplicates)), dtype=float)
        for location, duplicates in markers.items()
    }

    # The problems
    problems = {}

    with PyPlink(bfile) as bed:
        bim = bed.get_bim()

        for location, marker_names in markers.items():
            # Getting the genotypes of each of the duplicated markers
            genotypes = np.vstack(tuple(
                bed.get_geno_marker(marker_name)
                for marker_name in marker_names
            ))

            # Computing the missing masks once
            not_missing_mask = genotypes != -1

            # Computing the completion for this genomic locations
            completion[location] = (
                np.count_nonzero(not_missing_mask, axis=1) / genotypes.shape[1]
            )

            # Comparing markers with same genomic locations two by two
            full_status = []
            for i, j in itertools.combinations(range(genotypes.shape[0]), 2):
                # Getting the alleles for each marker
                alleles_i = set(bim.loc[marker_names[i], ["a1", "a2"]])
                alleles_j = set(bim.loc[marker_names[j], ["a1", "a2"]])

                # Getting only non-missing genotypes for both markers
                not_missing = not_missing_mask[i] & not_missing_mask[j]

                # Checking the allele status
                allele_status = check_allele_status(alleles_i, alleles_j)

                if allele_status:
                    full_status.append((
                        allele_status,
                        (marker_names[i], marker_names[j]),
                    ))

                # Skipping incompatible alleles
                if not compatible_alleles(alleles_i, alleles_j):
                    assert allele_status == "problem"
                    continue

                # The frequencies
                print(frequencies[marker_names[i]])
                print(frequencies[marker_names[j]])
                print()

                # Counting the number of same genotypes (and with flipped)
                # We only keep the maximum number found
                nb_same = max(
                    np.count_nonzero(
                        (genotypes[i] == genotypes[j]) & not_missing
                    ),
                    np.count_nonzero(
                        (genotypes[i] == (2 - genotypes[j])) & not_missing
                    ),
                )

                # Saving the concordance
                nb_not_missing = np.count_nonzero(not_missing)
                if nb_not_missing:
                    concordance[location][i, j] = nb_same / nb_not_missing
                    concordance[location][j, i] = concordance[location][i, j]

            # Saving problems
            if full_status:
                problems[location] = full_status

    # Saving the problem file
    logger.info("Saving problems to '%s'", out + ".problems")
    with open(out + ".problems", "w") as f:
        print("chr", "pos", "name", "status", sep="\t", file=f)
        for location in sorted(problems.keys()):
            # The status and name sets
            status_set = {status[0] for status in problems[location]}
            name_set = set(itertools.chain.from_iterable(
                status[1] for status in problems[location]
            ))
            print(*location, ";".join(sorted(name_set)),
                  ";".join(sorted(status_set)), sep="\t", file=f)

    # Diagonal should be one for the concordance
    for location in concordance:
        diag_mask = np.eye(concordance[location].shape[0], dtype=bool)
        concordance[location][diag_mask] = 1

    return completion, concordance


# def getIndexOfHeteroMen(genotypes, menIndex):
#     """Get the indexes of heterozygous men.

#     :param genotypes: the genotypes of everybody.
#     :param menIndex: the indexes of the men (for the genotypes).

#     :type genotypes: numpy.array
#     :type menIndex: numpy.array

#     :returns: a :py:class:`numpy.array` containing the indexes of the genotypes
#               to remove.

#     Finds the mean that have a heterozygous genotype for this current marker.
#     Usually used on sexual chromosomes.

#     """
#     toRemove = set()
#     for i in menIndex[0]:
#         for genotype in [set(j.split(" ")) for j in genotypes[:, i]]:
#             if len(genotype) != 1:
#                 # We have an heterozygous
#                 toRemove.add(i)

#     toRemove = list(toRemove)
#     toRemove.sort()

#     return (np.array(toRemove, dtype=int),)


# def flipGenotype(genotype):
#     """Flips a genotype.

#     :param genotype: the genotype to flip.

#     :type genotype: set

#     :returns: the new flipped genotype (as a :py:class:`set`)

#     .. testsetup::

#         from pyGenClean.DupSNPs.duplicated_snps import flipGenotype

#     .. doctest::

#         >>> flipGenotype({"A", "T"})
#         set(['A', 'T'])
#         >>> flipGenotype({"C", "T"})
#         set(['A', 'G'])
#         >>> flipGenotype({"T", "G"})
#         set(['A', 'C'])
#         >>> flipGenotype({"0", "0"})
#         Traceback (most recent call last):
#             ...
#         ProgramError: 0: unkown allele
#         >>> flipGenotype({"A", "N"})
#         Traceback (most recent call last):
#             ...
#         ProgramError: N: unkown allele

#     """
#     newGenotype = set()
#     for allele in genotype:
#         if allele == "A":
#             newGenotype.add("T")
#         elif allele == "C":
#             newGenotype.add("G")
#         elif allele == "T":
#             newGenotype.add("A")
#         elif allele == "G":
#             newGenotype.add("C")
#         else:
#             msg = "%(allele)s: unknown allele" % locals()
#             raise ProgramError(msg)

#     return newGenotype


# # def processTPED(uniqueSNPs, duplicatedSNPs, mapF, fileName, tfam, prefix):
# def processTPED(uniqueSNPs, mapF, fileName, tfam, prefix):
#     """Process the TPED file.

#     :param uniqueSNPs: the unique markers.
#     :param mapF: a representation of the ``map`` file.
#     :param fileName: the name of the ``tped`` file.
#     :param tfam: the name of the ``tfam`` file.
#     :param prefix: the prefix of all the files.

#     :type uniqueSNPs: dict
#     :type mapF: list
#     :type fileName: str
#     :type tfam: str
#     :type prefix: str

#     :returns: a tuple with the representation of the ``tped`` file
#               (:py:class:`numpy.array`) as first element, and the updated
#               position of the duplicated markers in the ``tped``
#               representation.

#     Copies the ``tfam`` file into ``prefix.unique_snps.tfam``. While reading
#     the ``tped`` file, creates a new one (``prefix.unique_snps.tped``)
#     containing only unique markers.

#     """
#     # Copying the tfam file
#     try:
#         shutil.copy(tfam, prefix + ".unique_snps.tfam")
#     except IOError:
#         msg = "%s: can't write file" % prefix + ".unique_snps.tfam"
#         raise ProgramError(msg)

#     tped = []
#     updatedSNPs = defaultdict(list)
#     outputFile = None
#     try:
#         outputFile = open(prefix + ".unique_snps.tped", "w")
#     except IOError:
#         msg = "%s: can't write to file" % prefix + ".unique_snps.tped"
#         raise ProgramError(msg)
#     nbSNP = 0
#     with open(fileName, 'r') as inputFile:
#         for line in inputFile:
#             nbSNP += 1
#             row = line.rstrip("\r\n").split("\t")
#             snpInfo = row[:4]
#             genotype = [i.upper() for i in row[4:]]

#             chromosome = snpInfo[0]
#             position = snpInfo[3]

#             if (chromosome, position) in uniqueSNPs:
#                 # Printing the new TPED file (unique SNPs only)
#                 print >>outputFile, "\t".join(snpInfo + genotype)
#             else:
#                 # Saving the TPED file (duplicated samples only)
#                 currPos = len(tped)
#                 tped.append(tuple(snpInfo + genotype))
#                 updatedSNPs[(chromosome, position)].append(currPos)
#     outputFile.close()

#     if len(mapF) != nbSNP:
#         msg = "%(fileName)s: no the same number of SNPs than MAP " \
#               "file" % locals()
#         raise ProgramError(msg)

#     tped = np.array(tped)

#     return tped, updatedSNPs


# def findUniques(mapF):
#     """Finds the unique markers in a MAP.

#     :param mapF: representation of a ``map`` file.

#     :type mapF: list

#     :returns: a :py:class:`dict` containing unique markers (according to their
#               genomic localisation).

#     """
#     uSNPs = {}
#     dSNPs = defaultdict(list)
#     for i, row in enumerate(mapF):
#         chromosome = row[0]
#         position = row[3]
#         snpID = (chromosome, position)
#         if snpID not in uSNPs:
#             # This is the first time we see this sample
#             uSNPs[snpID] = i
#         else:
#             # We have seen this sample at least once...
#             if snpID not in dSNPs:
#                 # This is the second time we see this sample...
#                 dSNPs[snpID].extend([uSNPs[snpID], i])
#             else:
#                 # We have seen this sample multiple times
#                 dSNPs[snpID].append(i)

#     # Removing the duplicates from the unique samples
#     for snpID in dSNPs.iterkeys():
#         if snpID in uSNPs:
#             del uSNPs[snpID]

#     return uSNPs


def parse_bim(filename: str, out: str) -> Dict[Tuple[int, int], List[str]]:
    """Reads the BIM file.

    Args:
        filename (str): the name of the ``BIM`` file.
        out (str): the output prefix.

    Returns:
        dict: the duplicated markers, according to genomic location.

    While reading the ``BIM`` file, it saves a file
    (``out.duplicated_marker_names``) containing the name of the unique
    duplicated markers (i.e. markers with same names, different positions).

    It also saves the list of unique and duplicated markers for further
    extraction.

    """
    # Sets to keep track of markers' name (and find duplicated ones, if any)
    marker_names = set()
    markers_with_same_name = set()

    # The unique and duplicated markers (according to genomic location)
    unique_markers: Dict[Tuple[int, int], str] = {}
    duplicated_markers = defaultdict(list)

    # Parsing the BIM file
    for marker in plink_utils.parse_bim(filename):
        # Checking for duplicates in names
        if marker.name in marker_names:
            markers_with_same_name.add(marker.name)
        else:
            marker_names.add(marker.name)

        # Checking for duplicates in genomic locations
        genomic_location = (marker.chrom, marker.pos)

        # Already a duplicated genomic location
        if genomic_location in duplicated_markers:
            duplicated_markers[genomic_location].append(marker.name)

        # A "new" duplicated genomic location
        elif genomic_location in unique_markers:
            duplicated_markers[genomic_location].append(marker.name)
            duplicated_markers[genomic_location].append(
                unique_markers.pop(genomic_location),
            )

        # First time we saw this genomic location
        else:
            unique_markers[genomic_location] = marker.name

    # Are there any markers with the same names?
    if markers_with_same_name:
        logger.info("  - found %d markers with same name",
                    len(markers_with_same_name))
        with open(out + ".duplicated_marker_names", "w") as f:
            print(*markers_with_same_name, sep="\n", file=f)

    # Saving the unique markers to file
    logger.info("  - found %d unique genomic locations",
                len(unique_markers))
    with open(out + ".unique_markers", "w") as f:
        print(*unique_markers.values(), sep="\n", file=f)

    # Saving the duplicated markers to file
    logger.info("  - found %d duplicated genomic locations",
                len(duplicated_markers))
    with open(out + ".duplicated_markers", "w") as f:
        for markers in duplicated_markers.values():
            print(*markers, sep="\n", file=f)

    return duplicated_markers


def check_args(args: argparse.Namespace) -> None:
    """Checks the arguments and options.

    Args:
        args (argparse.Namespace): the arguments and options.

    If there is a problem with an option, an exception is raised using the
    :py:class:`ProgramError` class, a message is printed to the
    :class:`sys.stderr` and the program exists with code 1.

    """
    if not plink_utils.check_files(args.bfile):
        raise ProgramError(f"{args.bfile}: no such fules.")

    # Checking the concordance threshold
    if not 0 <= args.snp_concordance_threshold <= 1:
        raise ProgramError(
            f"snp-concordance-threshold: must be between 0 and 1 "
            f"(not {args.snp_concordance_threshold})",
        )

    # Checking the completion threshold
    if not 0 <= args.snp_completion_threshold <= 1:
        raise ProgramError(
            f"snp-completion-threshold: must be between 0 and 1 "
            f"(not {args.snp_completion_threshold})",
        )

    # Checking the difference in frequency
    if not 0 <= args.frequency_difference <= 1:
        raise ProgramError(
            f"{args.frequency_difference}: maximal frequency difference: "
            f"value must be between 0 and 1 (inclusively)",
        )


def parse_args(argv: Optional[List[str]] = None) -> argparse.Namespace:
    """Parses the command line options and arguments."""
    parser = argparse.ArgumentParser(description=DESCRIPTION)

    parser.add_argument(
        "-v", "--version", action="version",
        version=f"pyGenClean {SCRIPT_NAME} {__version__}",
    )

    # Adding the arguments and options
    add_args(parser)

    return parser.parse_args(argv)


def add_args(parser: argparse.ArgumentParser) -> None:
    """Add arguments and options to the parser."""
    # The INPUT files
    group = parser.add_argument_group("Input File")
    group.add_argument(
        "--bfile", type=str, metavar="FILE", required=True,
        help="The input file prefix (will find the plink binary files by "
             "appending the prefix to the .bim, .bed and .fam files, "
             "respectively).",
    )

    # The options
    group = parser.add_argument_group("Options")
    group.add_argument(
        "--plink-1.07", dest="plink_107", action="store_true",
        help="Use original Plink (version 1.07)",
    )
    group.add_argument(
        "--snp-completion-threshold", type=float, metavar="FLOAT", default=0.9,
        help="The completion threshold to consider a replicate when choosing "
             "the best replicates and for composite creation. "
             "[default: %(default).1f]",
    )
    group.add_argument(
        "--snp-concordance-threshold", type=float, metavar="FLOAT",
        default=0.98,
        help="The concordance threshold to consider a replicate when choosing "
             "the best replicates and for composite creation. "
             "[default: %(default).2f]",
    )
    group.add_argument(
        "--frequency-difference", type=float, metavar="FLOAT", default=0.05,
        help="The maximum difference in frequency between duplicated markers "
             "[default: %(default).2f]",
    )

    # The OUTPUT files
    group = parser.add_argument_group("Output File")
    group.add_argument(
        "--out", type=str, metavar="FILE", default="dup_markers",
        help="The prefix of the output files. [default: %(default)s]",
    )
