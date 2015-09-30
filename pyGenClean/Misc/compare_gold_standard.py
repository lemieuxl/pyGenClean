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


import re
import os
import sys
import gzip
import shutil
import logging
import argparse
import subprocess

from .. import __version__
from ..PlinkUtils import createRowFromPlinkSpacedOutput
from ..DupSNPs.duplicated_snps import flipGenotype


logger = logging.getLogger("compare_gold_standard")


def main(argString=None):
    # Getting and checking the options
    args = parseArgs(argString)
    checkArgs(args)

    logger.info("Options used:")
    for key, value in vars(args).iteritems():
        logger.info("  --{} {}".format(key.replace("_", "-"), value))

    # Reading the same sample file
    logger.info("Reading the same samples file")
    same_samples = read_same_samples_file(args.same_samples, args.out)
    required_samples = set([i[0] for i in same_samples] +
                           [i[1] for i in same_samples])
    logger.info("  - Found a total of {} unique "
                "samples".format(len(required_samples)))

    # Check the fam files
    logger.info("Checking in fam files for required samples")
    if not check_fam_for_samples(required_samples, args.bfile + ".fam",
                                 args.gold_bfile + ".fam"):
        logger.info("  - Didn't find all required samples... STOP")
        sys.exit(0)

    # Find overlap markers with the gold standard file
    logger.info("Finding overlapping SNPs with gold standard")
    findOverlappingSNPsWithGoldStandard(args.bfile, args.gold_bfile, args.out,
                                        args.use_marker_names)

    # Extract the required SNPs using Plink
    logger.info("Extracting overlapping SNPs from the gold standard and the "
                "source panel")
    extractSNPs([args.gold_bfile, args.bfile],
                [args.out + ".gold_snp_to_extract",
                 args.out + ".source_snp_to_extract"],
                [args.out + ".gold_standard", args.out + ".source_panel"],
                args.sge)

    gold_input_prefix = args.out + ".gold_standard"
    # Renaiming the reference file, so that the SNP names are the same
    if not args.use_marker_names:
        logger.info("Renaming gold standard's SNPs to match source panel")
        renameSNPs(args.out + ".gold_standard", args.out + ".update_names",
                   args.out + ".gold_standard.rename")
        gold_input_prefix = args.out + ".gold_standard.rename"

    # Finding the SNPs to flip and flip them in the gold standard
    if not args.do_not_flip:
        # Computing the frequency
        logger.info("Computing gold standard frequencies")
        computeFrequency(gold_input_prefix,
                         gold_input_prefix + ".frequency")

        source_alleles = None
        if args.source_manifest is not None:
            logger.info("Reading the source manifest")
            source_alleles = read_source_manifest(args.source_manifest)
        else:
            logger.info("Reading the allele file")
            source_alleles = read_source_alleles(args.source_alleles)

        logger.info("Finding SNPs to flip or to exclude from gold standard")
        findFlippedSNPs(gold_input_prefix + ".frequency.frq", source_alleles,
                        args.out)

        # Excluding SNPs
        logger.info("Excluding SNPs and samples from the gold standard")
        exclude_SNPs_samples(gold_input_prefix,
                             gold_input_prefix + ".cleaned",
                             args.out + ".snp_to_remove",
                             args.out + ".gold_samples2keep")
        logger.info("Excluding SNPs and samples from source panel")
        exclude_SNPs_samples(args.out + ".source_panel",
                             args.out + ".source_panel.cleaned",
                             args.out + ".snp_to_remove",
                             args.out + ".source_panel_samples2keep")

        # Flipping the SNP that need to be flip in the gold standard
        logger.info("Flipping SNPs in gold standard")
        flipSNPs(gold_input_prefix + ".cleaned",
                 gold_input_prefix + ".cleaned.flipped",
                 args.out + ".snp_to_flip_in_reference")

        # Splitting the files, and running the duplicated samples script
        logger.info("Computing statistics")
        compute_statistics(args.out + ".duplicated_samples",
                           gold_input_prefix + ".cleaned.flipped",
                           args.out + ".source_panel.cleaned", same_samples,
                           args.sge, args.out)
    else:
        # We do not need to flip...
        # Splitting the files, and running the duplicated samples script
        logger.info("Keeping samples from the gold standard")
        exclude_SNPs_samples(gold_input_prefix,
                             gold_input_prefix + ".cleaned",
                             keepSample=args.out + ".gold_samples2keep")
        logger.info("Keeping samples from source panel")
        exclude_SNPs_samples(
            args.out + ".source_panel",
            args.out + ".source_panel.cleaned",
            keepSample=args.out + ".source_panel_samples2keep",
        )
        logger.info("Computing statistics")
        compute_statistics(args.out + ".duplicated_samples",
                           gold_input_prefix + ".cleaned",
                           args.out + ".source_panel.cleaned", same_samples,
                           args.sge, args.out)


def read_source_manifest(file_name):
    """Reads Illumina manifest."""
    alleles = {}
    open_function = open
    if file_name.endswith(".gz"):
        open_function = gzip.open
    with open_function(file_name, 'rb') as input_file:
        header_index = dict([
            (col_name, i) for i, col_name in
            enumerate(input_file.readline().rstrip("\r\n").split("\t"))
        ])
        for col_name in {"Name", "SNP", "IlmnStrand"}:
            if col_name not in header_index:
                msg = "{}: no column named {}".format(file_name, col_name)
                raise ProgramError(msg)

        for line in input_file:
            row = line.rstrip("\r\n").split("\t")
            marker_name = row[header_index["Name"]]
            the_alleles = row[header_index["SNP"]]
            strand = row[header_index["IlmnStrand"]]
            if (the_alleles == "[D/I]") or (the_alleles == "[I/D]"):
                # This is an indel, so we skip it
                continue
            top_alleles = illumina_to_snp(strand, the_alleles)
            allele1, allele2 = top_alleles.split(" ")
            alleles[marker_name] = {allele1, allele2}

    return alleles


def read_source_alleles(file_name):
    """Reads an allele file."""
    alleles = {}
    open_function = open
    if file_name.endswith(".gz"):
        open_function = gzip.open
    with open_function(file_name, 'rb') as input_file:
        for line in input_file:
            row = line.rstrip("\r\n").split("\t")
            marker_name = row[0]
            allele1, allele2 = row[1].split(" ")
            alleles[marker_name] = {allele1, allele2}

    return alleles


def illumina_to_snp(strand, snp):
    """Return the TOP strand of the marker.

    Function that takes a strand (TOP or BOT) and a SNP (e.g. : [A/C]) and
    returns a space separated AlleleA[space]AlleleB string.

    :param strand: Either "TOP" or "BOT"
    :type strand: str

    :param snp: [A/C], [A/T], [G/C], [T/C], [A/G], [C/G], [T/A] or [T/G].
    :type snp: str

    :returns: The nucleotide for allele A and the nucleotide for allele B
              (space separated)
    :rtype: str

    """

    # This correspondence matrix converts the alleles
    # from the annotation file to the appropriate
    # alleles for TOP strand.
    # This dictates most of the behavior for this script.
    # Modify it to customize allele conversion.
    illumina_correspondance = {
        "[A/C]": {
            "TOP": "A C",
        },
        "[A/G]": {
            "TOP": "A G",
        },
        "[T/C]": {
            "BOT": "A G",
        },
        "[T/G]": {
            "BOT": "A C",
        },
        "[A/T]": {
            "TOP": "A T",
            "BOT": "T A",
        },
        "[C/G]": {
            "TOP": "C G",
            "BOT": "G C",
        },
        "[G/C]": {
            "TOP": "G C",
            "BOT": "C G",
        },
        "[T/A]": {
            "TOP": "T A",
            "BOT": "A T",
        },
    }

    return illumina_correspondance[snp][strand]


def compute_statistics(out_dir, gold_prefix, source_prefix, same_samples,
                       use_sge, final_out_prefix):
    """Compute the statistics."""
    # Now, creating a temporary directory
    if not os.path.isdir(out_dir):
        try:
            os.mkdir(out_dir)
        except OSError:
            msg = "{}: file exists".format(out_dir)
            raise ProgramError(msg)

    # The out prefix
    out_prefix = os.path.join(out_dir, "tmp")

    # Subsetting the files
    logger.info("  - Subsetting the files")
    for k, (gold_sample, source_sample) in enumerate(same_samples):
        # Preparing the files
        try:
            filename = out_prefix + "_{}_gold.samples".format(k)
            with open(filename, 'w') as output_file:
                print >>output_file, "\t".join(gold_sample)

            filename = out_prefix + "_{}_source.samples".format(k)
            with open(filename, "w") as output_file:
                print >>output_file, "\t".join(source_sample)

        except IOError:
            msg = "can't write in dir {}".format(out_dir)
            raise ProgramError(msg)

    # Preparing the files
    nb = len(same_samples)
    keepSamples(
        [gold_prefix]*nb + [source_prefix]*nb,
        [out_prefix + "_{}_gold.samples".format(i) for i in range(nb)] +
        [out_prefix + "_{}_source.samples".format(i) for i in range(nb)],
        [out_prefix + "_{}_gold".format(i) for i in range(nb)] +
        [out_prefix + "_{}_source".format(i) for i in range(nb)],
        use_sge,
        transpose=True,
    )

    # Creating reports
    # The diff file
    diff_file = None
    try:
        diff_file = open(final_out_prefix + ".diff", 'w')
    except IOError:
        msg = "{}: can't write file".format(final_out_prefix + ".diff")
        raise ProgramError(msg)
    # The diff file header
    print >>diff_file, "\t".join(["name", "famID_source", "indID_source",
                                  "famID_gold", "indID_gold", "geno_source",
                                  "geno_gold"])

    # The concordance file
    concordance_file = None
    try:
        concordance_file = open(final_out_prefix + ".concordance", "w")
    except IOError:
        msg = "{}: can't write file".format(final_out_prefix + ".concordance")
        raise ProgramError(msg)
    # The concordance file header
    print >>concordance_file, "\t".join(["famID_source", "indID_source",
                                         "famID_gold", "indID_gold",
                                         "nb_geno", "nb_diff",
                                         "concordance"])

    for k in range(nb):
        # The samples
        gold_sample, source_sample = same_samples[k]

        # Reading the source freq file
        source_genotypes = {}
        filename = out_prefix + "_{}_source.tped".format(k)
        try:
            with open(filename, 'r') as input_file:
                for line in input_file:
                    row = line.rstrip("\r\n").split("\t")
                    marker_name = row[1]
                    geno = row[4]
                    if geno != "0 0":
                        geno = set(geno.split(" "))
                        source_genotypes[marker_name] = geno
        except IOError:
            msg = "{}: no such file".format(filename)
            raise ProgramError(msg)

        # Reading the gold freq file
        gold_genotypes = {}
        filename = out_prefix + "_{}_gold.tped".format(k)
        try:
            with open(filename, 'r') as input_file:
                for line in input_file:
                    row = line.rstrip("\r\n").split("\t")
                    marker_name = row[1]
                    geno = row[4]
                    if geno != "0 0":
                        geno = set(geno.split(" "))
                        gold_genotypes[marker_name] = geno
        except IOError:
            msg = "{}: no such file".format(filename)
            raise ProgramError(msg)

        # The genotyped markers in both, with their number
        genotyped_markers = source_genotypes.viewkeys() & \
            gold_genotypes.viewkeys()
        nb_genotyped = len(genotyped_markers)

        # Finding the number of differences, and creating the diff file
        nb_diff = 0
        for marker_name in genotyped_markers:
            # Getting the genotypes
            source_geno = source_genotypes[marker_name]
            gold_geno = gold_genotypes[marker_name]

            # Comparing the genotypes
            if source_geno != gold_geno:
                nb_diff += 1

                source_geno_print = list(source_geno - {"0"})
                source_geno_print.sort()
                if len(source_geno_print) == 1:
                    source_geno_print.append(source_geno_print[0])
                gold_geno_print = list(gold_geno - {"0"})
                gold_geno_print.sort()
                if len(gold_geno_print) == 1:
                    gold_geno_print.append(gold_geno_print[0])
                # We print in the diff file
                print >>diff_file, "\t".join([marker_name,
                                              "\t".join(source_sample),
                                              "\t".join(gold_sample),
                                              "/".join(source_geno_print),
                                              "/".join(gold_geno_print)])

        # Creating the concordance file
        concordance = 0.0
        if nb_genotyped != 0:
            concordance = (nb_genotyped - nb_diff) / float(nb_genotyped)
        print >>concordance_file, "\t".join(["\t".join(source_sample),
                                             "\t".join(gold_sample),
                                             str(nb_genotyped), str(nb_diff),
                                             str(concordance)])

    # Closing the output files
    diff_file.close()
    concordance_file.close()

    # Deleating the temporary directory
    try:
        shutil.rmtree(out_dir)
    except IOError:
        print >>sys.stderr, "   - Can't delete {}".format(out_dir)


def check_fam_for_samples(required_samples, source, gold):
    """Check fam files for required_samples."""
    # Checking the source panel
    source_samples = set()
    with open(source, 'r') as input_file:
        for line in input_file:
            sample = tuple(line.rstrip("\r\n").split(" ")[:2])
            if sample in required_samples:
                source_samples.add(sample)

    # Checking the gold standard
    gold_samples = set()
    with open(gold, 'r') as input_file:
        for line in input_file:
            sample = tuple(line.rstrip("\r\n").split(" ")[:2])
            if sample in required_samples:
                gold_samples.add(sample)

    # Checking if we got everything
    logger.info("  - Found {} samples in source panel".format(
        len(source_samples),
    ))
    logger.info("  - Found {} samples in gold standard".format(
        len(gold_samples),
    ))

    if len(required_samples - (source_samples | gold_samples)) != 0:
        return False
    else:
        return True


def read_same_samples_file(filename, out_prefix):
    """Reads a file containing same samples."""
    # The same samples
    same_samples = []

    # Creating the extraction files
    gold_file = None
    try:
        gold_file = open(out_prefix + ".gold_samples2keep", 'w')
    except IOError:
        msg = "{}: can't create file".format(out_prefix + ".gold_samples2keep")
        raise ProgramError(msg)

    source_file = None
    try:
        source_file = open(out_prefix + ".source_panel_samples2keep", 'w')
    except IOError:
        msg = ("{}: can't create "
               "file".format(out_prefix + ".source_panel_samples2keep"))
        raise ProgramError(msg)

    with open(filename, 'r') as input_file:
        for line in input_file:
            row = line.rstrip("\r\n").split("\t")

            # Getting the samples
            gold_sample = tuple(row[:2])
            source_sample = tuple(row[2:])
            same_samples.append((gold_sample, source_sample))

            # Printing files
            print >>gold_file, "\t".join(gold_sample)
            print >>source_file, "\t".join(source_sample)

    # Closing the files
    gold_file.close()
    source_file.close()

    return same_samples


def flipSNPs(inPrefix, outPrefix, flipFileName):
    """Flip SNPs using Plink."""
    plinkCommand = ["plink", "--noweb", "--bfile", inPrefix, "--flip",
                    flipFileName, "--make-bed", "--out", outPrefix]
    runCommand(plinkCommand)


def exclude_SNPs_samples(inPrefix, outPrefix, exclusionSNP=None,
                         keepSample=None, transpose=False):
    """Exclude some SNPs and keep some samples using Plink."""
    if (exclusionSNP is None) and (keepSample is None):
        msg = "Something wront with development... work on that source code..."
        raise ProgramError(msg)

    plinkCommand = ["plink", "--noweb", "--bfile", inPrefix, "--out",
                    outPrefix]
    if exclusionSNP is not None:
        # We want to exclude SNP
        plinkCommand.extend(["--exclude", exclusionSNP])
    if keepSample is not None:
        # We want to keep samples
        plinkCommand.extend(["--keep", keepSample])
    if transpose:
        # We want a transpose ped file
        plinkCommand.extend(["--recode", "--transpose", "--tab"])
    else:
        # We want a binary file
        plinkCommand.append("--make-bed")
    runCommand(plinkCommand)


def renameSNPs(inPrefix, updateFileName, outPrefix):
    """Updates the name of the SNPs using Plink."""
    plinkCommand = ["plink", "--noweb", "--bfile", inPrefix, "--update-map",
                    updateFileName, "--update-name", "--make-bed", "--out",
                    outPrefix]
    runCommand(plinkCommand)


def findFlippedSNPs(goldFrqFile1, sourceAlleles, outPrefix):
    """Find flipped SNPs and flip them in the data1."""
    goldAlleles = {}
    with open(goldFrqFile1, "r") as inputFile:
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

                goldAlleles[snpName] = alleles

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

    toRemoveOutputFileExplanation = None
    try:
        toRemoveOutputFileExplanation = open(
            outPrefix + ".snp_to_remove.explanation",
            "w",
        )
        print >>toRemoveOutputFileExplanation, "\t".join(["Name", "Reason",
                                                          "Alleles 1",
                                                          "Alleles 2"])
    except IOError:
        msg = "%(outPrefix)s.snp_to_remove: can't write file" % locals()
        raise ProgramError(msg)

    for snpName in goldAlleles.iterkeys():
        alleles1 = goldAlleles[snpName]
        alleles2 = sourceAlleles[snpName]

        if (len(alleles1) == 2) and (len(alleles2) == 2):
            # Both are heterozygous
            if (({"A", "T"} == alleles1 and {"A", "T"} == alleles2) or
                    ({"C", "G"} == alleles1 and {"C", "G"} == alleles2)):
                # We can't flip those..., so we remove them
                print >>toRemoveOutputFile, snpName
                print >>toRemoveOutputFileExplanation, "\t".join([
                    snpName,
                    "Undetermined",
                    "".join(alleles1),
                    "".join(alleles2),
                ])
            else:
                if alleles1 != alleles2:
                    # Let's try the flip one
                    if flipGenotype(alleles1) == alleles2:
                        # We need to flip it
                        print >>toFlipOutputFile, snpName
                    else:
                        # Those SNP are discordant...
                        print >>toRemoveOutputFile, snpName
                        print >>toRemoveOutputFileExplanation, "\t".join([
                            snpName,
                            "Invalid",
                            "".join(alleles1),
                            "".join(alleles2),
                        ])
        else:
            # We want to remove this SNP, because there is at least one
            # homozygous individual
            print >>toRemoveOutputFile, snpName
            tmp_allele1 = "".join(alleles1)
            if len(alleles1) == 1:
                tmp_allele1 += tmp_allele1
            tmp_allele2 = "".join(alleles1)
            if len(alleles1) == 1:
                tmp_allele2 += tmp_allele2
            print >>toRemoveOutputFileExplanation, "\t".join([snpName,
                                                              "Homozygous",
                                                              tmp_allele1,
                                                              tmp_allele2])

    # Closing output files
    toFlipOutputFile.close()
    toRemoveOutputFile.close()
    toRemoveOutputFileExplanation.close()


def computeFrequency(prefix, outPrefix):
    """Compute the frequency using Plink."""
    plinkCommand = ["plink", "--noweb", "--bfile", prefix, "--freq", "--out",
                    outPrefix]
    runCommand(plinkCommand)


def findOverlappingSNPsWithGoldStandard(prefix, gold_prefixe, out_prefix,
                                        use_marker_names=False):
    """Find the overlapping SNPs in 4 different data sets."""
    # Reading the main file
    sourceSnpToExtract = {}
    if use_marker_names:
        sourceSnpToExtract = set()
    duplicates = set()
    try:
        with open(prefix + ".bim", "r") as inputFile:
            for line in inputFile:
                row = line.rstrip("\r\n").split("\t")
                chromosome = row[0]
                position = row[3]
                snpName = row[1]

                if use_marker_names:
                    sourceSnpToExtract.add(snpName)
                else:
                    if (chromosome, position) not in sourceSnpToExtract:
                        sourceSnpToExtract[(chromosome, position)] = snpName
                    else:
                        # It's a duplicate
                        duplicates.add((chromosome, position))
    except IOError:
        msg = "%s.bim: no such file" % prefix
        raise ProgramError(msg)

    # Removing duplicates from the list
    if not use_marker_names:
        logger.info("  - There are {} duplicated markers "
                    "in {};".format(len(duplicates), prefix + ".bim"))
        logger.info("  - removing them for simplicity...")
        for snpID in duplicates:
            del sourceSnpToExtract[snpID]

    # Reading the Gold standard files
    goldSnpToExtract = {}
    if use_marker_names:
        goldSnpToExtract = set()
    with open(gold_prefixe + ".bim", "r") as inputFile:
        for line in inputFile:
            row = line.rstrip("\r\n").split("\t")
            chromosome = row[0]
            position = row[3]
            snpName = row[1]

            if use_marker_names:
                if snpName in sourceSnpToExtract:
                    goldSnpToExtract.add(snpName)
            else:
                if (chromosome, position) in sourceSnpToExtract:
                    # We want this SNP
                    goldSnpToExtract[(chromosome, position)] = snpName

    # Printing the names of the SNPs to extract
    goldOutputFile = None
    try:
        goldOutputFile = open(out_prefix + ".gold_snp_to_extract", "w")
    except IOError:
        msg = "%(out_prefix)s.goldSnpToExtract: can't write file" % locals()
        raise ProgramError(msg)

    sourceOutputFile = None
    try:
        sourceOutputFile = open(out_prefix + ".source_snp_to_extract", "w")
    except IOError:
        msg = "%(out_prefix)s.sourceSnpToExtract: can't write file" % locals()
        raise ProgramError(msg)

    changeNameOutputFile = None
    try:
        changeNameOutputFile = open(out_prefix + ".update_names", "w")
    except IOError:
        msg = "%(out_prefix)s.updateNames: can't write file" % locals()
        raise ProgramError(msg)

    # Writing the file
    if use_marker_names:
        for snpName in goldSnpToExtract:
            print >>sourceOutputFile, snpName
            print >>goldOutputFile, snpName
    else:
        for snpID in goldSnpToExtract.iterkeys():
            print >>sourceOutputFile, sourceSnpToExtract[snpID]
            print >>goldOutputFile, goldSnpToExtract[snpID]
            print >>changeNameOutputFile, "\t".join([
                goldSnpToExtract[snpID],
                sourceSnpToExtract[snpID],
            ])

    # Closing the output file
    goldOutputFile.close()
    sourceOutputFile.close()
    changeNameOutputFile.close()


def extractSNPs(prefixes, snpToExtractFileNames, outPrefixes, runSGE):
    """Extract a list of SNPs using Plink."""
    s = None
    jobIDs = []
    jobTemplates = []
    if runSGE:
        # Add the environment variable for DRMAA package
        if "DRMAA_LIBRARY_PATH" not in os.environ:
            t = "/shares/data/sge/lib/lx24-amd64/libdrmaa.so.1.0"
            os.environ["DRMAA_LIBRARY_PATH"] = t

        # Import the python drmaa library
        try:
            import drmaa
        except ImportError:
            raise ProgramError("drmaa is not install, install drmaa")

        # Initializing a session
        s = drmaa.Session()
        s.initialize()

    for k, prefix in enumerate(prefixes):
        plinkCommand = ["plink", "--noweb", "--bfile", prefix, "--extract",
                        snpToExtractFileNames[k], "--make-bed", "--out",
                        outPrefixes[k]]

        if runSGE:
            # We run using SGE
            # Creating the job template
            jt = s.createJobTemplate()
            jt.remoteCommand = plinkCommand[0]
            jt.workingDirectory = os.getcwd()
            jt.jobEnvironment = os.environ
            jt.args = plinkCommand[1:]
            jt.jobName = "_plink_extract_snps"

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
        # Deleating the job template, and exiting the session
        for jt in jobTemplates:
            s.deleteJobTemplate(jt)

        # Closing the connection
        s.exit()

        for hadProblem in hadProblems:
            if not hadProblem:
                msg = "Some SGE jobs had errors..."
                raise ProgramError(msg)


def keepSamples(prefixes, samplesToExtractFileNames, outPrefixes,
                runSGE, transpose=False):
    """Extract a list of SNPs using Plink."""
    s = None
    jobIDs = []
    jobTemplates = []
    if runSGE:
        # Add the environment variable for DRMAA package
        if "DRMAA_LIBRARY_PATH" not in os.environ:
            t = "/shares/data/sge/lib/lx24-amd64/libdrmaa.so.1.0"
            os.environ["DRMAA_LIBRARY_PATH"] = t

        # Import the python drmaa library
        import drmaa

        # Initializing a session
        s = drmaa.Session()
        s.initialize()

    for k, prefix in enumerate(prefixes):
        plinkCommand = ["plink", "--noweb", "--bfile", prefix, "--keep",
                        samplesToExtractFileNames[k], "--out", outPrefixes[k]]
        if transpose:
            plinkCommand.extend(["--recode", "--transpose", "--tab"])
        else:
            plinkCommand.append("--make-bed")

        if runSGE:
            # We run using SGE
            # Creating the job template
            jt = s.createJobTemplate()
            jt.remoteCommand = plinkCommand[0]
            jt.workingDirectory = os.getcwd()
            jt.jobEnvironment = os.environ
            jt.args = plinkCommand[1:]
            jt.jobName = "_plink_keep_samples"

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
        # Deleating the job template, and exiting the session
        for jt in jobTemplates:
            s.deleteJobTemplate(jt)

        # Closing the connection
        s.exit()

        for hadProblem in hadProblems:
            if not hadProblem:
                msg = "Some SGE jobs had errors..."
                raise ProgramError(msg)


def runCommand(command):
    """Run a command."""
    output = None
    try:
        output = subprocess.check_output(command,
                                         stderr=subprocess.STDOUT, shell=False)
    except subprocess.CalledProcessError:
        msg = "couldn't run command\n" + " ".join(command)
        raise ProgramError(msg)


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
    # Check if we have the binary files
    for prefix in [args.bfile, args.gold_bfile]:
        if prefix is None:
            msg = "no input file"
            raise ProgramError(msg)
        for fileName in [prefix + i for i in [".bed", ".bim", ".fam"]]:
            if not os.path.isfile(fileName):
                msg = "{}: no such file".format(fileName)
                raise ProgramError(msg)

    # Check for the same sample file
    if not os.path.isfile(args.same_samples):
        msg = "{}: no such file".format(args.same_samples)
        raise ProgramError(msg)

    # Check we have either a manifest or a allele file
    if (args.source_manifest is None) and (args.source_alleles is None):
        msg = ("need an allele file (either an Illumina manifest "
               "[--source-manifest] or a file containing alleles for each "
               "marker [--source-alleles]")
        raise ProgramError(msg)
    if ((args.source_manifest is not None) and
            (args.source_alleles is not None)):
        msg = ("use either --source-manifest or --source-alleles, not both")
        raise ProgramError(msg)
    # Check for the manifests
    if args.source_manifest is not None:
        if not os.path.isfile(args.source_manifest):
            msg = "{}: no such file".format(args.source_manifest)
            raise ProgramError(msg)
    if args.source_alleles is not None:
        if not os.path.isfile(args.source_alleles):
            msg = "{}: no such file".format(args.source_alleles)
            raise ProgramError(msg)

    return True


def parseArgs(argString=None):  # pragma: no cover
    """Parses the command line options and arguments.

    :returns: A :py:class:`numpy.Namespace` object created by
              the :py:mod:`argparse` module. It contains the values of the
              different options.

    ====================== ======= ============================================
        Options              Type                   Description
    ====================== ======= ============================================
    ``--bfile``            string  The input file prefix (will find the plink
                                   binary files by appending the prefix to the
                                   ``.bim``, ``.bed`` and ``.fam`` files,
                                   respectively).
    ``--gold-bfile``       string  The input file prefix (will find the plink
                                   binary files by appending the prefix to the
                                   ``.bim``, ``.bed`` and ``.fam`` files,
                                   respectively) for the Gold Standard.
    ``--same-samples``     string  A file containing samples which are present
                                   in both the gold standard and the source
                                   panel. One line by identity and tab
                                   separated. For each row, first sample is
                                   Gold Standard, second is source panel.
    ``--source-manifest``  string  The illumina marker manifest. This file
                                   should have tabs as field separator. There
                                   should be no lines before the main header
                                   line. There should be no lines after the
                                   last data line.
    ``--source-alleles``   string  A file containing the source alleles (TOP).
                                   Two columns (separated by tabulation, one
                                   with the marker name, the other with the
                                   alleles (separated by space). No header.
    ``--sge``              boolean Use SGE for parallelization.
    ``--do-not-flip``      boolean Do not flip SNPs. WARNING: only use this
                                   option only if the Gold Standard was
                                   generated using the same chip (hence,
                                   flipping is unnecessary).
    ``--use-marker-names`` boolean Use marker names instead of (chr, position).
                                   WARNING: only use this options only if the
                                   Gold Standard was generated using the same
                                   chip (hence, they have the same marker
                                   names).
    ``--out``              string  The prefix of the output files.
    ====================== ======= ============================================

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
pretty_name = "Gold standard comparison"
desc = "Compares the dataset with a gold standard."
long_desc = ("The script compares two dataset (the study dataset and a gold "
             "standard) to find genotyping discordance.")
parser = argparse.ArgumentParser(description=desc)
parser.add_argument("-v", "--version", action="version",
                    version="pyGenClean version {}".format(__version__))

# The INPUT files
group = parser.add_argument_group("Input File")
group.add_argument("--bfile", type=str, metavar="FILE", required=True,
                   help=("The input file prefix (will find the plink binary "
                         "files by appending the prefix to the .bim, .bed and "
                         ".fam files, respectively)."))
group.add_argument("--gold-bfile", type=str, metavar="FILE", required=True,
                   help=("The input file prefix (will find the plink binary "
                         "files by appending the prefix to the .bim, .bed and "
                         ".fam files, respectively) for the Gold Standard."))
group.add_argument("--same-samples", type=str, metavar="FILE", required=True,
                   help=("A file containing samples which are present in both "
                         "the gold standard and the source panel. One line by "
                         "identity and tab separated. For each row, first "
                         "sample is Gold Standard, second is source panel."))
group.add_argument("--source-manifest", type=str, metavar="FILE",
                   help=("The illumina marker manifest. This file should have "
                         "tabs as field separator. There should be no lines "
                         "before the main header line. There should be no "
                         "lines after the last data line."))
group.add_argument("--source-alleles", type=str, metavar="FILE",
                   help=("A file containing the source alleles (TOP). Two "
                         "columns (separated by tabulation, one with the "
                         "marker name, the other with the alleles (separated "
                         "by space). No header."))

# The options
group = parser.add_argument_group("Options")
group.add_argument("--sge", action="store_true",
                   help="Use SGE for parallelization.")
group.add_argument("--do-not-flip", action="store_true",
                   help=("Do not flip SNPs. WARNING: only use this option "
                         "only if the Gold Standard was generated using the "
                         "same chip (hence, flipping is unnecessary)."))
group.add_argument("--use-marker-names", action="store_true",
                   help=("Use marker names instead of (chr, position). "
                         "WARNING: only use this options only if the Gold "
                         "Standard was generated using the same chip (hence, "
                         "they have the same marker names)."))
# The OUTPUT files
group = parser.add_argument_group("Output File")
group.add_argument("--out", type=str, metavar="FILE",
                   default="compare_with_gold",
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
