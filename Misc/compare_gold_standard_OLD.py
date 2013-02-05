#!/usr/bin/env python2.7
## This file is part of Foobar.
## 
## Foobar is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
## 
## Foobar is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
## 
## You should have received a copy of the GNU General Public License
## along with Foobar.  If not, see <http://www.gnu.org/licenses/>.

import re
import os
import sys
import shutil
import argparse
import subprocess

from PlinkUtils import createRowFromPlinkSpacedOutput
from StatGenDataCleanUp.Step2.duplicated_snps import flipGenotype

def main(argString=None):
    # Getting and checking the options
    args = parseArgs(argString)
    checkArgs(args)

    print "   - Options used:"
    for key, value in vars(args).iteritems():
        print "      --{} {}".format(key, value)

    # Reading the same sample file
    print "   - Reading the same samples file"
    same_samples = read_same_samples_file(args.same_samples, args.out)
    required_samples = set([i[0] for i in same_samples] +
                           [i[1] for i in same_samples])
    print ("      - Found a total of {} unique "
           "samples".format(len(required_samples)))

    # Check the fam files
    print "   - Checking in fam files for required samples"
    if not check_fam_for_samples(required_samples, args.bfile + ".fam",
                                 args.gold_bfile + ".fam"):
        print "      - Didn't find all required samples... STOP"
        sys.exit(0)

    # Find overlap markers with the gold standard file
    print "   - Finding overlapping SNPs with gold standard"
    findOverlappingSNPsWithGoldStandard(args.bfile, args.gold_bfile, args.out)

    # Extract the required SNPs using Plink
    print ("   - Extracting overlapping SNPs from the gold standard and the "
           "source panel")
    extractSNPs([args.gold_bfile, args.bfile],
                [args.out + ".gold_snp_to_extract",
                 args.out + ".source_snp_to_extract"],
                [args.out + ".gold_standard", args.out + ".source_panel"],
                args.sge)

    # Renaiming the reference file, so that the SNP names are the same
    print "   - Renaiming gold standard's SNPs to match source panel"
    renameSNPs(args.out + ".gold_standard", args.out + ".update_names",
               args.out + ".gold_standard.rename")

    # Computing the frequency
    print "   - Computing gold standard frequencies"
    computeFrequency(args.out + ".gold_standard.rename",
                     args.out + ".gold_standard.rename.frequency")
    print "   - Computing source panel frequencies"
    computeFrequency(args.out + ".source_panel",
                     args.out + ".source_panel.frequency")

    # Finding the SNPs to flip and flip them in the gold standard
    print "   - Finding SNPs to flip or to exclude from gold standard"
    findFlippedSNPs(args.out + ".gold_standard.rename.frequency.frq",
                    args.out + ".source_panel.frequency.frq",
                    args.out)

    # Excluding SNPs
    print "   - Excluding SNPs and samples from the gold standard"
    exclude_SNPs_samples(args.out + ".gold_standard.rename",
                args.out + ".gold_standard.rename.cleaned",
                args.out + ".snp_to_remove",
                args.out + ".gold_samples2keep")
    print "   - Excluding SNPs and samples from source panel"
    exclude_SNPs_samples(args.out + ".source_panel",
                         args.out + ".source_panel.cleaned",
                         args.out + ".snp_to_remove",
                         args.out + ".source_panel_samples2keep")

    # Flipping the SNP that need to be flip in the gold standard
    print "   - Flipping SNPs in gold standard"
    flipSNPs(args.out + ".gold_standard.rename.cleaned",
             args.out + ".gold_standard.rename.cleaned.flipped",
             args.out + ".snp_to_flip_in_reference")

    # Splitting the files, and running the duplicated samples script
    print "   - Computing statistics"
    compute_statistics(args.out + ".duplicated_samples",
                       args.out + ".gold_standard.rename.cleaned.flipped",
                       args.out + ".source_panel.cleaned", same_samples,
                       args.sge, args.out)


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
    print "      - Subsetting the files"
    for k, (gold_sample, source_sample) in enumerate(same_samples):
        # Preparing the files
        try:
            with open(out_prefix + "_{}_gold.samples".format(k), 'w') as output_file:
                print >>output_file, "\t".join(gold_sample)
            with open(out_prefix + "_{}_source.samples".format(k), "w") as output_file:
                print >>output_file, "\t".join(source_sample)
        except IOError:
            msg = "can't write in dir {}".format(out_dir)
            raise ProgramError(msg)

    # Preparing the files
    nb = len(same_samples)
    keepSamples([gold_prefix]*nb + [source_prefix]*nb,
                [out_prefix + "_{}_gold.samples".format(i) for i in range(nb)] +
                [out_prefix + "_{}_source.samples".format(i) for i in range(nb)],
                [out_prefix + "_{}_gold".format(i) for i in range(nb)] +
                [out_prefix + "_{}_source".format(i) for i in range(nb)],
                use_sge)

    # Computing the frequencies
    print "      - Computing frequencies of source panel"
    computeParallelFreq([out_prefix + "_{}_source".format(i) for i in range(nb)],
                        [out_prefix + "_{}_source.freq".format(i) for i in range(nb)],
                        use_sge)

    # Renaiming the names in the files
    for k in range(nb):
        for l, panel in enumerate(["gold", "source"]):
            filename = out_prefix + "_{}_{}.fam".format(k, panel)
            try:
                data = None
                with open(filename, 'r') as input_file:
                    data = input_file.readlines()
                with open(filename, 'w') as output_file:
                    for line in data:
                        output_file.write(re.sub(" ".join(same_samples[k][l]),
                                                "{} {}".format(k, k), line))
            except IOError:
                msg = "{}: can't process file".format(filename)
                raise ProgramError(msg)

    # Merging the files
    print "      - Finding differences and concordance"
    mergeBFiles([out_prefix + "_{}_gold".format(i) for i in range(nb)],
                [(j + "bed", j + "bim", j + "fam") for j in [out_prefix + "_{}_source.".format(i) for i in range(nb)]],
                [out_prefix + "_{}_merged".format(i) for i in range(nb)],
                use_sge)

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
                                         "concordance"])

    for k in range(nb):
        # The samples
        gold_sample, source_sample = same_samples[k]

        # First, we read the frequency file to get the number of genotyped
        # markers.
        filename = out_prefix + "_{}_source.freq.frq".format(k)
        nb_genotyped = 0
        try:
            with open(filename, 'r') as input_file:
                header_index = None
                for i, line in enumerate(input_file):
                    row = createRowFromPlinkSpacedOutput(line)
                    if i == 0:
                        # This is the header
                        header_index = dict([(col, i) for i, col in enumerate(row)])
                        for col_name in ["NCHROBS"]:
                            if col_name not in header_index:
                                msg = "{}: no column named {}".format(filename,
                                                                      col_name)
                                raise ProgramError(msg)
                    else:
                        if row[header_index["NCHROBS"]] != "0":
                            nb_genotyped += 1
        except IOError:
            msg = "{}: no such file".format(filename)
            raise ProgramError(msg)

        # Reading the diff file
        filename = out_prefix + "_{}_merged.diff".format(k)
        nb_diff = 0
        try:
            with open(filename, 'r') as input_file:
                header_index = None
                for i, line in enumerate(input_file):
                    row = createRowFromPlinkSpacedOutput(line)
                    if i == 0:
                        # This is the header
                        header_index = dict([(col, i) for i, col in enumerate(row)])
                        for col_name in ["SNP", "FID", "IID", "NEW", "OLD"]:
                            if col_name not in header_index:
                                msg = "{}: no column named {}".format(filename,
                                                                      col_name)
                                raise ProgramError(msg)
                    else:
                        # This is data
                        # We keep the SNP, and NEW and OLD genotype
                        print >>diff_file, "\t".join([row[header_index["SNP"]],
                                                      "\t".join(source_sample),
                                                      "\t".join(gold_sample),
                                                      row[header_index["NEW"]],
                                                      row[header_index["OLD"]]])
                        nb_diff += 1
        except IOError:
            msg = "{}: no such file".format(filename)

        # Creating the concordance file
        concordance = (nb_genotyped - nb_diff) / float(nb_genotyped)
        print >>concordance_file, "\t".join(["\t".join(source_sample),
                                             "\t".join(gold_sample),
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
    print "      - Found {} samples in source panel".format(len(source_samples))
    print "      - Found {} samples in gold standard".format(len(gold_samples))

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

    return same_samples


def flipSNPs(inPrefix, outPrefix, flipFileName):
    """Flip SNPs using Plink."""
    plinkCommand = ["plink", "--noweb", "--bfile", inPrefix, "--flip",
                    flipFileName, "--make-bed", "--out", outPrefix]
    runCommand(plinkCommand)


def exclude_SNPs_samples(inPrefix, outPrefix, exclusionSNP, keepSample,
                         transpose=False):
    """Exclude some SNPs and keep some samples using Plink."""
    plinkCommand = ["plink", "--noweb", "--bfile", inPrefix, "--exclude",
                    exclusionSNP, "--keep", keepSample, "--out", outPrefix]
    if transpose:
        plinkCommand.extend(["--recode", "--transpose", "--tab"])
    else:
        plinkCommand.append("--make-bed")
    runCommand(plinkCommand)


def renameSNPs(inPrefix, updateFileName, outPrefix):
    """Updates the name of the SNPs using Plink."""
    plinkCommand = ["plink", "--noweb", "--bfile", inPrefix, "--update-map",
                    updateFileName, "--update-name", "--make-bed", "--out",
                    outPrefix]
    runCommand(plinkCommand)


def findFlippedSNPs(frqFile1, frqFile2, outPrefix):
    """Find flipped SNPs and flip them in the data1."""
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
                        headerIndex = dict([(row[j], j) \
                                                for j in xrange(len(row))])

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
            if ({"A", "T"} == alleles1) or ({"C", "G"} == alleles1) or ({"A", "T"} == alleles2) or ({"C", "G"} == alleles2):
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
    """Compute the frequency using Plink."""
    plinkCommand = ["plink", "--noweb", "--bfile", prefix, "--freq", "--out",
                    outPrefix]
    runCommand(plinkCommand)


def combinePlinkBinaryFiles(prefixes, outPrefix):
    """Combine plink binary files."""
    # The first file is the bfile, the others are the ones to merge
    outputFile = None
    try:
        outputFile = open(outPrefix + ".files_to_merge", "w")
    except IOError:
        msg = "%(outPrefix)s.filesToMerge: can't write file" % locals()
        raise ProgramError(msg)

    for prefix in prefixes[1:]:
        print >>outputFile, " ".join([prefix + i for i in [".bed", ".bim", ".fam"]])

    # Closing the output files
    outputFile.close()

    # Runing plink
    plinkCommand = ["plink", "--noweb", "--bfile", prefixes[0],
                    "--merge-list", outPrefix + ".files_to_merge",
                    "--make-bed", "--out", outPrefix]
    runCommand(plinkCommand)


def findOverlappingSNPsWithGoldStandard(prefix, gold_prefixe, out_prefix):
    """Find the overlapping SNPs in 4 different data sets."""
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
    print ("      - There are {} duplicated markers "
           "in {};".format(len(duplicates), prefix + ".bim"))
    print "        removing them for simplicity..."
    for snpID in duplicates:
        del sourceSnpToExtract[snpID]

    # Reading the Gold standard files
    goldSnpToExtract = {}
    with open(gold_prefixe + ".bim", "r") as inputFile:
        for line in inputFile:
            row = line.rstrip("\r\n").split("\t")
            chromosome = row[0]
            position = row[3]
            snpName = row[1]

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
    for snpID in goldSnpToExtract.iterkeys():
        print >>sourceOutputFile, sourceSnpToExtract[snpID]
        print >>goldOutputFile, goldSnpToExtract[snpID]
        print >>changeNameOutputFile, "\t".join([goldSnpToExtract[snpID],
                                                 sourceSnpToExtract[snpID]])

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
            os.environ["DRMAA_LIBRARY_PATH"] = "/shares/data/sge/lib/lx24-amd64/libdrmaa.so.1.0"
        
        # Import the python drmaa library
        import drmaa

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


def mergeBFiles(prefixes, mergedPrefixes, outPrefixes, runSGE):
    """Extract a list of SNPs using Plink."""
    s = None
    jobIDs = []
    jobTemplates = []
    if runSGE:
        # Add the environment variable for DRMAA package
        if "DRMAA_LIBRARY_PATH" not in os.environ:
            os.environ["DRMAA_LIBRARY_PATH"] = "/shares/data/sge/lib/lx24-amd64/libdrmaa.so.1.0"
        
        # Import the python drmaa library
        import drmaa

        # Initializing a session
        s = drmaa.Session()
        s.initialize()

    for k, prefix in enumerate(prefixes):
        plinkCommand = ["plink", "--noweb", "--bfile", prefix, "--bmerge",
                        mergedPrefixes[k][0], mergedPrefixes[k][1],
                        mergedPrefixes[k][2], "--merge-mode", "7", "--out",
                        outPrefixes[k]]

        if runSGE:
            # We run using SGE
            # Creating the job template
            jt = s.createJobTemplate()
            jt.remoteCommand = plinkCommand[0]
            jt.workingDirectory = os.getcwd()
            jt.jobEnvironment = os.environ
            jt.args = plinkCommand[1:]
            jt.jobName = "_plink_compute_freq"

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


def computeParallelFreq(prefixes, outPrefixes, runSGE):
    """Extract a list of SNPs using Plink."""
    s = None
    jobIDs = []
    jobTemplates = []
    if runSGE:
        # Add the environment variable for DRMAA package
        if "DRMAA_LIBRARY_PATH" not in os.environ:
            os.environ["DRMAA_LIBRARY_PATH"] = "/shares/data/sge/lib/lx24-amd64/libdrmaa.so.1.0"
        
        # Import the python drmaa library
        import drmaa

        # Initializing a session
        s = drmaa.Session()
        s.initialize()

    for k, prefix in enumerate(prefixes):
        plinkCommand = ["plink", "--noweb", "--bfile", prefix, "--freq",
                        "--out", outPrefixes[k]]

        if runSGE:
            # We run using SGE
            # Creating the job template
            jt = s.createJobTemplate()
            jt.remoteCommand = plinkCommand[0]
            jt.workingDirectory = os.getcwd()
            jt.jobEnvironment = os.environ
            jt.args = plinkCommand[1:]
            jt.jobName = "_plink_compute_freq"

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


def keepSamples(prefixes, samplesToExtractFileNames, outPrefixes, runSGE):
    """Extract a list of SNPs using Plink."""
    s = None
    jobIDs = []
    jobTemplates = []
    if runSGE:
        # Add the environment variable for DRMAA package
        if "DRMAA_LIBRARY_PATH" not in os.environ:
            os.environ["DRMAA_LIBRARY_PATH"] = "/shares/data/sge/lib/lx24-amd64/libdrmaa.so.1.0"
        
        # Import the python drmaa library
        import drmaa

        # Initializing a session
        s = drmaa.Session()
        s.initialize()

    for k, prefix in enumerate(prefixes):
        plinkCommand = ["plink", "--noweb", "--bfile", prefix, "--keep",
                        samplesToExtractFileNames[k], "--make-bed", "--out",
                        outPrefixes[k]]

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
prog = "check_ethnicity"
desc = """Removes samples according to IBS."""
parser = argparse.ArgumentParser(description=desc, prog=prog)

# The INPUT files
group = parser.add_argument_group("Input File")
group.add_argument("--bfile", type=str, metavar="FILE", required=True,
                   help=("The input file prefix (will find the plink binary "
                         "files by appending the prefix to the .bim, .bed and "
                         ".fam files, respectively."))
group.add_argument("--gold-bfile", type=str, metavar="FILE", required=True,
                   help=("The input file prefix (will find the plink binary "
                         "files by appending the prefix to the .bim, .bed and "
                         ".fam files, respectively.) for the Gold Standard"))
group.add_argument("--same-samples", type=str, metavar="FILE", required=True,
                   help=("A file containing samples which are present in both "
                         "the gold standard and the source panel. One line by "
                         "identity and tab separated. For each row, first "
                         "sample is Gold Standard, second is source panel."))
# The options
group = parser.add_argument_group("Options")
group.add_argument("--sge", action="store_true",
                    help="Use SGE for parallelization.")
# The OUTPUT files
group = parser.add_argument_group("Output File")
group.add_argument("--out", type=str, metavar="FILE",
                   default="compare_with_gold",
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
