#!/usr/bin/env python2.7

import os
import sys
import argparse
import subprocess

import PlinkUtils.plot_MDS as PlotMDS
from PlinkUtils import createRowFromPlinkSpacedOutput
import StatGenDataCleanUp.Step9.remove_IBS as TheStep9
import StatGenDataCleanUp.Step10.find_outliers as find_outliers
from StatGenDataCleanUp.Step2.duplicated_snps import flipGenotype
 
class Dummy(object):
    pass

def main(argString=None):
    # Getting and checking the options
    args = parseArgs(argString)
    checkArgs(args)

    print "   - Options used:"
    for key, value in vars(args).iteritems():
        print "      --{} {}".format(key, value)

    # Find overlap with the reference file
    print "   - Finding overlapping SNPs between reference and source panels"
    referencePrefixes = [args.ceu_bfile, args.yri_bfile, args.jpt_chb_bfile]
    popNames = ["CEU", "YRI", "JPT-CHB"]
    findOverlappingSNPsWithReference(args.bfile, referencePrefixes, args.out)

    # Extract the required SNPs using Plink
    print "   - Extracting overlapping SNPs from the reference panels"
    extractSNPs(args.out + ".ref_snp_to_extract", referencePrefixes, popNames,
                args.out + ".reference_panel", args.sge)
    print "   - Extracting overlapping SNPs from the source panel"
    extractSNPs(args.out + ".source_snp_to_extract", [args.bfile], ["ALL"],
                args.out + ".source_panel", False)

    # Combining the reference panel
    print "   - Combining the reference panels"
    combinePlinkBinaryFiles([args.out + ".reference_panel." + i \
                                for i in popNames],
                            args.out + ".reference_panel.ALL")

    # Renaiming the reference file, so that the SNP names are the same
    print "   - Renaiming reference panel's SNPs to match source panel"
    renameSNPs(args.out + ".reference_panel.ALL", args.out + ".update_names",
               args.out + ".reference_panel.ALL.rename")

    # Computing the frequency
    print "   - Computing reference panel frequencies"
    computeFrequency(args.out + ".reference_panel.ALL.rename",
                     args.out + ".reference_panel.ALL.rename.frequency")
    print "   - Computing source panel frequencies"
    computeFrequency(args.out + ".source_panel.ALL",
                     args.out + ".source_panel.ALL.frequency")

    # Finding the SNPs to flip and flip them in the reference panel
    print "   - Finding SNPs to flip or to exclude from reference panel"
    findFlippedSNPs(args.out + ".reference_panel.ALL.rename.frequency.frq",
                    args.out + ".source_panel.ALL.frequency.frq",
                    args.out)

    # Excluding SNPs
    print "   - Excluding SNPs from reference panel"
    excludeSNPs(args.out + ".reference_panel.ALL.rename",
                args.out + ".reference_panel.ALL.rename.cleaned",
                args.out + ".snp_to_remove")
    print "   - Excluding SNPs from source panel"
    excludeSNPs(args.out + ".source_panel.ALL",
                args.out + ".source_panel.ALL.cleaned",
                args.out + ".snp_to_remove")

    # Flipping the SNP that need to be flip in the reference
    print "   - Flipping SNPs in reference panel"
    flipSNPs(args.out + ".reference_panel.ALL.rename.cleaned",
             args.out + ".reference_panel.ALL.rename.cleaned.flipped",
             args.out + ".snp_to_flip_in_reference")

    # Combining the reference panel
    print "   - Combining reference and source panels"
    combinePlinkBinaryFiles([args.out + ".reference_panel.ALL.rename.cleaned.flipped",
                             args.out + ".source_panel.ALL.cleaned"],
                            args.out + ".final_dataset_for_genome")

    # Runing the Step 9
    print "   - Creating the genome file using Plink"
    newBfile = runTheStep9(args.out + ".final_dataset_for_genome", args.out,
                           args)

    # Creating the MDS file
    print "   - Creating the MDS file using Plink"
    createMDSFile(args.nb_components, newBfile,
                  args.out + ".mds", args.out + ".ibs.genome.genome")

    # Creating the population files
    print "   - Creating a population file"
    famFiles = [args.out + ".reference_panel." + i + ".fam" for i in popNames]
    famFiles.append(args.out + ".source_panel.ALL.fam")
    labels = popNames + ["SOURCE"]
    createPopulationFile(famFiles, labels, args.out + ".population_file")

    # Plot the MDS value
    print "   - Creating the MDS plot"
    plotMDS(args.out + ".mds.mds", args.out + ".mds",
            args.out + ".population_file", args)

    # Finding the outliers
    find_the_outliers(args.out + ".mds.mds", args.out + ".population_file",
                      args.outliers_of, args.multiplier, args.out)


def find_the_outliers(mds_file_name, population_file_name, ref_pop_name,
                  multiplier, out_prefix):
    """Finds the outliers of a given population."""
    options = ["--mds", mds_file_name, "--population-file",
               population_file_name, "--outliers-of", ref_pop_name,
               "--multiplier", multiplier, "--out", out_prefix]

    try:
        find_outliers.main(options)
    except find_outliers.ProgramError as e:
        msg = "find_outliers: {}".format(e)
        raise ProgramError(msg)


def createPopulationFile(inputFiles, labels, outputFileName):
    """Creates a population file."""
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
    """Plots the MDS value."""
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
    """Creates a MDS file using Plink."""
    plinkCommand = ["plink", "--noweb", "--bfile", inPrefix, "--read-genome",
                    genomeFileName, "--cluster", "--mds-plot",
                    str(nb_components), "--out", outPrefix]
    runCommand(plinkCommand)

def runTheStep9(inputPrefix, outPrefix, options):
    """Run the Step9 of the data clean up."""
    # The options
    new_options = ["--bfile", inputPrefix, "--genome-only",
                   "--min-nb-snp", str(options.min_nb_snp),
                   "--maf", options.maf,
                   "--line-per-file-for-sge", str(options.line_per_file_for_sge),
                   "--out", "{}.ibs".format(outPrefix)]
    new_options += ["--indep-pairwise"] + options.indep_pairwise
    if options.sge:
        new_options.append("--sge")

    # Checking the input file
    if not allFileExists([inputPrefix + i for i in [".bed", ".bim", ".fam"]]):
        msg = "{}: not a valid binary prefix".format(inputPrefix)
        raise ProgramError(msg)

    newBfile = None
    try:
        newBfile = TheStep9.main(new_options)
    except TheStep9.ProgramError as e:
        msg = "compute genome: {}".format(e)
        raise ProgramError(msg)

    return newBfile


def allFileExists(fileList):
    """Check that all file exists."""
    allExists = True
    for fileName in fileList:
        allExists = allExists and os.path.isfile(fileName)
    return allExists


def flipSNPs(inPrefix, outPrefix, flipFileName):
    """Flip SNPs using Plink."""
    plinkCommand = ["plink", "--noweb", "--bfile", inPrefix, "--flip",
                    flipFileName, "--make-bed", "--out", outPrefix]
    runCommand(plinkCommand)


def excludeSNPs(inPrefix, outPrefix, exclusionFileName):
    """Exclude some SNPs using Plink."""
    plinkCommand = ["plink", "--noweb", "--bfile", inPrefix, "--exclude",
                    exclusionFileName, "--make-bed", "--out", outPrefix]
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


def findOverlappingSNPsWithReference(prefix, referencePrefixes, outPrefix):
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
    for snpID in duplicates:
        del sourceSnpToExtract[snpID]

    # Reading each of the
    refSnpToExtract = {}
    for refPrefix in referencePrefixes:
        try:
            with open(refPrefix + ".bim", "r") as inputFile:
                for line in inputFile:
                    row = line.rstrip("\r\n").split("\t")
                    chromosome = row[0]
                    position = row[3]
                    snpName = row[1]

                    if (chromosome, position) in sourceSnpToExtract:
                        # We want this SNP
                        refSnpToExtract[(chromosome, position)] = snpName
        except IOError:
            msg = "%(refPrefix)s.bim: no such file" % locals()
            raise ProgramError(msg)

    # Printing the names of the SNPs to extract
    refOutputFile = None
    try:
        refOutputFile = open(outPrefix + ".ref_snp_to_extract", "w")
    except IOError:
        msg = "%(outPrefix)s.refSnpToExtract: can't write file" % locals()
        raise ProgramError(msg)

    sourceOutputFile = None
    try:
        sourceOutputFile = open(outPrefix + ".source_snp_to_extract", "w")
    except IOError:
        msg = "%(outPrefix)s.sourceSnpToExtract: can't write file" % locals()
        raise ProgramError(msg)

    changeNameOutputFile = None
    try:
        changeNameOutputFile = open(outPrefix + ".update_names", "w")
    except IOError:
        msg = "%(outPrefix)s.updateNames: can't write file" % locals()
        raise ProgramError(msg)

    # Writing the file
    for snpID in refSnpToExtract.iterkeys():
        print >>sourceOutputFile, sourceSnpToExtract[snpID]
        print >>refOutputFile, refSnpToExtract[snpID]
        print >>changeNameOutputFile, "\t".join([refSnpToExtract[snpID],
                                                 sourceSnpToExtract[snpID]])

    # Closing the output file
    refOutputFile.close()
    sourceOutputFile.close()
    changeNameOutputFile.close()


def extractSNPs(snpToExtractFileName, referencePrefixes, popNames, outPrefix,
                runSGE):
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


    for k, refPrefix in enumerate(referencePrefixes):
        plinkCommand = ["plink", "--noweb", "--bfile", refPrefix, "--extract",
                        snpToExtractFileName, "--make-bed", "--out",
                        outPrefix + "." + popNames[k]]

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
    # Check if we have the tped and the tfam files
    for prefix in [args.bfile, args.ceu_bfile, args.yri_bfile,
                   args.jpt_chb_bfile]:
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
group.add_argument("--ceu-bfile", type=str, metavar="FILE", required=True,
                   help=("The input file prefix (will find the plink binary "
                         "files by appending the prefix to the .bim, .bed and "
                         ".fam files, respectively.) for the CEU population"))
group.add_argument("--yri-bfile", type=str, metavar="FILE", required=True,
                   help=("The input file prefix (will find the plink binary "
                         "files by appending the prefix to the .bim, .bed and "
                         ".fam files, respectively.) for the YRI population"))
group.add_argument("--jpt-chb-bfile", type=str, metavar="FILE", required=True,
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
group.add_argument("--line-per-file-for-sge", type=int, metavar="INT",
                    default=100, help=("The number of line per file for "
                                       "SGE task array. [default: "
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
# The OUTPUT files
group = parser.add_argument_group("Output File")
group.add_argument("--out", type=str, metavar="FILE", default="ethnic",
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
