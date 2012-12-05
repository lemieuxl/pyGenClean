#!/usr/bin/env python2.7

import os
import sys
import glob
import gzip
import argparse
import subprocess
from collections import defaultdict

import numpy as npy
from PlinkUtils import createRowFromPlinkSpacedOutput
from StatGenDataCleanUp.Step9.merge_related_samples import merge_related_samples

def main(argString=None):
    # Getting and checking the options
    args = parseArgs(argString)
    checkArgs(args)

    print "   - Options used:"
    for key, value in vars(args).iteritems():
        print "      --{} {}".format(key, value)

    # Run plink
    print "   - Running Plink to extract SNPs according to LD"
    snpsToExtract = selectSNPsAccordingToLD(args)

    # Check there is enough SNP in the output file
    print "   - Checking if there are enough extracted SNP"
    if not checkNumberOfSNP(snpsToExtract, args.min_nb_snp):
        # There are not enough markers
        print "      - There are not enough SNPs"
        print "        Stopping now"

    else:
        # Remove the SNPs
        print "   - Extracting the SNPs using Plink"
        newBfile = extractSNPs(snpsToExtract, args)

        # Run the genome command from plink
        print "   - Creating the genome file using Plink"
        genomeFileName = runGenome(newBfile, args)

        if args.genome_only:
            # We just want the genome file
            return newBfile

        # Extract related individuals
        print "   - Finding related individuals from genome file"
        related_data = extractRelatedIndividuals(genomeFileName, args.out,
                                                 args.ibs2_ratio)

        # Plot the related data
        print "   - Plotting related individuals"
        plot_related_data(related_data["IBS2_RATIO"], related_data["Z1"],
                          related_data["CODE"], r"$Z_1$",
                          args.out + ".related_individuals_z1.png", args)
        plot_related_data(related_data["IBS2_RATIO"], related_data["Z2"],
                          related_data["CODE"], r"$Z_2$",
                          args.out + ".related_individuals_z2.png", args)


def plot_related_data(x, y, code, ylabel, fileName, options):
    """Plot Z1 and Z2 in function of IBS2* ratio."""
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
                    ["Full sibs (n={})".format(npy.sum(code == "1")),
                     ("Half sibs, grand-parent-child or uncle-nephew "
                      "(n={})".format(npy.sum(code == "2"))),
                     "Parent-child (n={})".format(npy.sum(code == "3")),
                     ("Twins or duplicated samples "
                      "(n={})".format(npy.sum(code == "4"))),
                     "Unknown (n={})".format(npy.sum(code == "5"))],
                    "best", numpoints=1, fancybox=True, prop=prop)
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
    """Extract related individuals according IBS2* ratio."""
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
            for columnName in ["FID1", "IID1", "FID2", "IID2", "Z0", "Z1", "Z2",
                               "HOMHOM", "HETHET"]:
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
                # Those paires might be related
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

                if (z0 >= 0.17 and z0 <= 0.33) and (z1 >=0.40 and z1 <= 0.60):
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

    # Merging the related individals
    merge_related_samples(outPrefix + ".related_individuals", outPrefix, False)

    # Creating the numpy array
    data = npy.array(data, dtype=[("IBS2_RATIO", float), ("Z1", float),
                                  ("Z2", float),
                                  ("CODE", "S{}".format(max([len(i[3]) for i in data])))])

    return data


def checkNumberOfSNP(fileName, minimumNumber):
    """Check there is enough SNPs in the file (with minimum)."""
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
    """Split a file."""
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
                        tmpFile = open(outPrefix + "_tmp.list%d" % nbTmpFile, "w")
                    except IOError:
                        msg = "tmp.list%d: can't write file" % nbTmpFile
                        raise ProgramError(msg)
                
                print >>tmpFile, " ".join(row[:2])

                if nbLine == linePerFile:
                    nbLine = 0
                    nbTmpFile += 1
                    tmpFile.close()
                    try:
                        tmpFile = open(outPrefix + "_tmp.list%d" % nbTmpFile, "w")
                    except IOError:
                        msg = "tmp.list%d: can't write file" % nbTmpFile
                        raise ProgramError(msg)
        tmpFile.close()
    except IOError:
        msg = "%s: no such file" % inputFileName
        raise ProgramError(msg)

    return nbTmpFile


def runGenome(bfile, options):
    """Run the genome command from plink."""
    outPrefix = options.out + ".genome"
    if options.sge:
        # We run genome using SGE
        # We need to create a frequency file using plink
        plinkCommand = ["plink", "--noweb", "--bfile", bfile, "--freq",
                        "--out", options.out + ".frequence"]
        runCommand(plinkCommand)

        # We need to split the .fam file
        nbJob = splitFile(bfile + ".fam", options.line_per_file_for_sge,
                          outPrefix)

        runGenomeSGE(bfile, options.out + ".frequence.frq", nbJob, outPrefix)

        # Merging genome files
        mergeGenomeLogFiles(outPrefix, nbJob)

    else:
        plinkCommand = ["plink", "--noweb", "--bfile", bfile, "--genome",
                        "--genome-full", "--out", outPrefix]
        runCommand(plinkCommand)

    return outPrefix + ".genome"


def mergeGenomeLogFiles(outPrefix, nbSet):
    """Merge genome and log files together."""
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

    # Removing the script
    try:
        os.remove(outPrefix + ".runPlinkGenomeSGE.sh")
    except IOError:
        msg = "runPlinkGenomeSGE.sh: can't delete"
        raise ProgramError(msg)

    # Closing the output files
    outputFile.close()
    outputLog.close()

    return outPrefix + ".genome"


def runGenomeSGE(bfile, freqFile, nbJob, outPrefix):
    """Runs the genome command from plink, on SGE."""
    # Add the environment variable for DRMAA package
    if "DRMAA_LIBRARY_PATH" not in os.environ:
        os.environ["DRMAA_LIBRARY_PATH"] = "/shares/data/sge/lib/lx24-amd64/libdrmaa.so.1.0"
    
    # Import the python drmaa library
    import drmaa

    # Print the script in a BASH file
    scriptFile = None
    try:
        scriptFile = open(outPrefix + ".runPlinkGenomeSGE.sh", "w")
    except IOError:
        msg = "runSGE.sh: can't write file"
        raise ProgramError(msg)
    print >>scriptFile, "#!/usr/bin/env bash"
    print >>scriptFile, "i=$SGE_TASK_ID"
    print >>scriptFile, "j=$SGE_TASK_ID"
    print >>scriptFile, "nbFile=%(nbJob)d" % locals()
    print >>scriptFile, ""
    print >>scriptFile, "while [ $i -le $nbFile ]"
    print >>scriptFile, "do"
    print >>scriptFile, "   plink \\"
    print >>scriptFile, "       --bfile %(bfile)s \\" % locals()
    print >>scriptFile, "       --read-freq %(freqFile)s \\" % locals()
    print >>scriptFile, "       --genome \\"
    print >>scriptFile, "       --genome-full \\"
    print >>scriptFile, "       --genome-lists %(outPrefix)s_tmp.list$j %(outPrefix)s_tmp.list$i \\" % locals()
    print >>scriptFile, "       --out %(outPrefix)s_output.sub.$j.$i" % locals()
    print >>scriptFile, "       i=$((i+1))"
    print >>scriptFile, "done"
    scriptFile.close()
    os.chmod(outPrefix + ".runPlinkGenomeSGE.sh", 0750)

    # Initializing a session
    s = drmaa.Session()
    s.initialize()

    # Creating the job template
    jt = s.createJobTemplate()
    jt.remoteCommand = outPrefix + ".runPlinkGenomeSGE.sh"
    jt.workingDirectory = os.getcwd()
    jt.jobEnvironment = os.environ
    jt.jobName = "_plink_genome"

    # Running the jobs and waiting for them
    jobList = s.runBulkJobs(jt, 1, nbJob, 1)
    s.synchronize(jobList, drmaa.Session.TIMEOUT_WAIT_FOREVER, False)
    jobFinished = True
    for curJob in jobList:
        retval = s.wait(curJob, drmaa.Session.TIMEOUT_WAIT_FOREVER)
        jobFinished = jobFinished and (retval.exitStatus == 0)

    # Deleating the job template, and exiting the session
    s.deleteJobTemplate(jt)
    s.exit()

    if not jobFinished:
        msg = "plink genome sge: can't run properly"
        raise ProgramError(msg)


def extractSNPs(snpsToExtract, options):
    """Remove SNPs using Plink."""
    outPrefix = options.out + ".pruned_data"
    plinkCommand = ["plink", "--noweb", "--bfile", options.bfile, "--extract",
                    snpsToExtract, "--make-bed", "--out", outPrefix]
    runCommand(plinkCommand)
    return outPrefix


def selectSNPsAccordingToLD(options):
    """Compute LD using Plink."""
    # The plink command
    outPrefix = options.out + ".pruning_" + options.indep_pairwise[2]
    plinkCommand = ["plink", "--noweb", "--bfile", options.bfile, "--maf",
                    options.maf, "--indep-pairwise"] + \
                    options.indep_pairwise + ["--out", outPrefix]

    runCommand(plinkCommand)

    return outPrefix + ".prune.in"


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
prog = "remove_IBS"
desc = """Removes samples according to IBS."""
parser = argparse.ArgumentParser(description=desc, prog=prog)

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
                   default=["50", "5", "0.1"],
                   help=("Three numbers: window size, window shift and the r2 "
                         "threshold. [default: %(default)s]"))
group.add_argument("--maf", type=str, metavar="FLOAT", default="0.05",
                   help=("Restrict to SNPs with MAF >= threshold. [default: "
                         "%(default)s]"))
group.add_argument("--ibs2-ratio", type=float, metavar="FLOAT", default=0.8,
                   help=("The initial IBS2* ratio (the minimum value to show "
                         "in the plot. [default: %(default).1f"))
group.add_argument("--sge", action="store_true",
                    help="Use SGE for parallelization.")
group.add_argument("--line-per-file-for-sge", type=int, metavar="INT",
                    default=100, help=("The number of line per file for SGE "
                                       "task array. [default: " "%(default)d]"))
# The OUTPUT files
group = parser.add_argument_group("Output File")
group.add_argument("--out", type=str, metavar="FILE", default="ibs",
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
