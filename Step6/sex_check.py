#!/usr/bin/env python2.7

import os
import sys
import argparse
import subprocess

import numpy as npy

from PlinkUtils import createRowFromPlinkSpacedOutput

def main(argString=None):
    # Getting and checking the options
    args = parseArgs(argString)
    checkArgs(args)

    print "   - Options used:"
    for key, value in vars(args).iteritems():
        print "      --{} {}".format(key, value)

    # Reads the bim file to see if chromosome 23 is there
    hasSexProblems = None
    if not checkBim("{}.bim".format(args.bfile), args.nbChr23, "23"):
        print ("   - There are not enough markers on chromosome 23\n"
               "     Stopping now!")
    else:
        # Run plink "check-sex"
        print "   - Running Plink for sex check"
        runPlinkSexCheck(args)

        # Reading plink "check-sex" output file
        print "   - Reading Plink's sex check output to find sex problems"
        hasSexProblems = readCheckSexFile(args.out + ".sexcheck",
                                        args.out + ".list_problem_sex",
                                        args.out + ".list_problem_sex_ids",
                                        args.femaleF, args.maleF)

    if hasSexProblems is not None and hasSexProblems:
        # Run plink to recode chr 23 in a ped format
        print "   - Creating recoded file for chr23 using Plink"
        createPedChr23UsingPlink(args)

        # Compute the hetero percentage
        print "   - Computing the heterozygous percentage"
        computeHeteroPercentage(args.out + ".chr23_recodeA.raw")

        # Run plink to get chr 24
        if checkBim("{}.bim".format(args.bfile), 1, "24"):
            print "   - Creating recoded file for chr24 using Plink"
            createPedChr24UsingPlink(args)

            # Compute the number of no call
            print "   - Computing the number of no calls"
            computeNoCall(args.out + ".chr24_recodeA.raw")
        else:
            print "   - Not enough markers on chr24"


def checkBim(fileName, minNumber, chromosome):
    """Checks the BIM file for chrN markers."""
    nbMarkers = 0
    with open(fileName, 'r') as inputFile:
        for line in inputFile:
            row = line.rstrip("\r\n").split("\t")
            if row[0] == chromosome:
                nbMarkers += 1
    if nbMarkers < minNumber:
        return False
    return True


def computeNoCall(fileName):
    """Computes the number of no call."""
    outputFile = None
    try:
        outputFile = open(fileName + ".noCall", "w")
    except IOError:
        msg = "%s: can't write file" % fileName + ".noCall"
        raise ProgramError(msg)
    print >>outputFile, "\t".join(["PED", "ID", "SEX", "nbGeno", "nbNoCall"])

    try:
        toPrint = []
        with open(fileName, "r") as inputFile:
            for i, line in enumerate(inputFile):
                row = line.rstrip("\r\n").split(" ")
                if i != 0:
                    # This is data
                    genotypes = npy.array(row[6:])
                    
                    nbMarker = len(genotypes)
                    nbNA = len(npy.where(genotypes == "NA")[0])

                    toPrint.append((row[0], row[1], row[4], str(nbMarker),
                                    str(nbNA)))

        toPrint.sort(reverse=True, key=lambda values: int(values[4]))

        for row in toPrint:
            print >>outputFile, "\t".join(row)

    except IOError:
        msg = "%(fileName)s: no such file" % locals()
        raise ProgramError(msg)

    # Closing the output file
    outputFile.close()


def computeHeteroPercentage(fileName):
    """Computes the heterozygosity percentage."""
    outputFile = None
    try:
        outputFile = open(fileName + ".hetero", "w")
    except IOError:
        msg = "%s: can't write file" % fileName + ".hetero"
        raise ProgramError(msg)
    print >>outputFile, "\t".join(["PED", "ID", "SEX", "HETERO"])

    try:
        toPrint = []
        with open(fileName, "r") as inputFile:
            for i, line in enumerate(inputFile):
                row = line.rstrip("\r\n").split(" ")
                if i != 0:
                    # This is data
                    genotypes = npy.array(row[6:])
                    genotypes = genotypes[npy.where(genotypes != "NA")]

                    nbMarker = len(genotypes)
                    nbHetero = len(npy.where(genotypes == "1")[0])
                    percentHetero = -9999
                    if nbMarker != 0:
                        percentHetero = nbHetero / float(nbMarker)

                    toPrint.append((row[0], row[1], row[4], percentHetero))

        # Sorting the data
        toPrint.sort(reverse=True, key=lambda values: values[3])

        # Printing the data
        for row in toPrint:
            value = list(row)
            if value[3] == -9999:
                value[3] = "ALL_NA"
            else:
                value[3] = str(value[3])
            print >>outputFile, "\t".join(value)

    except IOError:
        msg = "%(fileName)s: no such file" % locals()
        raise ProgramError(msg)

    # Closing the output file
    outputFile.close()


def readCheckSexFile(fileName, allProblemsFileName, idsFileName, femaleF,
                     maleF):
    """Reads the Plink check-sex output file."""
    allProblemsFile = None
    try:
        allProblemsFile = open(allProblemsFileName, 'w')
    except IOError:
        msg = "%(allProblemsFileName)s: can't write file" % locals()
        raise ProgramError(msg)

    idsFile = None
    try:
        idsFile = open(idsFileName, 'w')
    except IOError:
        msg = "%(idsFileName)s: can't write file" % locals()
        raise ProgramError(msg)

    try:
        with open(fileName, 'r') as inputFile:
            headerIndex = None
            nbProblems = 0
            nbTotalProblems = 0
            nbSexUnknown = 0
            nbFemaleThreshold = 0
            nbMaleThreshold = 0
            for line in inputFile:
                row = createRowFromPlinkSpacedOutput(line)

                if headerIndex is None:
                    # This is the header
                    headerIndex = dict([(row[i], i) \
                                            for i in xrange(len(row))])
                    for columnName in ["STATUS", "PEDSEX", "SNPSEX", "F",
                                       "FID", "IID"]:
                        if columnName not in headerIndex:
                            msg = "%(fileName)s: no culumn named " \
                                  "%(columnName)s" % locals()
                            raise ProgramError(msg)
                    print >>allProblemsFile, "\t".join(row)
                    continue

                # We have data
                status = row[headerIndex["STATUS"]]

                if status == "PROBLEM":
                    # We have a sex problem
                    nbTotalProblems += 1
                    pedsex = row[headerIndex["PEDSEX"]]

                    if pedsex == "0":
                        # The individual was "0", so we skip
                        nbSexUnknown += 1
                        continue

                    snpsex = row[headerIndex["SNPSEX"]]
                    if snpsex == "0":
                        # The new sex is unknown
                        f = None
                        try:
                            f = float(row[headerIndex["F"]])
                        except ValueError:
                            msg = "F=%s: not a float" % row[headerIndex["F"]]
                            raise ProgramError(msg)

                        if pedsex == "2":
                            # We have a female
                            if f < femaleF:
                                nbFemaleThreshold += 1
                                continue

                        if pedsex == "1":
                            # We have a male
                            if f > maleF:
                                nbMaleThreshold += 1
                                continue

                    print >>allProblemsFile, "\t".join(row)
                    
                    famID = row[headerIndex["FID"]]
                    indID = row[headerIndex["IID"]]
                    print >>idsFile, "\t".join([famID, indID])

                    nbProblems += 1
            
            print "   - Sex Check Summary"
            print "      - %(nbTotalProblems)d total problems" % locals()
            print "      - %(nbSexUnknown)d pedsex unknown" % locals()
            print "      - %(nbFemaleThreshold)d female F < %(femaleF)f" % locals()
            print "      - %(nbMaleThreshold)d male F > %(maleF)f" % locals()
            print "      - %(nbProblems)d problems kept" % locals()

    except IOError:
        msg = "%(fileName)s: no such file"

    if nbProblems == 0:
        # There are no sex problems to investigate
        print "   - There are no sex problem to investigate..."
        print "      - Nothing else to do..."
        return False
    return True

    # Closing the output files
    idsFile.close()
    allProblemsFile.close()


def runPlinkSexCheck(options):
    """Run Plink with the following options."""
    # The plink command
    plinkCommand = ["plink", "--noweb", "--bfile", options.bfile,
                    "--check-sex", "--out", options.out]
    runCommand(plinkCommand)


def createPedChr23UsingPlink(options):
    """Run Plink to create a ped format."""
    plinkCommand = ["plink", "--noweb", "--bfile", options.bfile, "--chr",
                    "23", "--recodeA", "--keep",
                    options.out + ".list_problem_sex_ids", "--out",
                    options.out + ".chr23_recodeA"]
    runCommand(plinkCommand)

def createPedChr24UsingPlink(options):
    """Run plink to create a ped format."""
    plinkCommand = ["plink", "--noweb", "--bfile", options.bfile, "--chr",
                    "24", "--recodeA", "--keep",
                    options.out + ".list_problem_sex_ids", "--out",
                    options.out + ".chr24_recodeA"]
    runCommand(plinkCommand)


def runCommand(command):
    """Run the command in Plink."""
    output = None
    try:
        output = subprocess.check_output(command,
                                         stderr=subprocess.STDOUT, shell=False)
    except subprocess.CalledProcessError:
        msg = "plink: couldn't run plink\n"
        msg += " ".join(command)
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

    # Ceck the number of markers on chromosome 23
    if args.nbChr23 < 0:
        msg = ("{}: number of markers on chr 23 must be "
               "positive".format(args.nbChr23))
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
prog = "sex_check"
desc = """Check sex using Plink"""
parser = argparse.ArgumentParser(description=desc, prog=prog)

# The INPUT files
group = parser.add_argument_group("Input File")
group.add_argument("--bfile", type=str, metavar="FILE", required=True,
                   help=("The input file prefix (will find the Plink binary "
                         "files by appending the prefix to the .bed, .bim, and "
                         ".fam files, respectively."))
# The options
group = parser.add_argument_group("Options")
group.add_argument("--femaleF", type=float, metavar="FLOAT", default=0.3,
                    help="The female F threshold. [default: < %(default)f]")
group.add_argument("--maleF", type=float, metavar="FLOAT", default=0.7,
                    help="The male F threshold. [default: > %(default)f]")
group.add_argument("--nbChr23", type=int, metavar="INT", default=50,
                    help=("The minimum number of markers on chromosome 23 "
                          "before computing Plink's sex check [default: "
                          "%(default)d]"))
# The OUTPUT files
group = parser.add_argument_group("Output File")
group.add_argument("--out", type=str, metavar="FILE", default="sexcheck",
                   help=("The prefix of the output files (which will be a "
                         "Plink binary file. [default: %(default)s]"))

if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print >>sys.stderr, "Cancelled by user"
        sys.exit(0)
    except ProgramError as e:
        parser.error(e.message)
