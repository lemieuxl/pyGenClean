#!/usr/bin/env python2.7

import os
import sys
import argparse
import itertools
from collections import Counter

import numpy as npy

def main(argString=None):
    # Getting and checking the options
    args = parseArgs(argString)
    checkArgs(args)

    # Reading the tfiles
    tfams_info = read_tfams(args.tfiles)
    
    # Reading the tpeds
    read_tpeds(args.tfiles, tfams_info, args)


def encode_genotype(genotypes):
    """Removes spaces and sort alleles."""
    return "".join(sorted(genotypes.upper())[1:])


def read_tpeds(in_prefix, tfams, options):
    """Reads the input tpeds."""
    tfams, tfam_order = tfams

    # The consensus tfiles
    try:
        with open("{}.consensus.tfam".format(options.out), 'w') as output_file:
            print >>output_file, "\n".join(tfams[0][tfam_order[0]])
    except IOError:
        msg = "{}.consensus.tfam: can't write file".format(options.out)
        raise ProgramError(msg)
    tped_output_file = None
    marker_probs = None
    marker_summary = None
    try:
        tped_output_file = open("{}.consensus.tped".format(options.out), 'w')
        marker_probs = open("{}.probs".format(options.out), 'w')
        marker_summary = open("{}.marker_summary".format(options.out), "w")
        marker_cohen_kappa = open("{}.cohen_kappa".format(options.out), "w")
        marker_fleiss_kappa = open("{}.fleiss_kappa".format(options.out), 'w')
        marker_agreement = open("{}.percent_agreement".format(options.out), "w")
        with open("{}.info".format(options.out), "w") as output_file:
            for i, name in enumerate(in_prefix):
                print >>output_file, "\t".join([str(i+1), name])
    except IOError:
        msg = "can't write file output files"
        raise ProgramError(msg)

    # The number of samples
    nb_samples = len(tfams[0])
    nb_tpeds = len(in_prefix)
    nb_markers = 0
    probs_by_samples = npy.zeros(nb_samples, dtype=float)

    # We read all the tpeds, one line at a time and compare the values for each
    # genotypes. The tped should be sorted by marker IDs.
    opened_files = [open("{}.tped".format(i), 'rb') for i in in_prefix]

    # Vectorizing the function
    encode_genotype_array = npy.vectorize(encode_genotype)

    # We read all files one at a time
    for values in itertools.izip(*opened_files):
        # We format values
        values = npy.array([value.rstrip("\r\n").split("\t")
                                 for value in values])

        # Checking the integrety of the tped files
        if len(values.shape) != 2:
            # Problem with the tped files
            msg = "problem with the tped files (different number of samples?)"
            raise ProgramError(msg)

        # We check that the marker ID is the same at each line
        marker_info = values[0,:4]
        marker_info[2] = "0"
        for i in xrange(values.shape[0]):
            tmp = values[i,:4]
            tmp[2] = "0"
            if not npy.all(marker_info == tmp):
                print marker_info
                print tmp
                msg = "{}: not sorted".format(", ".join(in_prefix))
                raise ProgramError(msg)

        # Increasing the number of markers
        nb_markers += 1

        # We now check that we have the same amount of samples than tfams
        if values.shape[1] - 4 != nb_samples:
            msg = "tpeds and tfams have different number of samples"
            raise ProgramError(msg)

        # Extracting the genotypes, and ordering the alleles
        genotypes = encode_genotype_array(values[:,4:])
        for i in xrange(genotypes.shape[0]):
            genotypes[i] = genotypes[i][tfam_order[i]]

        # Creating the new genotypes
        new_genotypes = npy.array(["00" for i in xrange(nb_samples)])
        new_genotypes_probs = npy.array([0.0 for i in xrange(nb_samples)])
        for i in xrange(nb_samples):
            geno_counter = Counter(genotypes[:,i])

            # Getting the number of no calls, and remove them
            nb_no_calls = geno_counter["00"]

            # If the number of no calls are greater than 50%, we continue
            if npy.true_divide(nb_no_calls, nb_tpeds) > 0.5:
                continue

            # Removing the no calls from the counter
            del geno_counter["00"]

            # Getting the most common genotype
            most_common_genotype, nb_appearance = geno_counter.most_common(1)[0]

            # Computing the probs
            probs = 0.0
            if nb_tpeds - nb_no_calls != 0:
                probs = npy.true_divide(nb_appearance, nb_tpeds - nb_no_calls)
            if probs > 0.5:
                new_genotypes[i] = most_common_genotype
                new_genotypes_probs[i] = probs

        # Computing Cohen's Kappa
        cohen_kappa = []
        agreement = []
        for i in xrange(genotypes.shape[0]):
            for j in xrange(i+1, genotypes.shape[0]):
                geno_1 = genotypes[i]
                geno_2 = genotypes[j]
                p0 = npy.true_divide(npy.sum(geno_1 == geno_2), nb_samples)
                c1 = Counter(geno_1)
                c2 = Counter(geno_2)
                pe = 0
                for g in c1.viewkeys() | c2.viewkeys():
                    pe += npy.true_divide(c1[g], nb_samples) * npy.true_divide(c2[g], nb_samples)
                k = "nan"
                if pe != 1:
                    k = (p0 - pe) / (1 - pe)
                cohen_kappa.append(str(k))
                agreement.append(str(p0))

        # Computing Fleiss' Kappa
        p_a = 0
        p_j = Counter()
        for i in xrange(nb_samples):
            geno_counter = Counter(genotypes.T[i])
            p_i = npy.sum(npy.array(geno_counter.values()) ** 2) - nb_tpeds
            p_i /= nb_tpeds * (nb_tpeds - 1.0)
            p_a += p_i
            p_j += geno_counter
        p_a /= float(nb_samples)
        p_e = npy.sum(npy.true_divide(npy.array(p_j.values()), nb_tpeds * nb_samples) ** 2)
        fleiss_kappa = "nan"
        if p_e != 1:
            fleiss_kappa = (p_a - p_e) / (1.0 - p_e)

        # Printing the cohen_kappa, fleiss_kappa and agreement files
        print >>marker_cohen_kappa, "\t".join([marker_info[1]] + cohen_kappa)
        print >>marker_agreement, "\t".join([marker_info[1]] + agreement)
        print >>marker_fleiss_kappa, "\t".join([marker_info[1]] + [str(fleiss_kappa)])

##         # This is an example to compute
##         nb_samples = 9
##         geno_1 = npy.array(["AA", "AT", "TT", "TT", "AT", "TT", "AA", "AT", "TT"])
##         geno_2 = npy.array(["AA", "00", "AT", "TT", "AT", "AT", "AA", "TT", "TT"])
##         print geno_1
##         print geno_2
##         print geno_1 == geno_2, "->", npy.sum(geno_1 == geno_2)
##         p0 = npy.true_divide(npy.sum(geno_1 == geno_2), nb_samples)
##         c1 = Counter(geno_1)
##         c2 = Counter(geno_2)
##         pe = 0
##         for g in c1.viewkeys() | c2.viewkeys():
##             pe += npy.true_divide(c1[g], nb_samples) * npy.true_divide(c2[g], nb_samples)
##         k = (p0 - pe) / (1 - pe)
##         print k

        # Printing the new genotypes
        print >>tped_output_file, "\t".join(list(marker_info) +
                                            ["{} {}".format(i[0], i[-1])
                                                for i in list(new_genotypes)])

        # Printing the probs
        print >>marker_probs, "\t".join([marker_info[1]] +
                                          map(str, list(new_genotypes_probs)))

        # Printing the marker summary
        print >>marker_summary, "\t".join([marker_info[1],
                                           str(npy.mean(new_genotypes_probs))])

        # Saving the probs by samples
        probs_by_samples += new_genotypes_probs

    # Closing the opened files
    for opened_file in opened_files:
        opened_file.close()
    tped_output_file.close()
    marker_probs.close()
    marker_summary.close()
    marker_cohen_kappa.close()
    marker_fleiss_kappa.close()
    marker_agreement.close()

    # Computes the mean per samples
    probs_by_samples = npy.true_divide(probs_by_samples, nb_markers)
    try:
        with open("{}.sample_summary".format(options.out), 'w') as output_file:
            for sample, sample_probs in itertools.izip(tfams[0][tfam_order[0]],
                                                       probs_by_samples):
                print >>output_file, "\t".join(sample.split("\t")[:2] +
                                               [str(sample_probs)])
    except IOError:
        msg = "{}.sample_summary: can't write file".format(options.out)
        raise ProgramError(msg)


def read_tfams(in_prefix):
    """Reads the input tfams."""
    tfams = [None for i in xrange(len(in_prefix))]
    tfam_order = [None for i in xrange(len(in_prefix))]
    for prefix_nb, file_name in enumerate(["{}.tfam".format(i) for i in in_prefix]):
        data = None
        with open(file_name, 'rb') as input_file:
            # We know there are no header, since it's a tfam file, so we just
            # read it, keeping only the first two columns, separated by a space,
            # since we know there are no spaces in sample names
            data = npy.array(["\t".join(j)
                        for j in [k.rstrip("\r\n").split("\t")
                        for k in input_file.readlines()]])

        # Saving the sample names and the sorted order
        tfams[prefix_nb] = data
        tfam_order[prefix_nb] = npy.argsort(data)

    # We check we have the same amount of samples
    nb_samples = len(tfams[0])
    for tfam in tfams:
        if len(tfam) != nb_samples:
            msg = "{}: different number of samples".format(", ".join(in_prefix))
            raise ProgramError(msg)

    # Checking the tfams
    tfam = tfams[0][tfam_order[0]]
    for i in xrange(len(tfams)):
        if not npy.all(tfam == tfams[i][tfam_order[i]]):
            msg = "tfams doesn't have same value (after ordering)"
            raise ProgramError(msg)

    return (tfams, tfam_order)


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
    # Checking the input tfiles
    if len(args.tfiles) < 2:
        msg = "needs at least two tfiles to compare"
        raise ProgramError(msg)

    # Checking that the file exists
    for in_prefix in args.tfiles:
        for file_name in ["{}.{}".format(in_prefix, i) for i in ["tfam", "tped"]]:
            if not os.path.isfile(file_name):
                msg = "{}: no such file".format(file_name)
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
prog = "compare_markers"
desc = """Compare markers using tped files"""
parser = argparse.ArgumentParser(description=desc, prog=prog)

# The INPUT files
group = parser.add_argument_group("Input File")
group.add_argument("--tfiles", type=str, metavar="FILE", required=True,
                   nargs="+",
                   help=("The prefix of the tfiles (tfam and tped). All the "
                         "tped files should be ordered by marker IDs only "
                         "(which should be unique) and of the same sizes "
                         "(samples, and markers)."))
# The OUTPUT files
group = parser.add_argument_group("Output File")
group.add_argument("--out", type=str, metavar="FILE",
                    default="compared_markers",
                    help="The prefix of the output files. [default: " \
                         "%(default)s]")

if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print >>sys.stderr, "Cancelled by user"
        sys.exit(0)
    except ProgramError as e:
        parser.error(e.message)
