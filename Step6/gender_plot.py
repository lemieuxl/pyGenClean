#!/usr/bin/env python2.7

import os
import sys
import gzip
import argparse

import numpy as npy

desc = """Plots the gender using intensities"""
parser = argparse.ArgumentParser(description=desc)

def main():
    # Getting and checking the options
    args = parseArgs()
    checkArgs(args)

    data = None
    if args.summarized_intensities is None:
        # Reading the BIM file
        marker_names_chr = read_bim(args.bfile + ".bim")

        # Reading the FAM file
        sample_names_gender = read_fam(args.bfile + ".fam")

        # Reading the sex problem file
        sex_problems = read_sex_problems(args.sex_problems)

        # Reading the intensity file
        data = read_intensities(args.intensities, marker_names_chr,
                                sample_names_gender, sex_problems)

    else:
        data = read_summarized_intensities(args.summarized_intensities)

    # Plot the gender
    plot_gender(data, args)


def read_sex_problems(file_name):
    """Reads the sex problem file."""
    if file_name is None:
        return frozenset()

    problems = None
    with open(file_name, 'r') as input_file:
        header_index = dict([(col_name, i) for i, col_name
                                in enumerate(input_file.readline().rstrip("\r\n").split("\t"))])
        problems = frozenset([i.rstrip("\r\n").split("\t")[header_index["IID"]] for i
                                in input_file.readlines()])
    return problems


def encode_chr(chromosome):
    """Encodes chromosomes."""
    if chromosome == "X":
        return 23
    if chromosome == "Y":
        return 24
    if chromosome == "XY":
        return 25
    if chromosome == "MT":
        return 26
    try:
        new_chromosome = int(chromosome)
        if new_chromosome < 1 or new_chromosome > 26:
            msg = "{}: invalid chromosome".format(chromosome)
        return new_chromosome
    except ValueError:
        msg = "{}: invalid chromosome".format(chromosome)
        raise ProgramError(msg)


def encode_gender(gender):
    """Encodes the gender."""
    if gender == "1":
        return "Male"
    if gender == "2":
        return "Female"
    return "Unknown"


def read_bim(file_name):
    """Reads the BIM file to gather marker names."""
    marker_names_chr = None
    with open(file_name, 'r') as input_file:
        marker_names_chr = dict([(i[1], encode_chr(i[0])) for i
                                    in [j.rstrip("\r\n").split("\t") for j
                                    in input_file.readlines()]
                                    if encode_chr(i[0]) in {23, 24}])
    return marker_names_chr


def read_fam(file_name):
    """Reads the FAM file to gather sample names."""
    sample_names_gender = None
    with open(file_name, 'r') as input_file:
        sample_names_gender = dict([(i[1], encode_gender(i[4])) for i
                                        in [j.rstrip("\r\n").split(" ") for j
                                        in input_file.readlines()]])
    return sample_names_gender


def plot_gender(data, options):
    """Plots the gender."""
    if data is None:
        # there is a problem...
        msg = ("no data: specify either '--bfile' and '--intensities', or "
               "'--summarized-intensities'")
        raise ProgramError(msg)

    import matplotlib as mpl
    if options.format != "X11":
        mpl.use("Agg")
    import matplotlib.pyplot as plt
    if options.format != "X11":
        plt.ioff()

    # The figure and axes
    fig = plt.figure()
    fig.subplots_adjust(top=0.84)
    ax = fig.add_subplot(111)

    # Changing the spines
    ax.xaxis.set_ticks_position("bottom")
    ax.yaxis.set_ticks_position("left")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    # Setting the title
    ax.set_xlabel(options.xlabel)
    ax.set_ylabel(options.ylabel)

    # For the legend
    plot_object = []
    labels = []

    # Plotting the OK males
    males = npy.logical_and(data["gender"] == "Male", data["status"] == "OK")
    tmp, = ax.plot(data["chr23"][males], data["chr24"][males], "o", ms=5,
                 mec="#0099CC", mfc="#0099CC")
    plot_object.append(tmp)
    labels.append("OK Males (n={})".format(sum(males)))
    if options.summarized_intensities is None:
        print_data_to_file(data[males], "{}.ok_males.txt".format(options.out))

    # Plotting the OK females
    females = npy.logical_and(data["gender"] == "Female",
                              data["status"] == "OK")
    tmp, = ax.plot(data["chr23"][females], data["chr24"][females], "o", ms=5,
                   mec="#CC0000", mfc="#CC0000")
    plot_object.append(tmp)
    labels.append("OK Females (n={})".format(sum(females)))
    if options.summarized_intensities is None:
        print_data_to_file(data[females],
                           "{}.ok_females.txt".format(options.out))

    # Plotting the OK unknowns
    unknowns = npy.logical_and(data["gender"] == "Unknown",
                               data["status"] == "OK")
    tmp, = ax.plot(data["chr23"][unknowns], data["chr24"][unknowns], "o", ms=5,
                   mec="#555555", mfc="#555555")
    plot_object.append(tmp)
    labels.append("OK Unknowns (n={})".format(sum(unknowns)))
    if options.summarized_intensities is None:
        print_data_to_file(data[unknowns],
                           "{}.ok_unknowns.txt".format(options.out))

    # Plotting the Problem males
    males = npy.logical_and(data["gender"] == "Male",
                            data["status"] == "Problem")
    tmp, = ax.plot(data["chr23"][males], data["chr24"][males], "^", ms=6,
                   mec="#000000", mfc="#669900")
    plot_object.append(tmp)
    labels.append("Problematic Males (n={})".format(sum(males)))
    if options.summarized_intensities is None:
        print_data_to_file(data[males],
                           "{}.problematic_males.txt".format(options.out))

    # Plotting the Problem females
    females = npy.logical_and(data["gender"] == "Female",
                              data["status"] == "Problem")
    tmp, = ax.plot(data["chr23"][females], data["chr24"][females], "v", ms=6,
                   mec="#000000", mfc="#9933CC")
    plot_object.append(tmp)
    labels.append("Problematic Females (n={})".format(sum(females)))
    if options.summarized_intensities is None:
        print_data_to_file(data[females],
                           "{}.problematic_females.txt".format(options.out))

    # Plotting the Problem unknowns
    unknowns = npy.logical_and(data["gender"] == "Unknown",
                               data["status"] == "Problem")
    tmp, = ax.plot(data["chr23"][unknowns], data["chr24"][unknowns], ">", ms=6,
                   mec="#000000", mfc="#555555")
    plot_object.append(tmp)
    labels.append("Problematic Unknown (n={})".format(sum(unknowns)))
    if options.summarized_intensities is None:
        print_data_to_file(data[unknowns],
                           "{}.problematic_unknowns.txt".format(options.out))

    # the legend
    prop = mpl.font_manager.FontProperties(size=10)
    leg = ax.legend(plot_object, labels, loc=8, numpoints=1, fancybox=True,
                    prop=prop, ncol=2, bbox_to_anchor=(0., 1.02, 1., .102),
                    borderaxespad=0.)

    # Setting the limit
    xlim = ax.get_xlim()
    ax.set_xlim((xlim[0]-0.01, xlim[1]+0.01))
    ylim = ax.get_ylim()
    ax.set_ylim((ylim[0]-0.01, ylim[1]+0.01))

    if options.format == "X11":
        plt.show()
    else:
        file_name = "{}.{}".format(options.out, options.format)
        try:
            plt.savefig(file_name)
        except IOError:
            msg = "{}: can't write file".format(file_name)
            raise ProgramError(msg)


def print_data_to_file(data, file_name):
    """Prints data to file."""
    try:
        with open(file_name, 'w') as output_file:
            print >>output_file, "\t".join(data.dtype.names)
            for row in data:
                print >>output_file, "\t".join(map(str, row))
    except IOError:
        msg = "{}: can't write file".format(file_name)
        raise ProgramError(msg)

def read_summarized_intensities(prefix):
    """Reads the summarized intensities from 6 files."""
    data = []
    for suffix in {".ok_females.txt", ".ok_males.txt", ".ok_unknowns.txt",
                   ".problematic_females.txt", ".problematic_males.txt",
                   ".problematic_unknowns.txt"}:
        with open(prefix + suffix, 'r') as input_file:
            header_index = None
            for line_nb, line in enumerate(input_file):
                row = line.rstrip("\r\n").split("\t")

                if line_nb == 0:
                    # This is the header
                    header_index = dict([(col_name, i) for i, col_name
                                            in enumerate(row)])
                    for col_name in {"sampleID", "chr23", "chr24", "gender",
                                     "status"}:
                        if col_name not in header_index:
                            msg = "{}: no column named {}".format(prefix+suffix,
                                                                  col_name)
                            raise ProgramError(msg)

                else:
                    sampleID = row[header_index["sampleID"]]
                    chr23 = row[header_index["chr23"]]
                    chr24 = row[header_index["chr24"]]
                    gender = row[header_index["gender"]]
                    status = row[header_index["status"]]

                    try:
                        chr23 = float(chr23)
                        chr24 = float(chr24)
                    except ValueError:
                        msg = ("{} and {}: bad summarized intensities for "
                               "sample {}".format(chr23, chr24, sampleID))
                        raise ProgramError(msg)

                    data.append((sampleID, chr23, chr24, gender, status))

    # Creating the data structure
    data = npy.array(data, dtype=[("sampleID", "a{}".format(max([len(i[0]) for i
                                                        in data]))),
                                  ("chr23", float), ("chr24", float),
                                  ("gender", "a{}".format(max([len(i[3]) for i
                                                        in data]))),
                                  ("status", "a{}".format(max([len(i[4]) for i
                                                        in data])))])
    return data


def read_intensities(file_name, needed_markers_chr, needed_samples_gender,
                     sex_problems):
    """Reads the intensities from a file."""
    input_file = None
    if file_name.endswith(".gz"):
        input_file = gzip.open(file_name, 'rb')
    else:
        input_file = open(file_name, 'r')

    # The intensities
    intensities = {}

    header_index = None
    for line_nb, line in enumerate(input_file):
        row = line.rstrip("\r\n").split("\t")

        if line_nb == 0:
            # This is the header
            header_index = dict([(col_name, i) for i, col_name
                                    in enumerate(row)])
            # Check the column names
            for col_name in {"SNP Name", "Sample ID", "X", "Y", "Chr"}:
                if col_name not in header_index:
                    msg = "{}: no column named {}".format(file_name, col_name)
                    raise ProgramError(msg)
        else:
            # This is the data
            # Check if we want this sample and this marker
            sampleID = row[header_index["Sample ID"]]
            markerID = row[header_index["SNP Name"]]
            chromosome = encode_chr(row[header_index["Chr"]])
            if chromosome not in {23, 24}:
                # Not good chromsoome
                continue
            if sampleID not in needed_samples_gender:
                # Sample not needed
                continue
            if markerID not in needed_markers_chr:
                # Marker not needed
                continue

            if sampleID not in intensities:
                # First time we see this sample
                intensities[sampleID] = [0, 0, 0, 0]

            # We get the intensities
            allele_a = row[header_index["X"]]
            allele_b = row[header_index["Y"]]
            
            # Check for NaN
            if allele_a == "NaN" or allele_b == "NaN":
                continue
            try:
                allele_a = float(allele_a)
                allele_b = float(allele_b)
            except ValueError:
                msg = "{}: {} and {}: wrong intensities".format(file_name,
                                                                allele_a,
                                                                allele_b)
                raise ProgramError(msg)

            if chromosome == 23:
                # Chromsoome 23
                intensities[sampleID][0] += allele_a + allele_b
                intensities[sampleID][1] += 1
            else:
                # Chromosome 24
                intensities[sampleID][2] += allele_a + allele_b
                intensities[sampleID][3] += 1

    # Closing the input file
    input_file.close()

    # Creating the data structure
    data = []
    for sampleID in intensities.iterkeys():
        sum_chr_23, nb_chr_23, sum_chr_24, nb_chr_24 = intensities[sampleID]
        status = "OK"
        if sampleID in sex_problems:
            status = "Problem"
        try:
            data.append((sampleID, sum_chr_23 / nb_chr_23, sum_chr_24 /
                         nb_chr_24, needed_samples_gender[sampleID], status))
        except ZeroDivisionError:
            msg = "0 marker on chr23 or chr24"
            raise ProgramError(msg)

    data = npy.array(data, dtype=[("sampleID", "a{}".format(max([len(i[0]) for i
                                                        in data]))),
                                  ("chr23", float), ("chr24", float),
                                  ("gender", "a{}".format(max([len(i[3]) for i
                                                        in data]))),
                                  ("status", "a{}".format(max([len(i[4]) for i
                                                        in data])))])
    return data


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
    # Checking if there are input files
    if args.summarized_intensities is None:
        if args.bfile is None and args.intensities is None:
            msg = ("need to specify either '--bfile' and '--intensities', or "
                   "'--summarized-intensities'")
            raise ProgramError(msg)

        # Checking the fam file bim
        for suffix in {".bim", ".fam"}:
            if not os.path.isfile(args.bfile + suffix):
                msg = "{}: no such file".format(args.bfile + suffix)
                raise ProgramError(msg)

        # Checking the intensity file
        if not os.path.isfile(args.intensities):
            msg = "{}: no such file".format(args.intensities)
            raise ProgramError(msg)
    else:
        if args.bfile is not None or args.intensities is not None:
            msg = ("need to specify either '--bfile' and '--intensities', or "
                   "'--summarized-intensities'")
            raise ProgramError(msg)

        for suffix in {".ok_females.txt", ".ok_males.txt", ".ok_unknowns.txt",
                       ".problematic_females.txt", ".problematic_males.txt",
                       ".problematic_unknowns.txt"}:
            if not os.path.isfile(args.summarized_intensities + suffix):
                msg = "{}: no such file".format(args.summarized_intensities +
                                                suffix)
                raise ProgramError(msg)

    # Checking the sex problem file
    if args.sex_problems is not None:
        if not os.path.isfile(args.sex_problems):
            msg = "{}: no such file".format(args.sex_problems)
            raise ProgramError(msg)

    return True


def parseArgs(): # pragma: no cover
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
    # The INPUT files
    group = parser.add_argument_group("Input File")
    group.add_argument("--bfile", type=str, metavar="FILE",
                       help=("The plink binary file containing information "
                             "about markers and individuals. Must be specified "
                             "if '--summarized-intensities' is not."))
    group.add_argument("--intensities", type=str, metavar="FILE",
                       help=("A file containing alleles intensities for each "
                             "of the markers located on the X and Y "
                             "chromosome. Must be specified if "
                             "'--summarized-intensities' is not."))
    group.add_argument("--summarized-intensities", type=str, metavar="FILE",
                       help=("The prefix of six files (prefix.ok_females.txt, "
                             "prefix.ok_males.txt, prefix.ok_unknowns.txt, "
                             "problematic_females.txt, problematic_males.txt "
                             "and problematic_unknowns.txt) containing "
                             "summarized chr23 and chr24 intensities. Must be "
                             "specified if '--bfile' and '--intensities' are "
                             "not."))
    group.add_argument("--sex-problems", type=str, metavar="FILE",
                       help=("The file containing indivuduals with sex "
                             "problems. This file is not read if the option "
                             "'summarized-intensities' is used."))

    # The options
    group = parser.add_argument_group("Options")
    addCustomOptions(group)

    # The OUTPUT files
    group = parser.add_argument_group("Output File")
    group.add_argument("--out", type=str, metavar="FILE",
                       default="sexcheck",
                       help="The prefix of the output files (which will be " \
                            "a Plink binary file. [default: %(default)s]")

    args = parser.parse_args()

    return args


def addCustomOptions(parser):
    """Add custom options to a parser."""
    parser.add_argument("--format", type=str, metavar="FORMAT", default="png",
                        choices=["png", "ps", "pdf", "X11"],
                        help=("The output file format (png, ps, pdf, or X11 "
                              "formats are available). [default: %(default)s]"))
    parser.add_argument("--xlabel", type=str, metavar="STRING",
                        default="X intensity",
                        help="The label of the X axis. [default: %(default)s]")
    parser.add_argument("--ylabel", type=str, metavar="STRING",
                        default="Y intensity",
                        help="The label of the Y axis. [default: %(default)s]")


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


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print >>sys.stderr, "Cancelled by user"
        sys.exit(0)
    except ProgramError as e:
        parser.error(e.message)
