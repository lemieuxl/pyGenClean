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
import gzip
import logging
import argparse

import numpy as np

from .. import __version__


logger = logging.getLogger("gender_plot")


def main(argString=None):
    """The main function of the module.

    :param argString: the options.

    :type argString: list

    These are the steps:

    1. Prints the options.
    2. If there are ``summarized_intensities`` provided, reads the files
       (:py:func:`read_summarized_intensities`) and skips to step 7.
    3. Reads the ``bim`` file to get markers on the sexual chromosomes
       (:py:func:`read_bim`).
    4. Reads the ``fam`` file to get gender (:py:func:`read_fam`).
    5. Reads the file containing samples with sex problems
       (:py:func:`read_sex_problems`).
    6. Reads the intensities and summarizes them (:py:func:`read_intensities`).
    7. Plots the summarized intensities (:py:func:`plot_gender`).

    """
    # Getting and checking the options
    args = parseArgs(argString)
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
    """Reads the sex problem file.

    :param file_name: the name of the file containing sex problems.

    :type file_name: str

    :returns: a :py:class:`frozenset` containing samples with sex problem.

    If there is no ``file_name`` (*i.e.* is ``None``), then an empty
    :py:class:`frozenset` is returned.

    """
    if file_name is None:
        return frozenset()

    problems = None
    with open(file_name, 'r') as input_file:
        header_index = dict([
            (col_name, i) for i, col_name in
            enumerate(input_file.readline().rstrip("\r\n").split("\t"))
        ])
        if "IID" not in header_index:
            msg = "{}: no column named IID".format(file_name)
            raise ProgramError(msg)

        problems = frozenset([
            i.rstrip("\r\n").split("\t")[header_index["IID"]]
            for i in input_file.readlines()
        ])
    return problems


def encode_chr(chromosome):
    """Encodes chromosomes.

    :param chromosome: the chromosome to encode.

    :type chromosome: str

    :returns: the encoded chromosome as :py:class:`int`.

    It changes ``X``, ``Y``, ``XY`` and ``MT`` to ``23``, ``24``, ``25`` and
    ``26``, respectively. It changes everything else as :py:class:`int`.

    If :py:class:`ValueError` is raised, then :py:class:`ProgramError` is
    also raised. If a chromosome as a value below 1 or above 26, a
    :py:class:`ProgramError` is raised.

    .. testsetup::

        from pyGenClean.SexCheck.gender_plot import encode_chr

    .. doctest::

        >>> [encode_chr(str(i)) for i in range(0, 11)]
        [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
        >>> [encode_chr(str(i)) for i in range(11, 21)]
        [11, 12, 13, 14, 15, 16, 17, 18, 19, 20]
        >>> [encode_chr(str(i)) for i in range(21, 27)]
        [21, 22, 23, 24, 25, 26]
        >>> [encode_chr(i) for i in ["X", "Y", "XY", "MT"]]
        [23, 24, 25, 26]
        >>> encode_chr("27")
        Traceback (most recent call last):
            ...
        ProgramError: 27: invalid chromosome
        >>> encode_chr("XX")
        Traceback (most recent call last):
            ...
        ProgramError: XX: invalid chromosome

    """
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
        if new_chromosome < 0 or new_chromosome > 26:
            msg = "{}: invalid chromosome".format(chromosome)
            raise ProgramError(msg)
        return new_chromosome
    except ValueError:
        msg = "{}: invalid chromosome".format(chromosome)
        raise ProgramError(msg)


def encode_gender(gender):
    """Encodes the gender.

    :param gender: the gender to encode.

    :type gender: str

    :returns: the encoded gender.

    It changes ``1`` and ``2`` to ``Male`` and ``Female`` respectively. It
    encodes everything else to ``Unknown``.

    .. testsetup::

        from pyGenClean.SexCheck.gender_plot import encode_gender

    .. doctest::

        >>> encode_gender("1")
        'Male'
        >>> encode_gender("2")
        'Female'
        >>> encode_gender("0")
        'Unknown'
        >>> encode_gender("This is not a gender code")
        'Unknown'

    """
    if gender == "1":
        return "Male"
    if gender == "2":
        return "Female"
    return "Unknown"


def read_bim(file_name):
    """Reads the BIM file to gather marker names.

    :param file_name: the name of the ``bim`` file.

    :type file_name: str

    :returns: a :py:class:`dict` containing the chromosomal location of each
              marker on the sexual chromosomes.

    It uses the :py:func:`encode_chr` to encode the chromosomes from ``X`` and
    ``Y`` to ``23`` and ``24``, respectively.

    """
    marker_names_chr = None
    with open(file_name, 'r') as input_file:
        marker_names_chr = dict([
            (i[1], encode_chr(i[0]))
            for i in [
                j.rstrip("\r\n").split("\t") for j in input_file.readlines()
            ] if encode_chr(i[0]) in {23, 24}
        ])
    return marker_names_chr


def read_fam(file_name):
    """Reads the FAM file to gather sample names.

    :param file_name: the ``fam`` file to read.

    :type file_name: str

    :returns: a :py:class:`dict` containing the gender of each samples.

    It uses the :py:func:`encode_gender` to encode the gender from ``1``and
    ``2`` to ``Male`` and ``Female``, respectively.

    """
    sample_names_gender = None
    with open(file_name, 'r') as input_file:
        sample_names_gender = dict([
            (i[1], encode_gender(i[4])) for i in [
                j.rstrip("\r\n").split(" ") for j in input_file.readlines()
            ]
        ])
    return sample_names_gender


def plot_gender(data, options):
    """Plots the gender.

    :param data: the data to plot.
    :param options: the options.

    :type data: numpy.recarray
    :type options: argparse.Namespace

    Plots the summarized intensities of the markers on the Y chromosomes in
    function of the markers on the X chromosomes, with problematic samples with
    different colors.

    Also uses :py:func:`print_data_to_file` to save the data, so that it is
    faster to rerun the analysis.

    """
    if data is None:
        # there is a problem...
        msg = ("no data: specify either '--bfile' and '--intensities', or "
               "'--summarized-intensities'")
        raise ProgramError(msg)

    import matplotlib as mpl
    if options.format != "X11" and mpl.get_backend() != "agg":
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
    males = np.logical_and(data["gender"] == "Male", data["status"] == "OK")
    tmp, = ax.plot(data["chr23"][males], data["chr24"][males], "o", ms=5,
                   mec="#0099CC", mfc="#0099CC")
    plot_object.append(tmp)
    labels.append("OK Males (n={})".format(sum(males)))
    if options.summarized_intensities is None:
        print_data_to_file(data[males], "{}.ok_males.txt".format(options.out))

    # Plotting the OK females
    females = np.logical_and(data["gender"] == "Female",
                             data["status"] == "OK")
    tmp, = ax.plot(data["chr23"][females], data["chr24"][females], "o", ms=5,
                   mec="#CC0000", mfc="#CC0000")
    plot_object.append(tmp)
    labels.append("OK Females (n={})".format(sum(females)))
    if options.summarized_intensities is None:
        print_data_to_file(data[females],
                           "{}.ok_females.txt".format(options.out))

    # Plotting the OK unknowns
    unknowns = np.logical_and(data["gender"] == "Unknown",
                              data["status"] == "OK")
    tmp, = ax.plot(data["chr23"][unknowns], data["chr24"][unknowns], "o", ms=5,
                   mec="#555555", mfc="#555555")
    plot_object.append(tmp)
    labels.append("OK Unknowns (n={})".format(sum(unknowns)))
    if options.summarized_intensities is None:
        print_data_to_file(data[unknowns],
                           "{}.ok_unknowns.txt".format(options.out))

    # Plotting the Problem males
    males = np.logical_and(data["gender"] == "Male",
                           data["status"] == "Problem")
    tmp, = ax.plot(data["chr23"][males], data["chr24"][males], "^", ms=6,
                   mec="#000000", mfc="#669900")
    plot_object.append(tmp)
    labels.append("Problematic Males (n={})".format(sum(males)))
    if options.summarized_intensities is None:
        print_data_to_file(data[males],
                           "{}.problematic_males.txt".format(options.out))

    # Plotting the Problem females
    females = np.logical_and(data["gender"] == "Female",
                             data["status"] == "Problem")
    tmp, = ax.plot(data["chr23"][females], data["chr24"][females], "v", ms=6,
                   mec="#000000", mfc="#9933CC")
    plot_object.append(tmp)
    labels.append("Problematic Females (n={})".format(sum(females)))
    if options.summarized_intensities is None:
        print_data_to_file(data[females],
                           "{}.problematic_females.txt".format(options.out))

    # Plotting the Problem unknowns
    unknowns = np.logical_and(data["gender"] == "Unknown",
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
    """Prints data to file.

    :param data: the data to print.
    :param file_name: the name of the output file.

    :type data: numpy.recarray
    :type file_name: str

    """
    try:
        with open(file_name, 'w') as output_file:
            print >>output_file, "\t".join(data.dtype.names)
            for row in data:
                print >>output_file, "\t".join(map(str, row))
    except IOError:
        msg = "{}: can't write file".format(file_name)
        raise ProgramError(msg)


def read_summarized_intensities(prefix):
    """Reads the summarized intensities from 6 files.

    :param prefix: the prefix of the six files.

    :type prefix: str

    :returns: a :py:class`numpy.recarray` containing the following columns (for
              each of the samples): ``sampleID``, ``chr23``, ``chr24``,
              ``gender`` and ``status``.

    Instead of reading a final report (like :py:func:`read_intensities`), this
    function reads six files previously created by this module to gather sample
    information. Here are the content of the six files:

    * ``prefix.ok_females.txt``: information about females without sex problem.
    * ``prefix.ok_males.txt``: information about males without sex problem.
    * ``prefix.ok_unknowns.txt``: information about unknown gender without sex
                                  problem.
    * ``prefix.problematic_females.txt``: information about females with sex
                                          problem.
    * ``prefix.problematic_males.txt``: information about males with sex
                                        problem.
    * ``prefix.problematic_unknowns.txt``: information about unknown gender
                                           with sex problem.

    Each file contains the following columns: ``sampleID``, ``chr23``,
    ``chr24``, ``gender`` and ``status``.

    The final data set contains the following information for each sample:

    * ``sampleID``: the sample ID.
    * ``chr23``: the summarized intensities for chromosome 23.
    * ``chr24``: the summarized intensities for chromosome 24.
    * ``gender``: the gender of the sample (``Male`` or ``Female``).
    * ``status``: the status of the sample (``OK`` or ``Problem``).

    The summarized intensities for a chromosome (:math:`S_{chr}`) is computed
    using this formula (where :math:`I_{chr}` is the set of all marker
    intensities on chromosome :math:`chr`):

    .. math::
        S_{chr} = \\frac{\\sum{I_{chr}}}{||I_{chr}||}

    """
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
                    header_index = dict([
                        (col_name, i) for i, col_name in enumerate(row)
                    ])
                    for col_name in {"sampleID", "chr23", "chr24", "gender",
                                     "status"}:
                        if col_name not in header_index:
                            msg = "{}: no column named {}".format(
                                prefix+suffix,
                                col_name,
                            )
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
    data = np.array(
        data,
        dtype=[("sampleID", "a{}".format(max([len(i[0]) for i in data]))),
               ("chr23", float), ("chr24", float),
               ("gender", "a{}".format(max([len(i[3]) for i in data]))),
               ("status", "a{}".format(max([len(i[4]) for i in data])))],
    )
    return data


def read_intensities(file_name, needed_markers_chr, needed_samples_gender,
                     sex_problems):
    """Reads the intensities from a file.

    :param file_name: the name of the input file.
    :param needed_markers_chr: the markers that are needed.
    :param needed_samples_gender: the gender of all the samples.
    :param sex_problems: the sample with sex problem.

    :type file_name: str
    :type needed_markers_chr: dict
    :type needed_samples_gender: dict
    :type sex_problems: frozenset

    :returns: a :py:class`numpy.recarray` containing the following columns (for
              each of the samples): ``sampleID``, ``chr23``, ``chr24``,
              ``gender`` and ``status``.

    Reads the normalized intensities from a final report. The file must contain
    the following columns: ``SNP Name``, ``Sample ID``, ``X``, ``Y`` and
    ``Chr``. It then keeps only the required markers (those that are on
    sexual chromosomes (``23`` or ``24``), encoding `NaN` intensities to zero.

    The final data set contains the following information for each sample:

    * ``sampleID``: the sample ID.
    * ``chr23``: the summarized intensities for chromosome 23.
    * ``chr24``: the summarized intensities for chromosome 24.
    * ``gender``: the gender of the sample (``Male`` or ``Female``).
    * ``status``: the status of the sample (``OK`` or ``Problem``).

    The summarized intensities for a chromosome (:math:`S_{chr}`) is computed
    using this formula (where :math:`I_{chr}` is the set of all marker
    intensities on chromosome :math:`chr`):

    .. math::
        S_{chr} = \\frac{\\sum{I_{chr}}}{||I_{chr}||}

    """
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
            header_index = dict([
                (col_name, i) for i, col_name in enumerate(row)
            ])
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
                # Chromosome 23
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

    data = np.array(
        data,
        dtype=[("sampleID", "a{}".format(max([len(i[0]) for i in data]))),
               ("chr23", float),
               ("chr24", float),
               ("gender", "a{}".format(max([len(i[3]) for i in data]))),
               ("status", "a{}".format(max([len(i[4]) for i in data])))],
    )
    return data


def checkArgs(args):
    """Checks the arguments and options.

    :param args: an object containing the options of the program.

    :type args: argparse.Namespace

    :returns: ``True`` if everything was OK.

    If there is a problem with an option, an exception is raised using the
    :py:class:`ProgramError` class, a message is printed to the
    :class:`sys.stderr` and the program exists with code 1.

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


def parseArgs(argString=None):  # pragma: no cover
    """Parses the command line options and arguments.

    :param argString: the options.

    :type argString: list

    :returns: A :py:class:`argparse.Namespace` object created by the
              :py:mod:`argparse` module. It contains the values of the
              different options.

    ============================ ====== =======================================
               Options            Type                    Description
    ============================ ====== =======================================
    ``--bfile``                  string The plink binary file containing
                                        information about markers and
                                        individuals.
    ``--intensities``            string A file containing alleles intensities
                                        for each of the markers located on the
                                        X and Y chromosome.
    ``--summarized-intensities`` string The prefix of six files containing
                                        summarized chr23 and chr24 intensities.
    ``--sex-problems``           string The file containing individuals with
                                        sex problems.
    ``--format``                 string The output file format (png, ps, pdf,
                                        or X11).
    ``--xlabel``                 string The label of the X axis.
    ``--ylabel``                 string The label of the Y axis.
    ``--out``                    string The prefix of the output files.
    ============================ ====== =======================================

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
desc = "Plots the gender using X and Y chromosomes' intensities"
parser = argparse.ArgumentParser(description=desc)
parser.add_argument("-v", "--version", action="version",
                    version="pyGenClean version {}".format(__version__))

# The INPUT files
group = parser.add_argument_group("Input File")
group.add_argument("--bfile", type=str, metavar="FILE", required=True,
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
group.add_argument("--sex-problems", type=str, metavar="FILE", required=True,
                   help=("The file containing individuals with sex "
                         "problems. This file is not read if the option "
                         "'summarized-intensities' is used."))
# The options
group = parser.add_argument_group("Options")
group.add_argument("--format", type=str, metavar="FORMAT", default="png",
                   choices=["png", "ps", "pdf", "X11"],
                   help=("The output file format (png, ps, pdf, or X11 "
                         "formats are available). [default: %(default)s]"))
group.add_argument("--xlabel", type=str, metavar="STRING",
                   default="X intensity",
                   help="The label of the X axis. [default: %(default)s]")
group.add_argument("--ylabel", type=str, metavar="STRING",
                   default="Y intensity",
                   help="The label of the Y axis. [default: %(default)s]")
# The OUTPUT files
group = parser.add_argument_group("Output File")
group.add_argument("--out", type=str, metavar="FILE",
                   default="sexcheck",
                   help=("The prefix of the output files (which will be "
                         "a Plink binary file. [default: %(default)s]"))


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
