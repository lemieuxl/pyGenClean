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
import shutil
import datetime
import argparse
import subprocess
import ConfigParser

import pyGenClean
import pyGenClean.LaTeX as latex_template
import pyGenClean.FlagHW.flag_hw as flag_hw
import pyGenClean.LaTeX.AutoReport as AutoReport
import pyGenClean.SexCheck.sex_check as sex_check
import pyGenClean.PlateBias.plate_bias as plate_bias
import pyGenClean.FlagMAF.flag_maf_zero as flag_maf_zero
import pyGenClean.DupSNPs.duplicated_snps as duplicated_snps
import pyGenClean.Ethnicity.check_ethnicity as check_ethnicity
import pyGenClean.Misc.compare_gold_standard as compare_gold_standard
import pyGenClean.DupSamples.duplicated_samples as duplicated_samples
import pyGenClean.MarkerMissingness.snp_missingness as snp_missingness
import pyGenClean.SampleMissingness.sample_missingness as sample_missingness
import pyGenClean.NoCallHetero.clean_noCall_hetero_snps as noCall_hetero_snps
import pyGenClean.RelatedSamples.find_related_samples as find_related_samples
import pyGenClean.HeteroHap.remove_heterozygous_haploid \
                                                as remove_heterozygous_haploid

import PlinkUtils.subset_data as subset_data


# Getting the version
prog_version = pyGenClean.__version__


def main():
    """The main function.

    These are the steps performed for the data clean up:

    1. Prints the version number.
    2. Reads the configuration file (:py:func:`read_config_file`).
    3. Creates a new directory with ``data_clean_up`` as prefix and the date
       and time as suffix. In the improbable event that the directory already
       exists, asks the user the permission to overwrite it.
    4. Check the input file type (``bfile``, ``tfile`` or ``file``).
    5. Creates an intermediate directory with the section as prefix and the
       script name as suffix (inside the previous directory).
    6. Runs the required script in order (according to the configuration file
       section).

    .. note::
        The main function is not responsible to check if the required files
        exist. This should be done in the ``run`` functions.

    """
    # Getting and checking the options
    args = parse_args()
    check_args(args)

    print "Data Clean Up version {}".format(prog_version)

    # Reading the configuration file
    order, conf = read_config_file(args.conf)

#     # The directory name
#     dirname = "data_clean_up."
#     dirname += datetime.datetime.today().strftime("%Y-%m-%d_%H.%M.%S")
#     if os.path.isdir(dirname):
#         answer = "N"
#         if not args.overwrite:
#             # The directory already exists...
#             print >>sys.stderr, ("WARNING: {}: directory already "
#                                  "exists".format(dirname))
#             print >>sys.stderr, "Overwrite [Y/N]? ",
#             answer = raw_input()
#         if args.overwrite or answer.upper() == "Y":
#             # Delete everything with the directory
#             shutil.rmtree(dirname)
#         elif answer.upper() == "N":
#             print >>sys.stderr, "STOPING NOW"
#             sys.exit(0)
#         else:
#             msg = "{}: not a valid answer (Y or N)".format(answer)
#             raise ProgramError(msg)
#
#     # Creating the output directory
#     os.mkdir(dirname)
    dirname = "data_clean_up.2014-09-03_15.18.17"

    # Executing the data clean up
    current_input = None
    current_input_type = None
    if args.tfile is not None:
        current_input = args.tfile
        current_input_type = "tfile"
    elif args.bfile is not None:
        current_input = args.bfile
        current_input_type = "bfile"
    else:
        current_input = args.file
        current_input_type = "file"

    latex_summaries = []
    steps = []
    descriptions = []
    for number in order:
        # Getting the script name and its options
        script_name, options = conf[number]

        # Getting the output prefix
        output_prefix = os.path.join(dirname,
                                     "{}_{}".format(number, script_name))

        # Getting the function to use
        function_to_use = available_functions[script_name]

        # Executing the function
        print "\nRunning {} {}".format(number, script_name)
        print ("   - Using {} as prefix for input "
               "files".format(current_input))
        print "   - Results will be in [ {} ]".format(output_prefix)
        current_input, current_input_type, summary, desc = function_to_use(
            current_input,
            current_input_type,
            output_prefix,
            dirname,
            options,
        )

        # Saving what's necessary for the LaTeX report
        latex_summaries.append(summary)
        steps.append(script_name)
        descriptions.append(desc)

    # A dummy background section content
    dummy_background = (
        "Lorem ipsum dolor sit amet, consectetur adipiscing elit. "
        "Suspendisse lectus ligula, volutpat eget convallis a, porttitor "
        "vitae est. Pellentesque ornare ipsum vitae odio sodales, eu "
        "elementum urna pretium. Donec luctus non leo sed euismod. Phasellus "
        "in diam et leo fringilla adipiscing ullamcorper nec sapien. Sed "
        "condimentum metus at lacus vehicula vulputate. Nam fermentum "
        "faucibus ipsum ut gravida. In sed felis tellus. Aliquam imperdiet, "
        "augue et eleifend cursus, elit risus accumsan justo, eu aliquam "
        "quam massa id risus. Donec sagittis orci lorem, a vulputate lacus "
        "sodales ut. Proin massa massa, aliquet vitae felis et, porttitor "
        "ornare enim."
    )

    # We create the automatic report
    project_name = "Dummy Project"
    logo_path = os.path.join(os.environ["HOME"], "Pictures",
                             "statgen_logo.png")
    report_name = os.path.join(dirname, "automatic_report.tex")
    AutoReport.create_report(dirname, report_name, project_name=project_name,
                             logo_path=logo_path, steps=steps,
                             descriptions=descriptions,
                             summaries=latex_summaries,
                             background=dummy_background)


def run_duplicated_samples(in_prefix, in_type, out_prefix, options):
    """Runs step1 (duplicated samples).

    :param in_prefix: the prefix of the input files.
    :param in_type: the type of the input files.
    :param out_prefix: the output prefix.
    :param options: the options needed.

    :type in_prefix: string
    :type in_type: string
    :type out_prefix: string
    :type options: list of strings

    :returns: a tuple containing the prefix of the output files (the input
              prefix for the next script) and the type of the output files
              (``tfile``).

    This function calls the :py:mod:`DupSamples.duplicated_samples` module. The
    required file type for this module is ``tfile``, hence the need to use the
    :py:func:`check_input_files` to check if the file input file type is the
    good one, or to create it if needed.

    """
    # Creating the output directory
    os.mkdir(out_prefix)

    # We know we need tfile
    required_type = "tfile"
    check_input_files(in_prefix, in_type, required_type)

    # We need to inject the name of the input file and the name of the output
    # prefix
    options += ["--{}".format(required_type), in_prefix,
                "--out", os.path.join(out_prefix, "dup_samples")]

    # We run the script
    try:
        duplicated_samples.main(options)
    except duplicated_samples.ProgramError as e:
        msg = "duplicated_samples: {}".format(e)
        raise ProgramError(msg)

    # We know this step does produce a new data set (tfile), so we return it
    return os.path.join(out_prefix, "dup_samples.final"), "tfile"


def run_duplicated_snps(in_prefix, in_type, out_prefix, options):
    """Runs step2 (duplicated snps).

    :param in_prefix: the prefix of the input files.
    :param in_type: the type of the input files.
    :param out_prefix: the output prefix.
    :param options: the options needed.

    :type in_prefix: string
    :type in_type: string
    :type out_prefix: string
    :type options: list of strings

    :returns: a tuple containing the prefix of the output files (the input
              prefix for the next script) and the type of the output files
              (``tfile``).

    This function calls the :py:mod:`DupSNPs.duplicated_snps` module. The
    required file type for this module is ``tfile``, hence the need to use the
    :py:func:`check_input_files` to check if the file input file type is the
    good one, or to create it if needed.

    .. note::
        This function creates a ``map`` file, needed for the
        :py:mod:`DupSNPs.duplicated_snps` module.

    """
    # Creating the output directory
    os.mkdir(out_prefix)

    # We know we need a tfile
    required_type = "tfile"
    check_input_files(in_prefix, in_type, required_type)

    # This step require a map file (we now have a tfile)
    if not os.path.isfile(in_prefix + ".map"):
        outputFile = None
        try:
            outputFile = open(in_prefix + ".map", "w")
        except IOError:
            msg = "{}: can't write file".format(in_prefix + ".map")
            raise ProgramError(msg)
        try:
            with open(in_prefix + ".tped", 'r') as inputFile:
                for line in inputFile:
                    row = line.rstrip("\r\n").split("\t")
                    print >>outputFile, "\t".join(row[:4])
        except IOError:
            msg = "{}: no such file".format(in_prefix + ".tped")
            raise ProgramError(msg)
        outputFile.close()

    # We need to inject the name of the input file and the name of the output
    # prefix
    options += ["--{}".format(required_type), in_prefix,
                "--out", os.path.join(out_prefix, "dup_snps")]

    # We run the script
    try:
        duplicated_snps.main(options)
    except duplicated_snps.ProgramError as e:
        msg = "duplicated_snps: {}".format(e)
        raise ProgramError(msg)

    # We know this step does produce a new data set (tfile), so we return it
    return os.path.join(out_prefix, "dup_snps.final"), "tfile"


def run_noCall_hetero_snps(in_prefix, in_type, out_prefix, base_dir, options):
    """Runs step 3 (clean no call and hetero).

    :param in_prefix: the prefix of the input files.
    :param in_type: the type of the input files.
    :param out_prefix: the output prefix.
    :param options: the options needed.

    :type in_prefix: string
    :type in_type: string
    :type out_prefix: string
    :type options: list of strings

    :returns: a tuple containing the prefix of the output files (the input
              prefix for the next script) and the type of the output files
              (``tfile``).

    This function calls the :py:mod:`NoCallHetero.clean_noCall_hetero_snps`
    module. The required file type for this module is ``tfile``, hence the need
    to use the :py:func:`check_input_files` to check if the file input file
    type is the good one, or to create it if needed.

    """
#     # Creating the output directory
#     os.mkdir(out_prefix)

    # We know we need a tfile
    required_type = "tfile"
    check_input_files(in_prefix, in_type, required_type)

    # We need to inject the name of the input file and the name of the output
    # prefix
    script_prefix = os.path.join(out_prefix, "clean_noCall_hetero")
    options += ["--{}".format(required_type), in_prefix,
                "--out", script_prefix]

    # We run the script
    try:
        noCall_hetero_snps.main(options)
    except noCall_hetero_snps.ProgramError as e:
        msg = "noCall_hetero_snps: {}".format(e)
        raise ProgramError(msg)

    # We want to save in a file the markers and samples that were removed
    # There are two files to look at, which contains only one row, the name of
    # the markers:
    #     - prefix.allFailed
    #     - prefix.allHetero
    nb_all_failed = 0
    nb_all_hetero = 0
    o_filename = os.path.join(base_dir, "excluded_markers.txt")
    try:
        with open(o_filename, "a") as o_file:
            # The first file
            i_filename = script_prefix + ".allFailed"
            if os.path.isfile(i_filename):
                with open(i_filename, "r") as i_file:
                    for line in i_file:
                        nb_all_failed += 1
                        print >>o_file, line.rstrip("\r\n") + "\tall_failed"

            # The second file
            i_filename = os.path.join(script_prefix + ".allHetero")
            if os.path.isfile(i_filename):
                with open(i_filename, "r") as i_file:
                    for line in i_file:
                        nb_all_hetero += 1
                        print >>o_file, line.rstrip("\r\n") + "\tall_hetero"

    except IOError:
        msg = "{}: can't write to file".format(o_filename)
        raise ProgramError(msg)

    # We write a LaTeX summary
    latex_file = os.path.join(script_prefix + ".summary.tex")
    try:
        with open(latex_file, "w") as o_file:
            print >>o_file, latex_template.subsection(
                noCall_hetero_snps.pretty_name,
            )
            text = (
                "After scrutiny, {:,d} marker{} were excluded from the "
                "dataset because of a call rate of 0. Also, {:,d} marker{} "
                "were excluded from the dataset because all samples were "
                "heterozygous (excluding the mitochondrial "
                "chromosome)".format(nb_all_failed,
                                     "s" if nb_all_failed > 0 else "",
                                     nb_all_hetero,
                                     "s" if nb_all_hetero > 0 else "")
            )
            print >>o_file, latex_template.wrap_lines(text, 80)

    except IOError:
        msg = "{}: cannot write LaTeX summary".format(latex_file)
        raise ProgramError(msg)

    # We know this step does produce a new data set (tfile), so we return it
    # along with the report name
    return (os.path.join(out_prefix, "clean_noCall_hetero"), "tfile",
            latex_file, noCall_hetero_snps.desc)


def run_sample_missingness(in_prefix, in_type, out_prefix, base_dir, options):
    """Runs step4 (clean mind).

    :param in_prefix: the prefix of the input files.
    :param in_type: the type of the input files.
    :param out_prefix: the output prefix.
    :param options: the options needed.

    :type in_prefix: string
    :type in_type: string
    :type out_prefix: string
    :type options: list of strings

    :returns: a tuple containing the prefix of the output files (the input
              prefix for the next script) and the type of the output files
              (``bfile``).

    This function calls the :py:mod:`SampleMissingness.sample_missingness`
    module. The required file type for this module is either a ``bfile`` or a
    ``tfile``, hence the need to use the :py:func:`check_input_files` to check
    if the file input file type is the good one, or to create it if needed.

    """
#     # Creating the output directory
#     os.mkdir(out_prefix)

    # We are looking at what we have
    required_type = "tfile"
    if in_type == "bfile":
        required_type = "bfile"
    check_input_files(in_prefix, in_type, required_type)

    # We need to inject the name of the input file and the name of the output
    # prefix
    script_prefix = os.path.join(out_prefix, "clean_mind")
    options += ["--ifile", in_prefix,
                "--out", script_prefix]
    if required_type == "bfile":
        options.append("--is-bfile")

    # We run the script
    try:
        sample_missingness.main(options)
    except sample_missingness.ProgramError as e:
        msg = "sample_missingness: {}".format(e)
        raise ProgramError(msg)

    # We want to modify the description, so that it contains the option used
    desc = sample_missingness.desc
    mind_value = sample_missingness.parser.get_default("mind")
    if "--mind" in options:
        # Since we already run the script, we know that the mind option is a
        # valid float
        mind_value = options[options.index("--mind") + 1]
    if desc[-1] == ".":
        desc = desc[:-1] + r" (${}={}$).".format(latex_template.texttt("mind"),
                                                 mind_value)

    # We want to save in a file the samples that were removed
    # There is one file to look at, which contains only one row, the name of
    # the samples:
    #     - prefix.irem (file will exists only if samples were removed)
    nb_samples = 0
    o_filename = os.path.join(base_dir, "excluded_samples.txt")
    try:
        # Checking if the file exists
        i_filename = script_prefix + ".irem"
        if os.path.isfile(i_filename):
            # True, so sample were removed
            with open(i_filename, "r") as i_file, \
                    open(o_filename, "a") as o_file:
                for line in i_file:
                    nb_samples += 1
                    print >>o_file, line.rstrip("\r\n") + "\tmind {}".format(
                        mind_value,
                    )

    except IOError:
        msg = "{}: can't write to file".format(o_filename)
        raise ProgramError(msg)

    # We write a LaTeX summary
    latex_file = os.path.join(script_prefix + ".summary.tex")
    try:
        with open(latex_file, "w") as o_file:
            print >>o_file, latex_template.subsection(
                sample_missingness.pretty_name,
            )
            text = ("Using a {} threshold of {} ({} keeping only samples with "
                    r"a missing rate $\leq {}$), {:,d} sample{} were excluded "
                    "from the dataset.".format(latex_template.texttt("mind"),
                                               mind_value,
                                               latex_template.textit("i.e."),
                                               mind_value, nb_samples,
                                               "s" if nb_samples > 0 else ""))
            print >>o_file, latex_template.wrap_lines(text)

    except IOError:
        msg = "{}: cannot write LaTeX summary".format(latex_file)
        raise ProgramError(msg)

    # We know this step does produce a new data set (bfile), so we return it
    return os.path.join(out_prefix, "clean_mind"), "bfile", latex_file, desc


def run_snp_missingness(in_prefix, in_type, out_prefix, base_dir, options):
    """Run step5 (clean geno).

    :param in_prefix: the prefix of the input files.
    :param in_type: the type of the input files.
    :param out_prefix: the output prefix.
    :param options: the options needed.

    :type in_prefix: string
    :type in_type: string
    :type out_prefix: string
    :type options: list of strings

    :returns: a tuple containing the prefix of the output files (the input
              prefix for the next script) and the type of the output files
              (``bfile``).

    This function calls the :py:mod:`MarkerMissingness.snp_missingness` module.
    The required file type for this module is ``bfile``, hence the need to use
    the :py:func:`check_input_files` to check if the file input file type is
    the good one, or to create it if needed.

    """
#     # Creating the output directory
#     os.mkdir(out_prefix)

    # We know we need a bfile
    required_type = "bfile"
    check_input_files(in_prefix, in_type, required_type)

    # We need to inject the name of the input file and the name of the output
    # prefix
    script_prefix = os.path.join(out_prefix, "clean_geno")
    options += ["--{}".format(required_type), in_prefix,
                "--out", script_prefix]

    # We run the script
    try:
        snp_missingness.main(options)
    except snp_missingness.ProgramError as e:
        msg = "snp_missingness: {}".format(e)
        raise ProgramError(msg)

    # We want to modify the description, so that it contains the option used
    desc = snp_missingness.desc
    geno_value = snp_missingness.parser.get_default("geno")
    if "--geno" in options:
        # Since we already run the script, we know that the mind option is a
        # valid float
        geno_value = options[options.index("--geno") + 1]
    if desc[-1] == ".":
        desc = desc[:-1] + r" (${}={}$).".format(latex_template.texttt("mind"),
                                                 geno_value)

    # We want to save in a file the samples that were removed
    # There is one file to look at, which contains only one row, the name of
    # the samples:
    #     - prefix.removed_snps
    nb_markers = 0
    o_filename = os.path.join(base_dir, "excluded_markers.txt")
    try:
        # Checking if the file exists
        i_filename = script_prefix + ".removed_snps"
        if os.path.isfile(i_filename):
            # True, so sample were removed
            with open(i_filename, "r") as i_file, \
                    open(o_filename, "a") as o_file:
                for line in i_file:
                    nb_markers += 1
                    print >>o_file, line.rstrip("\r\n") + "\tgeno {}".format(
                        geno_value,
                    )

    except IOError:
        msg = "{}: can't write to file".format(o_filename)
        raise ProgramError(msg)

    # We write a LaTeX summary
    latex_file = os.path.join(script_prefix + ".summary.tex")
    try:
        with open(latex_file, "w") as o_file:
            print >>o_file, latex_template.subsection(
                snp_missingness.pretty_name,
            )
            text = ("Using a {} threshold of {} ({} keeping only markers with "
                    r"a missing rate $\leq {}$), {:,d} marker{} were excluded "
                    "from the dataset.".format(latex_template.texttt("geno"),
                                               geno_value,
                                               latex_template.textit("i.e."),
                                               geno_value, nb_markers,
                                               "s" if nb_markers > 0 else ""))
            print >>o_file, latex_template.wrap_lines(text)

    except IOError:
        msg = "{}: cannot write LaTeX summary".format(latex_file)
        raise ProgramError(msg)

    # We know this step does produce a new data set (bfile), so we return it
    return os.path.join(out_prefix, "clean_geno"), "bfile", latex_file, desc


def run_sex_check(in_prefix, in_type, out_prefix, base_dir, options):
    """Runs step6 (sexcheck).

    :param in_prefix: the prefix of the input files.
    :param in_type: the type of the input files.
    :param out_prefix: the output prefix.
    :param options: the options needed.

    :type in_prefix: string
    :type in_type: string
    :type out_prefix: string
    :type options: list of strings

    :returns: a tuple containing the prefix of the output files (the input
              prefix for the next script) and the type of the output files
              (``bfile``).

    This function calls the :py:mod:`SexCheck.sex_check` module. The required
    file type for this module is ``bfile``, hence the need to use the
    :py:func:`check_input_files` to check if the file input file type is the
    good one, or to create it if needed.

    .. note::
        The :py:mod:`SexCheck.sex_check` module doesn't return usable output
        files. Hence, this function returns the input file prefix and its type.

    """
#     # Creating the output directory
#     os.mkdir(out_prefix)

    # We know we need a bfile
    required_type = "bfile"
    check_input_files(in_prefix, in_type, required_type)

    # We need to inject the name of the input file and the name of the output
    # prefix
    script_prefix = os.path.join(out_prefix, "sexcheck")
    options += ["--{}".format(required_type), in_prefix,
                "--out", script_prefix]

    # We run the script
    try:
        sex_check.main(options)
    except sex_check.ProgramError as e:
        msg = "sex_check {}".format(e)
        raise ProgramError(msg)

    # We write a LaTeX summary
    latex_file = os.path.join(script_prefix + ".summary.tex")
    try:
        with open(latex_file, "w") as o_file:
            print >>o_file, latex_template.subsection(sex_check.pretty_name)
            text = (
                "Lorem ipsum dolor sit amet, consectetur adipiscing elit. "
                "Aenean imperdiet libero id laoreet vulputate. Sed pretium "
                "malesuada sapien nec blandit. Maecenas iaculis metus "
                "ultricies volutpat varius. Etiam vulputate nisi augue, a "
                "dapibus turpis convallis sit amet. Fusce tempor dolor sed "
                "nulla varius malesuada. Phasellus euismod lectus sed "
                "velit auctor, quis viverra nulla consequat. Donec "
                "tincidunt viverra nisi ut efficitur. Sed ultrices nisl "
                "diam, quis efficitur neque maximus et. Vestibulum commodo "
                "mi sit amet euismod congue."
            )
            print >>o_file, latex_template.wrap_lines(text)

    except IOError:
        msg = "{}: cannot write LaTeX summary".format(latex_file)
        raise ProgramError(msg)

    # We know this step does not produce a new data set, so we return the
    # original one
    return in_prefix, required_type, latex_file, sex_check.desc


def run_plate_bias(in_prefix, in_type, out_prefix, options):
    """Runs step7 (plate bias).

    :param in_prefix: the prefix of the input files.
    :param in_type: the type of the input files.
    :param out_prefix: the output prefix.
    :param options: the options needed.

    :type in_prefix: string
    :type in_type: string
    :type out_prefix: string
    :type options: list of strings

    :returns: a tuple containing the prefix of the output files (the input
              prefix for the next script) and the type of the output files
              (``bfile``).

    This function calls the :py:mod:`PlateBias.plate_bias` module. The required
    file type for this module is ``bfile``, hence the need to use the
    :py:func:`check_input_files` to check if the file input file type is the
    good one, or to create it if needed.

    .. note::
        The :py:mod:`PlateBias.plate_bias` module doesn't return usable output
        files. Hence, this function returns the input file prefix and its type.

    """
    # Creating the output directory
    os.mkdir(out_prefix)

    # We know we need bfile
    required_type = "bfile"
    check_input_files(in_prefix, in_type, required_type)

    # We need to inject the name of the input file and the name of the output
    # prefix
    options += ["--{}".format(required_type), in_prefix,
                "--out", os.path.join(out_prefix, "plate_bias")]

    # We run the script
    try:
        plate_bias.main(options)
    except plate_bias.ProgramError as e:
        msg = "plate_bias: {}".format(e)
        raise ProgramError(msg)

    # We know this step doesn't produce an new data set, so we return the old
    # prefix and the old in_type
    return in_prefix, required_type


def run_remove_heterozygous_haploid(in_prefix, in_type, out_prefix, options):
    """Runs step8 (remove heterozygous haploid).

    :param in_prefix: the prefix of the input files.
    :param in_type: the type of the input files.
    :param out_prefix: the output prefix.
    :param options: the options needed.

    :type in_prefix: string
    :type in_type: string
    :type out_prefix: string
    :type options: list of strings

    :returns: a tuple containing the prefix of the output files (the input
              prefix for the next script) and the type of the output files
              (``bfile``).

    This function calls the :py:mod:`HeteroHap.remove_heterozygous_haploid`
    module. The required file type for this module is ``bfile``, hence the need
    to use the :py:func:`check_input_files` to check if the file input file
    type is the good one, or to create it if needed.

    """
    # Creating the output directory
    os.mkdir(out_prefix)

    # We know we need bfile
    required_type = "bfile"
    check_input_files(in_prefix, in_type, required_type)

    # We need to inject the name of the input file and the name of the output
    # prefix
    options += ["--{}".format(required_type), in_prefix,
                "--out", os.path.join(out_prefix, "without_hh_genotypes")]

    # We run the script
    try:
        remove_heterozygous_haploid.main(options)
    except remove_heterozygous_haploid.ProgramError as e:
        msg = "remove_heterozygous_haploid: {}".format(e)
        raise ProgramError(msg)

    # We know this step produces an new data set (bfile), so we return it
    return os.path.join(out_prefix, "without_hh_genotypes"), "bfile"


def run_find_related_samples(in_prefix, in_type, out_prefix, options):
    """Runs step9 (find related samples).

    :param in_prefix: the prefix of the input files.
    :param in_type: the type of the input files.
    :param out_prefix: the output prefix.
    :param options: the options needed.

    :type in_prefix: string
    :type in_type: string
    :type out_prefix: string
    :type options: list of strings

    :returns: a tuple containing the prefix of the output files (the input
              prefix for the next script) and the type of the output files
              (``bfile``).

    This function calls the :py:mod:`RelatedSamples.find_related_samples`
    module. The required file type for this module is ``bfile``, hence the need
    to use the :py:func:`check_input_files` to check if the file input file
    type is the good one, or to create it if needed.

    .. note::
        The :py:mod:`RelatedSamples.find_related_samples` module doesn't return
        usable output files. Hence, this function returns the input file prefix
        and its type.

    """
    # Creating the output directory
    os.mkdir(out_prefix)

    # We know we need bfile
    required_type = "bfile"
    check_input_files(in_prefix, in_type, required_type)

    # We need to inject the name of the input file and the name of the output
    # prefix
    options += ["--{}".format(required_type), in_prefix,
                "--out", os.path.join(out_prefix, "ibs")]

    # We run the script
    try:
        find_related_samples.main(options)
    except find_related_samples.ProgramError as e:
        msg = "find_related_samples: {}".format(e)
        raise ProgramError(msg)

    # We know this step doesn't produce an new data set, so we return the old
    # prefix and the old in_type
    return in_prefix, required_type


def run_check_ethnicity(in_prefix, in_type, out_prefix, options):
    """Runs step10 (check ethnicity).

    :param in_prefix: the prefix of the input files.
    :param in_type: the type of the input files.
    :param out_prefix: the output prefix.
    :param options: the options needed.

    :type in_prefix: string
    :type in_type: string
    :type out_prefix: string
    :type options: list of strings

    :returns: a tuple containing the prefix of the output files (the input
              prefix for the next script) and the type of the output files
              (``bfile``).

    This function calls the :py:mod:`Ethnicity.check_ethnicity` module. The
    required file type for this module is ``bfile``, hence the need to use the
    :py:func:`check_input_files` to check if the file input file type is the
    good one, or to create it if needed.

    .. note::
        The :py:mod:`Ethnicity.check_ethnicity` module doesn't return usable
        output files. Hence, this function returns the input file prefix and
        its type.

    """
    # Creating the output directory
    os.mkdir(out_prefix)

    # We know we need bfile
    required_type = "bfile"
    check_input_files(in_prefix, in_type, required_type)

    # We need to inject the name of the input file and the name of the output
    # prefix
    options += ["--{}".format(required_type), in_prefix,
                "--out", os.path.join(out_prefix, "ethnicity")]

    # We run the script
    try:
        check_ethnicity.main(options)
    except check_ethnicity.ProgramError as e:
        msg = "check_ethnicity: {}".format(e)
        raise ProgramError(msg)

    # We know this step doesn't produce an new data set, so we return the old
    # prefix and the old in_type
    return in_prefix, required_type


def run_flag_maf_zero(in_prefix, in_type, out_prefix, options):
    """Runs step11 (flag MAF zero).

    :param in_prefix: the prefix of the input files.
    :param in_type: the type of the input files.
    :param out_prefix: the output prefix.
    :param options: the options needed.

    :type in_prefix: string
    :type in_type: string
    :type out_prefix: string
    :type options: list of strings

    :returns: a tuple containing the prefix of the output files (the input
              prefix for the next script) and the type of the output files
              (``bfile``).

    This function calls the :py:mod:`FlagMAF.flag_maf_zero` module. The
    required file type for this module is ``bfile``, hence the need to use the
    :py:func:`check_input_files` to check if the file input file type is the
    good one, or to create it if needed.

    .. note::
        The :py:mod:`FlagMAF.flag_maf_zero` module doesn't return usable output
        files. Hence, this function returns the input file prefix and its type.

    """
    # Creating the output directory
    os.mkdir(out_prefix)

    # We know we need bfile
    required_type = "bfile"
    check_input_files(in_prefix, in_type, required_type)

    # We need to inject the name of the input file and the name of the output
    # prefix
    options += ["--{}".format(required_type), in_prefix,
                "--out", os.path.join(out_prefix, "flag_maf_0")]

    # We run the script
    try:
        flag_maf_zero.main(options)
    except flag_maf_zero.ProgramError as e:
        msg = "flag_maf_zero: {}".format(e)
        raise ProgramError(msg)

    # We know this step doesn't produce an new data set, so we return the old
    # prefix and the old in_type
    return in_prefix, required_type


def run_flag_hw(in_prefix, in_type, out_prefix, options):
    """Runs step12 (flag HW).

    :param in_prefix: the prefix of the input files.
    :param in_type: the type of the input files.
    :param out_prefix: the output prefix.
    :param options: the options needed.

    :type in_prefix: string
    :type in_type: string
    :type out_prefix: string
    :type options: list of strings

    :returns: a tuple containing the prefix of the output files (the input
              prefix for the next script) and the type of the output files
              (``bfile``).

    This function calls the :py:mod:`FlagHW.flag_hw` module. The required file
    type for this module is ``bfile``, hence the need to use the
    :py:func:`check_input_files` to check if the file input file type is the
    good one, or to create it if needed.

    .. note::
        The :py:mod:`FlagHW.flag_hw` module doesn't return usable output files.
        Hence, this function returns the input file prefix and its type.

    """
    # Creating the output directory
    os.mkdir(out_prefix)

    # We know we need bfile
    required_type = "bfile"
    check_input_files(in_prefix, in_type, required_type)

    # We need to inject the name of the input file and the name of the output
    # prefix
    options += ["--{}".format(required_type), in_prefix,
                "--out", os.path.join(out_prefix, "flag_hw")]

    # We run the script
    try:
        flag_hw.main(options)
    except flag_hw.ProgramError as e:
        msg = "flag_hw: {}".format(e)
        raise ProgramError(msg)

    # We know this step doesn't produce an new data set, so we return the old
    # prefix and the old in_type
    return in_prefix, required_type


def run_compare_gold_standard(in_prefix, in_type, out_prefix, options):
    """Compares with a gold standard data set (compare_gold_standard.

    :param in_prefix: the prefix of the input files.
    :param in_type: the type of the input files.
    :param out_prefix: the output prefix.
    :param options: the options needed.

    :type in_prefix: string
    :type in_type: string
    :type out_prefix: string
    :type options: list of strings

    :returns: a tuple containing the prefix of the output files (the input
              prefix for the next script) and the type of the output files
              (``bfile``).

    This function calls the :py:mod:`Misc.compare_gold_standard` module. The
    required file type for this module is ``bfile``, hence the need to use the
    :py:func:`check_input_files` to check if the file input file type is the
    good one, or to create it if needed.

    .. note::
        The :py:mod:`Misc.compare_gold_standard` module doesn't return usable
        output files. Hence, this function returns the input file prefix and
        its type.

    """
    # Creating the output directory
    os.mkdir(out_prefix)

    # We know we need bfile
    required_type = "bfile"
    check_input_files(in_prefix, in_type, required_type)

    # We need to inject the name of the input file and the name of the output
    # prefix
    options += ["--{}".format(required_type), in_prefix,
                "--out", os.path.join(out_prefix, "compare_with_gold")]

    # We run the script
    try:
        compare_gold_standard.main(options)
    except compare_gold_standard.ProgramError as e:
        msg = "compare_gold_standard: {}".format(e)
        raise ProgramError(msg)

    # We know this step doesn't produce an new data set, so we return the old
    # prefix and the old in_type
    return in_prefix, required_type


def run_subset_data(in_prefix, in_type, out_prefix, options):
    """Subsets the data.

    :param in_prefix: the prefix of the input files.
    :param in_type: the type of the input files.
    :param out_prefix: the output prefix.
    :param options: the options needed.

    :type in_prefix: string
    :type in_type: string
    :type out_prefix: string
    :type options: list of strings

    :returns: a tuple containing the prefix of the output files (the input
              prefix for the next script) and the type of the output files
              (``bfile``).

    This function calls the :py:mod:`PlinkUtils.subset_data` module. The
    required file type for this module is ``bfile``, hence the need to use the
    :py:func:`check_input_files` to check if the file input file type is the
    good one, or to create it if needed.

    .. note::
        The output file type is the same as the input file type.

    """
    # Creating the output directory
    os.mkdir(out_prefix)

    # Looking at what we have
    required_type = None
    if in_type == "bfile":
        required_type = "bfile"
    elif in_type == "tfile":
        required_type = "tfile"
    else:
        required_type = "file"
    if "--is-bfile" in set(options):
        required_type = "bfile"
    check_input_files(in_prefix, in_type, required_type)

    # We need to inject the name of the input file and the name of the output
    # prefix
    options += ["--ifile", in_prefix,
                "--out", os.path.join(out_prefix, "subset")]
    if required_type == "bfile":
        options.append("--is-bfile")
    elif required_type == "tfile":
        options.append("--is-tfile")
    else:
        options.append("--is-file")

    # We run the script
    try:
        subset_data.main(options)
    except subset_data.ProgramError as e:
        msg = "subset_data: {}".format(e)
        raise ProgramError(msg)

    # We know this step does produce a new data set (bfile), so we return it
    return os.path.join(out_prefix, "subset"), required_type


def run_command(command):
    """Run a command using subprocesses.

    :param command: the command to run.

    :type command: list of strings

    Tries to run a command. If it fails, raise a :py:class:`ProgramError`.

    .. warning::
        The variable ``command`` should be a list of strings (no other type).

    """
    output = None
    try:
        output = subprocess.check_output(command, stderr=subprocess.STDOUT,
                                         shell=False)
    except subprocess.CalledProcessError:
        msg = "couldn't run command\n{}".format(command)
        raise ProgramError(msg)


def check_input_files(prefix, the_type, required_type):
    """Check that the file is of a certain file type.

    :param prefix: the prefix of the input files.
    :param the_type: the type of the input files (bfile, tfile or file).
    :param required_type: the required type of the input files (bfile, tfile or
                          file).

    :type prefix: string
    :type the_type: string
    :type required_type: string

    :returns: ``True`` if everything is OK.

    Checks if the files are of the required type, according to their current
    type. The available types are ``bfile`` (binary), ``tfile`` (transposed)
    and ``file`` (normal).

    """
    # The files required for each type
    bfile_type = {".bed", ".bim", ".fam"}
    tfile_type = {".tped", ".tfam"}
    file_type = {".ped", ".map"}

    # Check if of bfile, tfile and file
    plink_command = ["plink", "--noweb", "--out", prefix]
    if required_type == "bfile":
        # We need bfile
        plink_command += ["--make-bed"]
        if the_type == "bfile":
            return True
        elif the_type == "tfile":
            # We have tfile, we need to create bfile from tfile
            plink_command += ["--tfile", prefix]
        elif the_type == "file":
            # We have file, we need to create bfile from file
            plink_command += ["--file", prefix]
        else:
            msg = "{}: no suitable input format...".format(prefix)
            raise ProgramError(msg)

        # We create the required files
        if os.path.isfile(prefix + ".log"):
            # There is a log file... we need to copy it
            shutil.copyfile(prefix + ".log", prefix + ".olog")
        run_command(plink_command)

        # Everything is now fine
        return True

    elif required_type == "tfile":
        # We need a tfile
        plink_command += ["--recode", "--transpose", "--tab"]
        if the_type == "tfile":
            return True
        elif the_type == "bfile":
            # We have bfile, we need to create tfile from bfile
            plink_command += ["--bfile", prefix]
        elif the_type == "file":
            # We have file, we need to create tfile from file
            plink_command += ["--file", prefix]
        else:
            msg = "{}: no suitable input format...".format(prefix)
            raise ProgramError(msg)

        # We create the required files
        if os.path.isfile(prefix + ".log"):
            # There is a log file... we need to copy it
            shutil.copyfile(prefix + ".log", prefix + ".olog")
        run_command(plink_command)

        # Everything is now fine
        return True

    elif required_type == "file":
        # We need a file
        plink_command += ["--recode", "--tab"]
        if the_type == "file":
            return True
        elif the_type == "bfile":
            # We have bfile, we need to create file from bfile
            plink_command += ["--bfile", prefix]
        elif the_type == "tfile":
            # We have tfile, we need to create file from tfile
            plink_command += ["--tfile", prefix]
        else:
            msg = "{}: no suitable input format...".format(prefix)
            raise ProgramError(msg)

        # We create the required files
        if os.path.isfile(prefix + ".log"):
            # There is a log file... we need to copy it
            shutil.copyfile(prefix + ".log", prefix + ".olog")
        run_command(plink_command)

        # Everything is now fine
        return True

    else:
        msg = "{}: unknown file format".format(required_type)
        raise ProgramError(msg)


def all_files_exist(file_list):
    """Check if all files exist.

    :param file_list: the names of files to check.

    :type file_list: list of strings

    :returns: ``True`` if all files exist, ``False`` otherwise.

    """
    all_exist = True
    for filename in file_list:
        all_exist = all_exist and os.path.isfile(filename)
    return all_exist


def read_config_file(filename):
    """Reads the configuration file.

    :param filename: the name of the file containing the configuration.

    :type filename: string

    :returns: A tuple where the first element is a list of sections, and the
              second element is a map containing the configuration (options and
              values).

    The structure of the configuration file is important. Here is an example of
    a configuration file::

        [1] # Computes statistics on duplicated samples
        script = duplicated_samples

        [2] # Removes samples according to missingness
        script = sample_missingness

        [3] # Removes markers according to missingness
        script = snp_missingness

        [4] # Removes samples according to missingness (98%)
        script = sample_missingness
        mind = 0.02

        [5] # Performs a sex check
        script = sex_check

        [6] # Flags markers with MAF=0
        script = flag_maf_zero

        [7] # Flags markers according to Hardy Weinberg
        script = flag_hw

        [8] # Subset the dataset (excludes markers and remove samples)
        script = subset
        exclude = .../filename
        rempove = .../filename

    Sections are in square brackets and must be ``integer``. The section number
    represent the step at which the script will be run (*i.e.* from the
    smallest number to the biggest). The sections must be continuous.

    Each section contains the script names (``script`` variable) and options of
    the script (all other variables) (*e.g.* section 4 runs the
    ``sample_missingness`` script (:py:func:`run_sample_missingness`) with
    option ``mind`` sets to 0.02).

    Here is a list of the available scripts:

    * ``duplicated_samples`` (:py:func:`run_duplicated_samples`)
    * ``duplicated_snps`` (:py:func:`run_duplicated_snps`)
    * ``noCall_hetero_snps`` (:py:func:`run_noCall_hetero_snps`)
    * ``sample_missingness`` (:py:func:`run_sample_missingness`)
    * ``snp_missingness`` (:py:func:`run_snp_missingness`)
    * ``sex_check`` (:py:func:`run_sex_check`)
    * ``plate_bias`` (:py:func:`run_plate_bias`)
    * ``remove_heterozygous_haploid``
      (:py:func:`run_remove_heterozygous_haploid`)
    * ``find_related_samples`` (:py:func:`run_find_related_samples`)
    * ``check_ethnicity`` (:py:func:`run_check_ethnicity`)
    * ``flag_maf_zero`` (:py:func:`run_flag_maf_zero`)
    * ``flag_hw`` (:py:func:`run_flag_hw`)
    * ``subset`` (:py:func:`run_subset_data`)
    * ``compare_gold_standard`` (:py:func:`run_compare_gold_standard`)

    """
    # Creating the config parser
    config = ConfigParser.RawConfigParser(allow_no_value=True)
    config.optionxform = str
    config.read(filename)

    # Checking the section names
    sections = None
    try:
        sections = sorted([int(i) for i in config.sections()])
    except ValueError:
        # Section not integer
        msg = ("{}: sections must be integers: "
               "{}".format(filename, config.sections()))
        raise ProgramError(msg)
    if sections != range(min(sections), max(sections)+1):
        # Missing a section
        msg = "{}: maybe a section is missing: {}".format(filename, sections)
        raise ProgramError(msg)
    sections = [str(i) for i in sections]

    # Reading the configuration for each sections
    configuration = {}
    for section in sections:
        # Getting the script variable (and check it)
        script_name = None
        try:
            script_name = config.get(section, "script")
        except ConfigParser.NoOptionError:
            msg = ("{}: section {}: no variable called "
                   "'script'".format(filename, section))
            raise ProgramError(msg)
        if script_name not in available_modules:
            msg = ("{}: section {}: script {}: invalid script "
                   "name".format(filename, section, script_name))
            raise ProgramError(msg)

        # Getting the variables
        options = []
        for variable_name, variable_value in config.items(section):
            unwanted_options = {"bfile", "tfile", "file", "out"}
            if script_name in {"sample_missingness", "subset"}:
                unwanted_options |= {"ifile", "is-bfile", "is-tfile"}
            if script_name == "subset":
                unwanted_options.add("is-file")
            for unwanted in unwanted_options:
                if variable_name == unwanted:
                    msg = ("{}: section {}: do not use {} as an option for "
                           "{}".format(filename, section, unwanted,
                                       script_name))
                    raise ProgramError(msg)
            if variable_name != "script":
                options.append("--" + variable_name)
                if variable_value is not None:
                    if variable_name in {"indep-pairwise", "sge-nodes",
                                         "ibs-sge-nodes"}:
                        # This is a special option
                        options.extend(variable_value.split(" "))
                    else:
                        options.append(variable_value)

        # Saving the configuration
        configuration[section] = (script_name, options)

    return sections, configuration


def check_args(args):
    """Checks the arguments and options.

    :param args: an object containing the options and arguments of the program.

    :type args: :py:class:`argparse.Namespace`

    :returns: ``True`` if everything was OK.

    If there is a problem with an option, an exception is raised using the
    :py:class:`ProgramError` class, a message is printed to the
    :class:`sys.stderr` and the program exits with error code 1.

    """
    # Checking the configuration file
    if not os.path.isfile(args.conf):
        msg = "{}: no such file".format(args.conf)
        raise ProgramError(msg)

    # Check the input files
    if args.bfile is None and args.tfile is None and args.file is None:
        msg = "needs one input file prefix (--bfile, --tfile or --file)"
        raise ProgramError(msg)
    if args.bfile is not None and args.tfile is None and args.file is None:
        for fileName in [args.bfile + i for i in [".bed", ".bim", ".fam"]]:
            if not os.path.isfile(fileName):
                msg = "{}: no such file".format(fileName)
                raise ProgramError(msg)
    elif args.tfile is not None and args.bfile is None and args.file is None:
        for fileName in [args.tfile + i for i in [".tped", ".tfam"]]:
            if not os.path.isfile(fileName):
                msg = "{}: no such file".format(fileName)
                raise ProgramError(msg)
    elif args.file is not None and args.bfile is None and args.tfile is None:
        for fileName in [args.file + i for i in [".ped", ".map"]]:
            if not os.path.isfile(fileName):
                msg = "{}: no such file". format(fileName)
                raise ProgramError(msg)
    else:
        msg = "needs only one input file prefix (--bfile, --tfile or --file)"
        raise ProgramError(msg)

    return True


def parse_args():
    """Parses the command line options and arguments.

    :returns: A :py:class:`argparse.Namespace` object created by the
              :py:mod:`argparse` module. It contains the values of the
              different options.

    ===============   =======  ================================================
        Options        Type                      Description
    ===============   =======  ================================================
    ``--bfile``       String   The input binary file prefix from Plink.
    ``--tfile``       String   The input transposed file prefix from Plink.
    ``--file``        String   The input file prefix from Plink.
    ``--conf``        String   The parameter file for the data clean up.
    ``--overwrite``   Boolean  Overwrites output directories without asking the
                               user.
    ===============   =======  ================================================

    .. note::
        No option check is done here (except for the one automatically done by
        :py:mod:`argparse`). Those need to be done elsewhere (see
        :py:func:`checkArgs`).

    """
    return parser.parse_args()


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
desc = """Runs the data clean up (version {}).""".format(prog_version)
parser = argparse.ArgumentParser(description=desc)
parser.add_argument("-v", "--version", action="version",
                    version="%(prog)s {}".format(prog_version))

group = parser.add_argument_group("Input File")
group.add_argument("--bfile", type=str, metavar="FILE",
                   help=("The input file prefix (will find the plink "
                         "binary files by appending the prefix to the "
                         ".bim, .bed and .fam files, respectively)."))
group.add_argument("--tfile", type=str, metavar="FILE",
                   help=("The input file prefix (will find the plink "
                         "transposed files by appending the prefix to the "
                         ".tped and .tfam files, respectively)."))
group.add_argument("--file", type=str, metavar="FILE",
                   help=("The input file prefix (will find the plink "
                         "files by appending the prefix to the "
                         ".ped and .fam files)."))

group = parser.add_argument_group("Configuration File")
group.add_argument("--conf", type=str, metavar="FILE", required=True,
                   help="The parameter file for the data clean up.")

group = parser.add_argument_group("General Options")
group.add_argument("--overwrite", action="store_true",
                   help=("Overwrites output directories without asking the "
                         "user. [DANGEROUS]"))

# The available modules
available_modules = {
    "duplicated_samples": duplicated_samples,
    "duplicated_snps": duplicated_snps,
    "noCall_hetero_snps": noCall_hetero_snps,
    "sample_missingness": sample_missingness,
    "snp_missingness": snp_missingness,
    "sex_check": sex_check,
    "plate_bias": plate_bias,
    "remove_heterozygous_haploid": remove_heterozygous_haploid,
    "find_related_samples": find_related_samples,
    "check_ethnicity": check_ethnicity,
    "flag_maf_zero": flag_maf_zero,
    "flag_hw": flag_hw,
    "subset": subset_data,
    "compare_gold_standard": compare_gold_standard,
}
available_functions = {
    "duplicated_samples": run_duplicated_samples,
    "duplicated_snps": run_duplicated_snps,
    "noCall_hetero_snps": run_noCall_hetero_snps,
    "sample_missingness": run_sample_missingness,
    "snp_missingness": run_snp_missingness,
    "sex_check": run_sex_check,
    "plate_bias": run_plate_bias,
    "remove_heterozygous_haploid": run_remove_heterozygous_haploid,
    "find_related_samples": run_find_related_samples,
    "check_ethnicity": run_check_ethnicity,
    "flag_maf_zero": run_flag_maf_zero,
    "flag_hw": run_flag_hw,
    "subset": run_subset_data,
    "compare_gold_standard": run_compare_gold_standard,
}

# Calling the main, if necessary
if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print >>sys.stderr, "Cancelled by user"
        sys.exit(0)
    except ProgramError as e:
        parser.error(e.message)
