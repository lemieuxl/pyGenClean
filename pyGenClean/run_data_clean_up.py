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
import re
import sys
import time
import shutil
import logging
import datetime
import argparse
import subprocess
import ConfigParser
from glob import glob
from collections import namedtuple, Counter, defaultdict

from . import __version__, add_file_handler_to_root

from .pipeline_error import ProgramError

from .FlagHW import flag_hw
from .SexCheck import sex_check
from .PlateBias import plate_bias
from .FlagMAF import flag_maf_zero
from .DupSNPs import duplicated_snps
from .Ethnicity import check_ethnicity
from .Misc import compare_gold_standard
from .Contamination import contamination
from .DupSamples import duplicated_samples
from .MarkerMissingness import snp_missingness
from .RelatedSamples import find_related_samples
from .SampleMissingness import sample_missingness
from .HeteroHap import remove_heterozygous_haploid
from .NoCallHetero import clean_noCall_hetero_snps as noCall_hetero_snps

from .LaTeX import auto_report
from .LaTeX import utils as latex_template
from .LaTeX.merge_reports import add_custom_options as report_options

from .PlinkUtils import subset_data
from .PlinkUtils import get_plink_version
from .PlinkUtils import createRowFromPlinkSpacedOutput


# Configuring logging
logger = logging.getLogger("pyGenClean")


# A namedtuple that will represent what each steps returns
_StepResult = namedtuple(
    "_StepResult",
    ["next_file", "next_file_type", "latex_summary", "description",
     "long_description", "graph_path"],
)


def main():
    """The main function.

    These are the steps performed for the data clean up:

    1. Prints the version number.
    2. Reads the configuration file (:py:func:`read_config_file`).
    3. Creates a new directory with ``data_clean_up`` as prefix and the date
       and time as suffix.
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

    # The directory name
    dirname = "data_clean_up."
    dirname += datetime.datetime.today().strftime("%Y-%m-%d_%H.%M.%S")
    while os.path.isdir(dirname):
        time.sleep(1)
        dirname = "data_clean_up."
        dirname += datetime.datetime.today().strftime("%Y-%m-%d_%H.%M.%S")

    # Creating the output directory
    os.mkdir(dirname)

    # Configuring the root logger
    add_file_handler_to_root(os.path.join(dirname, "pyGenClean.log"))

    logger.info("pyGenClean version {}".format(__version__))
    plink_version = get_plink_version()
    logger.info("Using Plink version {}".format(plink_version))

    # Reading the configuration file
    logger.info("Reading configuration file [ {} ]".format(args.conf))
    order, conf = read_config_file(args.conf)

    # Executing the data clean up
    current_in = None
    current_in_type = None
    suffixes = None
    if args.tfile is not None:
        current_in = args.tfile
        current_in_type = "tfile"
        suffixes = (".tped", ".tfam")
    elif args.bfile is not None:
        current_in = args.bfile
        current_in_type = "bfile"
        suffixes = (".bed", ".bim", ".fam")
    else:
        current_in = args.file
        current_in_type = "file"
        suffixes = (".ped", ".map")

    # Creating the excluded files
    try:
        with open(os.path.join(dirname, "excluded_markers.txt"), "w") as o_f:
            pass
        with open(os.path.join(dirname, "excluded_samples.txt"), "w") as o_f:
            pass
        with open(os.path.join(dirname, "initial_files.txt"), "w") as o_file:
            for s in suffixes:
                print >>o_file, current_in + s

    except IOError:
        msg = "{}: cannot write summary".format(dirname)
        raise ProgramError(msg)

    # Counting the number of markers and samples in the datafile
    logger.info("Counting initial number of samples and markers")
    nb_markers, nb_samples = count_markers_samples(current_in,
                                                   current_in_type)
    logger.info("  - {:,d} samples".format(nb_samples))
    logger.info("  - {:,d} markers".format(nb_markers))

    # Creating the result summary file containing the initial numbers
    try:
        with open(os.path.join(dirname, "results_summary.txt"), "w") as o_file:
            print >>o_file, "# initial"
            print >>o_file, ("Initial number of markers\t"
                             "{:,d}".format(nb_markers))
            print >>o_file, ("Initial number of samples\t"
                             "{:,d}".format(nb_samples))
            print >>o_file, "---"
    except IOError:
        msg = "{}: cannot write summary".format(dirname)
        raise ProgramError(msg)

    latex_summaries = []
    steps = []
    descriptions = []
    long_descriptions = []
    graphic_paths = set()
    for number in order:
        # Getting the script name and its options
        script_name, options = conf[number]

        # Getting the output prefix
        output_prefix = os.path.join(dirname,
                                     "{}_{}".format(number, script_name))

        # Getting the function to use
        function_to_use = available_functions[script_name]

        # Executing the function
        logger.info("Running {} {}".format(number, script_name))
        logger.info("  - Using {} as prefix for input "
                    "files".format(current_in))
        logger.info("  - Results will be in [ {} ]".format(output_prefix))

        # Executing the function
        step_results = function_to_use(
            in_prefix=current_in,
            in_type=current_in_type,
            out_prefix=output_prefix,
            base_dir=dirname,
            options=options,
        )

        # Updating the input files and input file types
        current_in = step_results.next_file
        current_in_type = step_results.next_file_type

        # Saving what's necessary for the LaTeX report
        latex_summaries.append(step_results.latex_summary)
        steps.append(script_name)
        descriptions.append(step_results.description)
        long_descriptions.append(step_results.long_description)
        if step_results.graph_path is not None:
            graphic_paths.update(step_results.graph_path)

    # Counting the final number of samples and markers
    logger.info("Counting final number of samples and markers")
    nb_markers, nb_samples = count_markers_samples(current_in,
                                                   current_in_type)
    logger.info("  - {:,d} samples".format(nb_samples))
    logger.info("  - {:,d} markers".format(nb_markers))

    # Getting the final suffixes
    suffixes = None
    if current_in_type == "tfile":
        suffixes = ((".tped", nb_markers), (".tfam", nb_samples))
    elif current_in_type == "bfile":
        suffixes = ((".bed", None), (".bim", nb_markers), (".fam", nb_samples))
    else:
        suffixes = ((".ped", nb_samples), (".map", nb_markers))

    with open(os.path.join(dirname, "final_files.txt"), "w") as o_file:
        for s, nb in suffixes:
            if nb:
                print >>o_file, current_in + s + "\t{:,d}".format(nb)
            else:
                print >>o_file, current_in + s

    # Generating the graphics paths file
    graphic_paths_fn = None
    if len(graphic_paths) > 0:
        try:
            graphic_paths_fn = os.path.join(dirname, "graphic_paths.txt")
            with open(graphic_paths_fn, "w") as o_file:
                for path in sorted(graphic_paths):
                    print >>o_file, path
        except IOError:
            msg = "{}: cannot write summary".format(dirname)
            raise ProgramError(msg)

    # We create the automatic report
    logger.info("Generating automatic report")
    report_name = os.path.join(dirname, "automatic_report.tex")
    auto_report.create_report(
        dirname,
        report_name,
        project_name=args.report_number,
        steps=steps,
        descriptions=descriptions,
        graphic_paths_fn=graphic_paths_fn,
        long_descriptions=long_descriptions,
        summaries=latex_summaries,
        background=args.report_background,
        summary_fn=os.path.join(dirname, "results_summary.txt"),
        report_title=args.report_title,
        report_author=args.report_author,
        initial_files=os.path.join(dirname, "initial_files.txt"),
        final_files=os.path.join(dirname, "final_files.txt"),
        final_nb_markers="{:,d}".format(nb_markers),
        final_nb_samples="{:,d}".format(nb_samples),
        plink_version=plink_version,
    )


def run_duplicated_samples(in_prefix, in_type, out_prefix, base_dir, options):
    """Runs step1 (duplicated samples).

    :param in_prefix: the prefix of the input files.
    :param in_type: the type of the input files.
    :param out_prefix: the output prefix.
    :param base_dir: the output directory.
    :param options: the options needed.

    :type in_prefix: str
    :type in_type: str
    :type out_prefix: str
    :type base_dir: str
    :type options: list

    :returns: a tuple containing the prefix of the output files (the input
              prefix for the next script) and the type of the output files
              (``tfile``).

    This function calls the :py:mod:`pyGenClean.DupSamples.duplicated_samples`
    module. The required file type for this module is ``tfile``, hence the need
    to use the :py:func:`check_input_files` to check if the file input file
    type is the good one, or to create it if needed.

    """
    # Creating the output directory
    os.mkdir(out_prefix)

    # We know we need tfile
    required_type = "tfile"
    check_input_files(in_prefix, in_type, required_type)

    # We need to inject the name of the input file and the name of the output
    # prefix
    script_prefix = os.path.join(out_prefix, "dup_samples")
    options += ["--{}".format(required_type), in_prefix,
                "--out", script_prefix]

    # We run the script
    try:
        duplicated_samples.main(options)
    except duplicated_samples.ProgramError as e:
        msg = "duplicated_samples: {}".format(e)
        raise ProgramError(msg)

    # Reading the number of duplicated samples
    duplicated_count = defaultdict(int)
    if os.path.isfile(script_prefix + ".duplicated_samples.tfam"):
        with open(script_prefix + ".duplicated_samples.tfam", "r") as i_file:
            duplicated_count = Counter([
                tuple(createRowFromPlinkSpacedOutput(line)[:2])
                for line in i_file
            ])

    # Counting the number of zeroed out genotypes per duplicated sample
    zeroed_out = defaultdict(int)
    if os.path.isfile(script_prefix + ".zeroed_out"):
        with open(script_prefix + ".zeroed_out", "r") as i_file:
            zeroed_out = Counter([
                tuple(line.rstrip("\r\n").split("\t")[:2])
                for line in i_file.read().splitlines()[1:]
            ])
    nb_zeroed_out = sum(zeroed_out.values())

    # Checking the not good enough samples
    not_good_enough = set()
    if os.path.isfile(script_prefix + ".not_good_enough"):
        with open(script_prefix + ".not_good_enough", "r") as i_file:
            not_good_enough = {
                tuple(line.rstrip("\r\n").split("\t")[:4])
                for line in i_file.read().splitlines()[1:]
            }

    # Checking which samples were chosen
    chosen_sample = set()
    if os.path.isfile(script_prefix + ".chosen_samples.info"):
        with open(script_prefix + ".chosen_samples.info", "r") as i_file:
            chosen_sample = {
                tuple(line.rstrip("\r\n").split("\t"))
                for line in i_file.read().splitlines()[1:]
            }

    # Finding if some 'not_good_enough' samples were chosen
    not_good_still = {s[2:] for s in chosen_sample & not_good_enough}

    # We create a LaTeX summary
    latex_file = os.path.join(script_prefix + ".summary.tex")
    try:
        with open(latex_file, "w") as o_file:
            print >>o_file, latex_template.subsection(
                duplicated_samples.pretty_name
            )
            text = (
                "A total of {:,d} duplicated sample{} {} found.".format(
                    len(duplicated_count),
                    "s" if len(duplicated_count) > 1 else "",
                    "were" if len(duplicated_count) > 1 else "was",
                )
            )
            print >>o_file, latex_template.wrap_lines(text)

            if len(duplicated_count) > 0:
                text = (
                    "While merging duplicates, a total of {:,d} genotype{} {} "
                    "zeroed out. A total of {:,d} sample{} {} found to be not "
                    "good enough for duplicate completion.".format(
                        nb_zeroed_out,
                        "s" if nb_zeroed_out > 1 else "",
                        "were" if nb_zeroed_out > 1 else "was",
                        len(not_good_enough),
                        "s" if len(not_good_enough) > 1 else "",
                        "were" if len(not_good_enough) > 1 else "was",
                    )
                )
                print >>o_file, latex_template.wrap_lines(text)

                table_label = re.sub(
                    r"[/\\]",
                    "_",
                    script_prefix,
                ) + "_dup_samples"
                text = (
                    r"Table~\ref{" + table_label + "} summarizes the number "
                    "of each duplicated sample with some characteristics."
                )
                print >>o_file, latex_template.wrap_lines(text)

                if len(not_good_still) > 0:
                    text = latex_template.textbf(
                        "There {} {:,d} sample{} that {} not good due to low "
                        "completion or concordance, but {} still selected as "
                        "the best duplicate (see Table~{}).".format(
                            "were" if len(not_good_still) > 1 else "was",
                            len(not_good_still),
                            "s" if len(not_good_still) > 1 else "",
                            "were" if len(not_good_still) > 1 else "was",
                            "were" if len(not_good_still) > 1 else "was",
                            r"~\ref{" + table_label + "}",
                        )
                    )
                    print >>o_file, latex_template.wrap_lines(text)

                # Getting the template
                longtable_template = latex_template.jinja2_env.get_template(
                    "longtable_template.tex",
                )

                # The table caption
                table_caption = (
                    "Summary of the {:,d} duplicated sample{}. The number of "
                    "duplicates and the total number of zeroed out genotypes "
                    "are shown.".format(
                        len(duplicated_count),
                        "s" if len(duplicated_count) > 1 else "",
                    )
                )

                if len(not_good_still) > 0:
                    table_caption += (
                        " A total of {:,d} sample{} (highlighted) {} not good "
                        "enough for completion, but {} chosen as the best "
                        "duplicate, and {} still in the final "
                        "dataset).".format(
                            len(not_good_still),
                            "s" if len(not_good_still) > 1 else "",
                            "were" if len(not_good_still) > 1 else "was",
                            "were" if len(not_good_still) > 1 else "was",
                            "are" if len(not_good_still) > 1 else "is",
                        )
                    )

                duplicated_samples_list = duplicated_count.most_common()
                print >>o_file, longtable_template.render(
                    table_caption=table_caption,
                    table_label=table_label,
                    nb_col=4,
                    col_alignments="llrr",
                    text_size="scriptsize",
                    header_data=[("FID", 1), ("IID", 1), ("Nb Duplicate", 1),
                                 ("Nb Zeroed", 1)],
                    tabular_data=[
                        [latex_template.sanitize_tex(fid),
                         latex_template.sanitize_tex(iid),
                         "{:,d}".format(nb),
                         "{:,d}".format(zeroed_out[(fid, iid)])]
                        for (fid, iid), nb in duplicated_samples_list
                    ],
                    highlighted=[
                        (fid, iid) in not_good_still
                        for fid, iid in [i[0] for i in duplicated_samples_list]
                    ],
                )

    except IOError:
        msg = "{}: cannot write LaTeX summary".format(latex_file)
        raise ProgramError(msg)

    # Writing the summary results
    with open(os.path.join(base_dir, "results_summary.txt"), "a") as o_file:
        print >>o_file, "# {}".format(script_prefix)
        counter = Counter(duplicated_count.values()).most_common()
        if counter:
            print >>o_file, "Number of replicated samples"
        else:
            print >>o_file, "Number of replicated samples\t0"
        for rep_type, rep_count in counter:
            print >>o_file, "  - x{}\t{:,d}\t\t-{:,d}".format(
                rep_type,
                rep_count,
                (rep_count * rep_type) - rep_count,
            )
        print >>o_file, ("Poorly chosen replicated "
                         "samples\t{:,d}".format(len(not_good_still)))
        print >>o_file, "---"

    # We know this step does produce a new data set (tfile), so we return it
    return _StepResult(
        next_file=os.path.join(out_prefix, "dup_samples.final"),
        next_file_type="tfile",
        latex_summary=latex_file,
        description=duplicated_samples.desc,
        long_description=duplicated_samples.long_desc,
        graph_path=None,
    )


def run_duplicated_snps(in_prefix, in_type, out_prefix, base_dir, options):
    """Runs step2 (duplicated snps).

    :param in_prefix: the prefix of the input files.
    :param in_type: the type of the input files.
    :param out_prefix: the output prefix.
    :param base_dir: the output directory.
    :param options: the options needed.

    :type in_prefix: str
    :type in_type: str
    :type out_prefix: str
    :type base_dir: str
    :type options: list

    :returns: a tuple containing the prefix of the output files (the input
              prefix for the next script) and the type of the output files
              (``tfile``).

    This function calls the :py:mod:`pyGenClean.DupSNPs.duplicated_snps`
    module. The required file type for this module is ``tfile``, hence the need
    to use the :py:func:`check_input_files` to check if the file input file
    type is the good one, or to create it if needed.

    .. note::
        This function creates a ``map`` file, needed for the
        :py:mod:`pyGenClean.DupSNPs.duplicated_snps` module.

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
                    row = createRowFromPlinkSpacedOutput(line)
                    print >>outputFile, "\t".join(row[:4])
        except IOError:
            msg = "{}: no such file".format(in_prefix + ".tped")
            raise ProgramError(msg)
        outputFile.close()

    # We need to inject the name of the input file and the name of the output
    # prefix
    script_prefix = os.path.join(out_prefix, "dup_snps")
    options += ["--{}".format(required_type), in_prefix,
                "--out", script_prefix]

    # We run the script
    try:
        duplicated_snps.main(options)
    except duplicated_snps.ProgramError as e:
        msg = "duplicated_snps: {}".format(e)
        raise ProgramError(msg)

    # Reading the number of duplicated markers
    duplicated_count = defaultdict(int)
    if os.path.isfile(script_prefix + ".duplicated_snps.tped"):
        with open(script_prefix + ".duplicated_snps.tped", "r") as i_file:
            duplicated_count = Counter(
                (i[0], i[3]) for i in [
                    tuple(createRowFromPlinkSpacedOutput(line)[:4])
                    for line in i_file
                ]
            )

    # Counting the number of zeroed out genotypes per duplicated markers
    zeroed_out = defaultdict(int)
    if os.path.isfile(script_prefix + ".zeroed_out"):
        with open(script_prefix + ".zeroed_out", "r") as i_file:
            zeroed_out = Counter([
                tuple(line.rstrip("\r\n").split("\t")[:2])
                for line in i_file.read().splitlines()[1:]
            ])
    nb_zeroed_out = sum(zeroed_out.values())

    # Checking the not good enough markers
    not_good_enough = set()
    if os.path.isfile(script_prefix + ".not_good_enough"):
        with open(script_prefix + ".not_good_enough", "r") as i_file:
            not_good_enough = {
                line.rstrip("\r\n").split("\t")[0]
                for line in i_file.read().splitlines()[1:]
            }

    # Checking which markers were chosen
    chosen_markers = set()
    if os.path.isfile(script_prefix + ".chosen_snps.info"):
        with open(script_prefix + ".chosen_snps.info", "r") as i_file:
            chosen_markers = set(i_file.read().splitlines())

    # Finding if some 'not_good_enough' samples were chosen
    not_good_still = chosen_markers & not_good_enough

    # Adding the 'not chosen markers' to the list of excluded markers
    removed_markers = set()
    o_filename = os.path.join(base_dir, "excluded_markers.txt")
    if os.path.isfile(script_prefix + ".removed_duplicates"):
        with open(script_prefix + ".removed_duplicates", "r") as i_file:
            removed_markers = set(i_file.read().splitlines())
            with open(o_filename, "a") as o_file:
                for marker_id in removed_markers:
                    print >>o_file, marker_id + "\t" + "removed duplicate"

    # Writing the summary results
    total_remaining = 0
    with open(os.path.join(base_dir, "results_summary.txt"), "a") as o_file:
        print >>o_file, "# {}".format(script_prefix)
        rep_counter = Counter(duplicated_count.values()).most_common()
        if rep_counter:
            print >>o_file, "Number of replicated markers"
        else:
            print >>o_file, "Number of replicated markers\t0"
        total_nb_removed_rep = 0
        for rep_type, rep_count in rep_counter:
            nb_removed_rep = (rep_count * rep_type) - rep_count
            print >>o_file, "  - x{}\t{:,d}\t-{:,d}".format(
                rep_type,
                rep_count,
                nb_removed_rep,
            )
            total_nb_removed_rep += nb_removed_rep
        total_remaining = total_nb_removed_rep - len(removed_markers)
        print >>o_file, (
            "Number of replicated markers kept\t{nb:,d}\t+{nb:,d}".format(
                nb=total_remaining,
            )
        )
        print >>o_file, ("Poorly chosen replicated markers\t"
                         "{nb:,d}".format(nb=len(not_good_still)))
        print >>o_file, ("Final number of excluded markers\t"
                         "{nb:,d}".format(nb=len(removed_markers)))
        print >>o_file, "---"

    # We create a LaTeX summary
    latex_file = os.path.join(script_prefix + ".summary.tex")
    try:
        with open(latex_file, "w") as o_file:
            print >>o_file, latex_template.subsection(
                duplicated_snps.pretty_name
            )

            text = (
                "A total of {:,d} duplicated marker{} {} found.".format(
                    len(duplicated_count),
                    "s" if len(duplicated_count) > 1 else "",
                    "were" if len(duplicated_count) > 1 else "was",
                )
            )
            print >>o_file, latex_template.wrap_lines(text)

            if len(duplicated_count) > 0:
                text = (
                    "While merging duplicates, a total of {:,d} genotype{} {} "
                    "zeroed out. A total of {:,d} marker{} {} found to be not "
                    "good enough for duplicate completion.".format(
                        nb_zeroed_out,
                        "s" if nb_zeroed_out > 1 else "",
                        "were" if nb_zeroed_out > 1 else "was",
                        len(not_good_enough),
                        "s" if len(not_good_enough) > 1 else "",
                        "were" if len(not_good_enough) > 1 else "was",
                    )
                )
                print >>o_file, latex_template.wrap_lines(text)

                text = (
                    "A total of {:,d} marker{} {} excluded while creating the "
                    "final dataset.".format(
                        len(removed_markers),
                        "s" if len(removed_markers) > 1 else "",
                        "were" if len(removed_markers) > 1 else "was",
                    )
                )
                print >>o_file, latex_template.wrap_lines(text)

                if total_remaining > 0:
                    text = latex_template.textbf(
                        "In total, {:,d} maker{} {} not merged for different "
                        "reasons (low completion rate, discordant allele, "
                        "discordant MAF, etc) and {} still present in the "
                        "dataset.".format(
                            total_remaining,
                            "s" if total_remaining > 1 else "",
                            "were" if total_remaining > 1 else "was",
                            "are" if total_remaining > 1 else "is",
                        )
                    )
                    print >>o_file, latex_template.wrap_lines(text)

                if len(not_good_still) > 0:
                    start = "A total of"
                    end = " and {} still present in the final dataset.".format(
                        "are" if len(not_good_still) > 1 else "is",
                    )
                    if total_remaining > 0:
                        start = "Out of these,"
                        end = "."
                    text = latex_template.textbf(
                        start + " {:,d} marker{} {} not good enough for "
                        "completion, but {} still selected as the best "
                        "duplicate{}".format(
                            len(not_good_still),
                            "s" if len(not_good_still) > 1 else "",
                            "were" if len(not_good_still) > 1 else "was",
                            "were" if len(not_good_still) > 1 else "was",
                            end,
                        )
                    )
                    print >>o_file, latex_template.wrap_lines(text)

    except IOError:
        msg = "{}: cannot write LaTeX summary".format(latex_file)
        raise ProgramError(msg)

    # We know this step does produce a new data set (tfile), so we return it
    return _StepResult(
        next_file=os.path.join(out_prefix, "dup_snps.final"),
        next_file_type="tfile",
        latex_summary=latex_file,
        description=duplicated_snps.desc,
        long_description=duplicated_snps.long_desc,
        graph_path=None,
    )


def run_noCall_hetero_snps(in_prefix, in_type, out_prefix, base_dir, options):
    """Runs step 3 (clean no call and hetero).

    :param in_prefix: the prefix of the input files.
    :param in_type: the type of the input files.
    :param out_prefix: the output prefix.
    :param base_dir: the output directory.
    :param options: the options needed.

    :type in_prefix: str
    :type in_type: str
    :type out_prefix: str
    :type base_dir: str
    :type options: list

    :returns: a tuple containing the prefix of the output files (the input
              prefix for the next script) and the type of the output files
              (``tfile``).

    This function calls the
    :py:mod:`pyGenClean.NoCallHetero.clean_noCall_hetero_snps` module. The
    required file type for this module is ``tfile``, hence the need to use the
    :py:func:`check_input_files` to check if the file input file type is the
    good one, or to create it if needed.

    """
    # Creating the output directory
    os.mkdir(out_prefix)

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
    with open(o_filename, "a") as o_file:
        # The first file
        i_filename = script_prefix + ".allFailed"
        if os.path.isfile(i_filename):
            with open(i_filename, "r") as i_file:
                for line in i_file:
                    nb_all_failed += 1
                    print >>o_file, line.rstrip("\r\n") + "\tall failed"

        # The second file
        i_filename = os.path.join(script_prefix + ".allHetero")
        if os.path.isfile(i_filename):
            with open(i_filename, "r") as i_file:
                for line in i_file:
                    nb_all_hetero += 1
                    print >>o_file, line.rstrip("\r\n") + "\tall hetero"

    # We write a LaTeX summary
    latex_file = os.path.join(script_prefix + ".summary.tex")
    try:
        with open(latex_file, "w") as o_file:
            print >>o_file, latex_template.subsection(
                noCall_hetero_snps.pretty_name,
            )
            text = (
                "After scrutiny, {:,d} marker{} {} excluded from the "
                "dataset because of a call rate of 0. Also, {:,d} marker{} "
                "{} excluded from the dataset because all samples were "
                "heterozygous (excluding the mitochondrial "
                "chromosome)".format(nb_all_failed,
                                     "s" if nb_all_failed > 1 else "",
                                     "were" if nb_all_failed > 1 else "was",
                                     nb_all_hetero,
                                     "s" if nb_all_hetero > 1 else "",
                                     "were" if nb_all_hetero > 1 else "was")
            )
            print >>o_file, latex_template.wrap_lines(text, 80)

    except IOError:
        msg = "{}: cannot write LaTeX summary".format(latex_file)
        raise ProgramError(msg)

    # Writing the summary results
    with open(os.path.join(base_dir, "results_summary.txt"), "a") as o_file:
        print >>o_file, "# {}".format(script_prefix)
        print >>o_file, ("Number of completely failed markers\t"
                         "{nb:,d}\t-{nb:,d}".format(nb=nb_all_failed))
        print >>o_file, "---"
        print >>o_file, ("Number of all heterozygous markers\t"
                         "{nb:,d}\t-{nb:,d}".format(nb=nb_all_hetero))
        print >>o_file, "---"

    # We know this step does produce a new data set (tfile), so we return it
    # along with the report name
    return _StepResult(
        next_file=os.path.join(out_prefix, "clean_noCall_hetero"),
        next_file_type="tfile",
        latex_summary=latex_file,
        description=noCall_hetero_snps.desc,
        long_description=noCall_hetero_snps.long_desc,
        graph_path=None,
    )


def run_sample_missingness(in_prefix, in_type, out_prefix, base_dir, options):
    """Runs step4 (clean mind).

    :param in_prefix: the prefix of the input files.
    :param in_type: the type of the input files.
    :param out_prefix: the output prefix.
    :param base_dir: the output directory.
    :param options: the options needed.

    :type in_prefix: str
    :type in_type: str
    :type out_prefix: str
    :type base_dir: str
    :type options: list

    :returns: a tuple containing the prefix of the output files (the input
              prefix for the next script) and the type of the output files
              (``bfile``).

    This function calls the
    :py:mod:`pyGenClean.SampleMissingness.sample_missingness` module. The
    required file type for this module is either a ``bfile`` or a ``tfile``,
    hence the need to use the :py:func:`check_input_files` to check if the file
    input file type is the good one, or to create it if needed.

    """
    # Creating the output directory
    os.mkdir(out_prefix)

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
    if desc.endswith("."):
        desc = desc[:-1] + r" (${}={}$).".format(latex_template.texttt("mind"),
                                                 mind_value)
    long_description = sample_missingness.long_desc.format(
        success_rate="{:.1f}".format((1-float(mind_value)) * 100),
    )

    # We want to save in a file the samples that were removed
    # There is one file to look at, which contains only one row, the name of
    # the samples:
    #     - prefix.irem (file will exists only if samples were removed)
    nb_samples = 0
    o_filename = os.path.join(base_dir, "excluded_samples.txt")
    i_filename = script_prefix + ".irem"
    if os.path.isfile(i_filename):
        # True, so sample were removed
        with open(i_filename, "r") as i_file, open(o_filename, "a") as o_file:
            for line in i_file:
                nb_samples += 1
                print >>o_file, line.rstrip("\r\n") + "\tmind {}".format(
                    mind_value,
                )

    # We write a LaTeX summary
    latex_file = os.path.join(script_prefix + ".summary.tex")
    try:
        with open(latex_file, "w") as o_file:
            print >>o_file, latex_template.subsection(
                sample_missingness.pretty_name,
            )
            text = ("Using a {} threshold of {} ({} keeping only samples with "
                    r"a missing rate $\leq {}$), {:,d} sample{} {} excluded "
                    "from the dataset.".format(
                        latex_template.texttt("mind"),
                        mind_value,
                        latex_template.textit("i.e."),
                        mind_value,
                        nb_samples,
                        "s" if nb_samples > 1 else "",
                        "were" if nb_samples > 1 else "was",
                    ))
            print >>o_file, latex_template.wrap_lines(text)

    except IOError:
        msg = "{}: cannot write LaTeX summary".format(latex_file)
        raise ProgramError(msg)

    # Writing the summary results
    with open(os.path.join(base_dir, "results_summary.txt"), "a") as o_file:
        print >>o_file, "# {}".format(script_prefix)
        print >>o_file, ("Number of samples with missing rate higher "
                         "than {t}\t{nb:,d}\t\t-{nb:,d}".format(
                            t=mind_value,
                            nb=nb_samples,
                         ))
        print >>o_file, "---"

    # We know this step does produce a new data set (bfile), so we return it
    return _StepResult(
        next_file=os.path.join(out_prefix, "clean_mind"),
        next_file_type="bfile",
        latex_summary=latex_file,
        description=desc,
        long_description=long_description,
        graph_path=None,
    )


def run_snp_missingness(in_prefix, in_type, out_prefix, base_dir, options):
    """Run step5 (clean geno).

    :param in_prefix: the prefix of the input files.
    :param in_type: the type of the input files.
    :param out_prefix: the output prefix.
    :param base_dir: the output directory.
    :param options: the options needed.

    :type in_prefix: str
    :type in_type: str
    :type out_prefix: str
    :type base_dir: str
    :type options: list

    :returns: a tuple containing the prefix of the output files (the input
              prefix for the next script) and the type of the output files
              (``bfile``).

    This function calls the
    :py:mod:`pyGenClean.MarkerMissingness.snp_missingness` module. The required
    file type for this module is ``bfile``, hence the need to use the
    :py:func:`check_input_files` to check if the file input file type is the
    good one, or to create it if needed.

    """
    # Creating the output directory
    os.mkdir(out_prefix)

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
    if desc.endswith("."):
        desc = desc[:-1] + r" (${}={}$).".format(latex_template.texttt("geno"),
                                                 geno_value)
    long_description = sample_missingness.long_desc.format(
        success_rate="{:.1f}".format((1-float(geno_value)) * 100),
    )

    # We want to save in a file the samples that were removed
    # There is one file to look at, which contains only one row, the name of
    # the samples:
    #     - prefix.removed_snps
    nb_markers = 0
    o_filename = os.path.join(base_dir, "excluded_markers.txt")
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

    # We write a LaTeX summary
    latex_file = os.path.join(script_prefix + ".summary.tex")
    try:
        with open(latex_file, "w") as o_file:
            print >>o_file, latex_template.subsection(
                snp_missingness.pretty_name,
            )
            text = ("Using a {} threshold of {} ({} keeping only markers with "
                    r"a missing rate $\leq {}$), {:,d} marker{} {} excluded "
                    "from the dataset.".format(
                        latex_template.texttt("geno"),
                        geno_value,
                        latex_template.textit("i.e."),
                        geno_value, nb_markers,
                        "s" if nb_markers > 1 else "",
                        "were" if nb_markers > 1 else "was",
                    ))
            print >>o_file, latex_template.wrap_lines(text)

    except IOError:
        msg = "{}: cannot write LaTeX summary".format(latex_file)
        raise ProgramError(msg)

    # Writing the summary results
    with open(os.path.join(base_dir, "results_summary.txt"), "a") as o_file:
        print >>o_file, "# {}".format(script_prefix)
        print >>o_file, ("Number of markers with missing rate higher "
                         "than {t}\t{nb:,d}\t-{nb:,d}".format(
                            t=geno_value,
                            nb=nb_markers,
                         ))
        print >>o_file, "---"

    # We know this step does produce a new data set (bfile), so we return it
    return _StepResult(
        next_file=os.path.join(out_prefix, "clean_geno"),
        next_file_type="bfile",
        latex_summary=latex_file,
        description=desc,
        long_description=long_description,
        graph_path=None,
    )


def run_contamination(in_prefix, in_type, out_prefix, base_dir, options):
    """Runs the contamination check for samples.

    :param in_prefix: the prefix of the input files.
    :param in_type: the type of the input files.
    :param out_prefix: the output prefix.
    :param base_dir: the output directory.
    :param options: the options needed.

    :type in_prefix: str
    :type in_type: str
    :type out_prefix: str
    :type base_dir: str
    :type options: list

    :returns: a tuple containing the prefix of the output files (the input
              prefix for the next script) and the type of the output files
              (``bfile``).

    """
    # Creating the output directory
    os.mkdir(out_prefix)

    # We know we need a bfile
    required_type = "bfile"
    check_input_files(in_prefix, in_type, required_type)

    # We need to inject the name of the input file and the name of the output
    # prefix
    script_prefix = os.path.join(out_prefix, "contamination")
    options += ["--{}".format(required_type), in_prefix,
                "--out", script_prefix]

    # We run the script
    try:
        contamination.main(options)
    except contamination.ProgramError as e:
        msg = "contamination: {}".format(e)
        raise ProgramError(msg)

    # Counting the number of markers used for contamination
    nb_autosomal = 0
    with open(script_prefix + ".to_extract", "r") as i_file:
        nb_autosomal = len(i_file.read().splitlines())

    # Reading the "contamination" file
    nb_tested_samples = 0
    contaminated_table = []
    nb_contaminated_samples = 0
    with open(script_prefix + ".bafRegress", "r") as i_file:
        header = None
        for i, line in enumerate(i_file):
            row = line.rstrip("\r\n").split("\t")
            if i == 0:
                contaminated_table.append(("sample", "estimate$^1$",
                                           "stderr$^1$", "tval", "pval",
                                           "callrate", "Nhom$^2$"))
                header = {name: i for i, name in enumerate(row)}

                for name in ("sample", "estimate", "stderr", "tval", "pval",
                             "callrate", "Nhom"):
                    if name not in header:
                        raise ProgramError("{}: missing column {}".format(
                            script_prefix + ".bafRegress",
                            name,
                        ))
                continue

            # One more sample
            nb_tested_samples += 1

            # Reading the data
            estimate = float(row[header["estimate"]])
            if estimate > 0.01:
                # This might be a contaminated samples
                contaminated_table.append((
                    latex_template.sanitize_tex(row[header["sample"]]),
                    "{:.4f}".format(float(row[header["estimate"]])),
                    "{:.4f}".format(float(row[header["stderr"]])),
                    "{:.4f}".format(float(row[header["tval"]])),
                    latex_template.format_numbers("{:.3e}".format(
                        float(row[header["pval"]]),
                    )),
                    "{:.4f}".format(float(row[header["callrate"]])),
                    "{:,d}".format(int(row[header["Nhom"]])),
                ))
                nb_contaminated_samples += 1

    # We write a LaTeX summary
    latex_file = os.path.join(script_prefix + ".summary.tex")
    try:
        with open(latex_file, "w") as o_file:
            print >>o_file, latex_template.subsection(
                contamination.pretty_name,
            )
            text = ("A total of {:,d} sample{} {} analyzed for contamination "
                    r"using \textit{bafRegress}\cite{bafRegress}. The "
                    "analysis was performed using {:,d} autosomal marker{}. "
                    "Using a threshold of 0.01, {:,d} sample{} {} estimated "
                    "to be contaminated.".format(
                        nb_tested_samples,
                        "s" if nb_tested_samples > 1 else "",
                        "were" if nb_tested_samples > 1 else "was",
                        nb_autosomal,
                        "s" if nb_autosomal > 1 else "",
                        nb_contaminated_samples,
                        "s" if nb_contaminated_samples > 1 else "",
                        "were" if nb_contaminated_samples > 1 else "was",
                        bafRegress="{bafRegress}",
                    ))
            print >>o_file, latex_template.wrap_lines(text)

            if nb_contaminated_samples > 0:
                label = re.sub(r"[/\\]", "_", script_prefix)
                text = (
                    r"Table~\ref{" + label + "} list all the samples that "
                    r"were estimated to be contaminated (\textit{i.e.} with "
                    r"an estimate $>0.01$)."
                )
                print >>o_file, latex_template.wrap_lines(text)

                # Getting the template
                longtable_template = latex_template.jinja2_env.get_template(
                    "longtable_template.tex",
                )

                # Rendering
                print >>o_file, longtable_template.render(
                    table_caption=r"List of all possible contaminated samples "
                                  r"(\textit{i.e.} with an estimate computed "
                                  r"by \textit{bafRegress} $>0.01$).",
                    table_label=label,
                    nb_col=len(contaminated_table[1]),
                    col_alignments="lrrrrrr",
                    text_size="scriptsize",
                    header_data=zip(contaminated_table[0],
                                    [1 for i in contaminated_table[0]]),
                    tabular_data=sorted(
                        contaminated_table[1:],
                        key=lambda item: item[0],
                    ),
                    footnotes=[
                        "$^1$Contamination estimate (with standard error)",
                        "$^2$Number of homozygous genotypes used in the model",
                    ],
                )

    except IOError:
        msg = "{}: cannot write LaTeX summary".format(latex_file)
        raise ProgramError(msg)

    # Writing the summary results
    with open(os.path.join(base_dir, "results_summary.txt"), "a") as o_file:
        print >>o_file, "# {}".format(script_prefix)
        print >>o_file, ("Number of possibly contaminated "
                         "samples\t{:,d}".format(nb_contaminated_samples))
        print >>o_file, "---"

    # We know this step doesn't produce an new data set, so we return the old
    # prefix and the old in_type
    return _StepResult(
        next_file=in_prefix,
        next_file_type=required_type,
        latex_summary=latex_file,
        description=contamination.desc,
        long_description=contamination.long_desc,
        graph_path=None,
    )


def run_sex_check(in_prefix, in_type, out_prefix, base_dir, options):
    """Runs step6 (sexcheck).

    :param in_prefix: the prefix of the input files.
    :param in_type: the type of the input files.
    :param out_prefix: the output prefix.
    :param base_dir: the output directory.
    :param options: the options needed.

    :type in_prefix: str
    :type in_type: str
    :type out_prefix: str
    :type base_dir: str
    :type options: list

    :returns: a tuple containing the prefix of the output files (the input
              prefix for the next script) and the type of the output files
              (``bfile``).

    This function calls the :py:mod:`pyGenClean.SexCheck.sex_check` module. The
    required file type for this module is ``bfile``, hence the need to use the
    :py:func:`check_input_files` to check if the file input file type is the
    good one, or to create it if needed.

    .. note::
        The :py:mod:`pyGenClean.SexCheck.sex_check` module doesn't return
        usable output files. Hence, this function returns the input file prefix
        and its type.

    """
    # Creating the output directory
    os.mkdir(out_prefix)

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

    # Reading the hetero file on X
    hetero = {}
    if os.path.isfile(script_prefix + ".chr23_recodeA.raw.hetero"):
        with open(script_prefix + ".chr23_recodeA.raw.hetero", "r") as i_file:
            header = {
                name: i for i, name in
                enumerate(createRowFromPlinkSpacedOutput(i_file.readline()))
            }
            for required_col in ("PED", "ID", "HETERO"):
                if required_col not in header:
                    msg = "{}: no column named {}".format(
                        script_prefix + ".chr23_recodeA.raw.hetero",
                        required_col,
                    )
                    raise ProgramError(msg)

            # Reading the data
            for line in i_file:
                row = line.rstrip("\r\n").split("\t")
                famid = row[header["PED"]]
                indid = row[header["ID"]]

                # Formatting the hetero value
                het = None
                try:
                    het = "{:.4f}".format(float(row[header["HETERO"]]))
                except:
                    het = "N/A"

                hetero[(famid, indid)] = het

    # Reading the number of no call on Y
    nb_no_call = {}
    if os.path.isfile(script_prefix + ".chr24_recodeA.raw.noCall"):
        with open(script_prefix + ".chr24_recodeA.raw.noCall", "r") as i_file:
            header = {
                name: i for i, name in
                enumerate(createRowFromPlinkSpacedOutput(i_file.readline()))
            }
            for required_col in ("PED", "ID", "nbGeno", "nbNoCall"):
                if required_col not in header:
                    msg = "{}: no column named {}".format(
                        script_prefix + ".chr24_recodeA.raw.noCall",
                        required_col,
                    )
                    raise ProgramError(msg)

            # Reading the data
            for line in i_file:
                row = line.rstrip("\r\n").split("\t")
                famid = row[header["PED"]]
                indid = row[header["ID"]]

                # Getting the statistics
                nb_geno = row[header["nbGeno"]]
                nb_nocall = row[header["nbNoCall"]]

                percent = None
                try:
                    percent = "{:.4f}".format(
                        float(nb_nocall) / float(nb_geno),
                    )
                except:
                    percent = "N/A"
                nb_no_call[(famid, indid)] = percent

    # Reading the problem file to gather statistics. Note that dataset without
    # problem will only have the header line (and no data)
    nb_problems = 0
    table = []
    nb_no_genetic = 0
    nb_discordant = 0
    with open(script_prefix + ".list_problem_sex", "r") as i_file:
        # Reading the header
        header = i_file.readline().rstrip("\r\n").split("\t")
        table.append(header)
        header = {name: i for i, name in enumerate(header)}
        for required_col in ("FID", "IID", "SNPSEX"):
            if required_col not in header:
                msg = "{}: no column named {}".format(
                    script_prefix + ".list_problem_sex",
                    required_col,
                )
                raise ProgramError(msg)

        # Adding the missing column name
        table[-1].append("HET")
        table[-1].append(r"\%NOCALL")

        # Reading the rest of the data
        for line in i_file:
            nb_problems += 1

            # Creating the row
            row = line.rstrip("\r\n").split("\t")

            # Counting
            if row[header["SNPSEX"]] == "0":
                nb_no_genetic += 1
            else:
                nb_discordant += 1

            table.append([
                latex_template.sanitize_tex(row[header[name]])
                for name in ("FID", "IID", "PEDSEX", "SNPSEX", "STATUS", "F")
            ])
            table[-1].append(
                hetero.get((row[header["FID"]], row[header["IID"]]), "N/A"),
            )
            table[-1].append(
                nb_no_call.get((row[header["FID"]], row[header["IID"]]), "N/A")
            )

    # Getting the value for the maleF option
    male_f = sex_check.parser.get_default("maleF")
    if "--maleF" in options:
        male_f = options[options.index("--maleF") + 1]

    # Getting the value for the femaleF option
    female_f = sex_check.parser.get_default("femaleF")
    if "--femaleF" in options:
        female_f = options[options.index("--femaleF") + 1]

    # We write a LaTeX summary
    latex_file = os.path.join(script_prefix + ".summary.tex")
    graphics_paths = set()
    try:
        with open(latex_file, "w") as o_file:
            print >>o_file, latex_template.subsection(sex_check.pretty_name)
            text = (
                "Using $F$ thresholds of {male_f} and {female_f} for males "
                "and females respectively, {nb_problems:,d} sample{plural} "
                "had gender problem according to Plink.".format(
                    male_f=male_f,
                    female_f=female_f,
                    nb_problems=nb_problems,
                    plural="s" if nb_problems > 1 else "",
                )
            )
            print >>o_file, latex_template.wrap_lines(text)

            # The float template
            float_template = latex_template.jinja2_env.get_template(
                "float_template.tex",
            )

            if nb_problems > 0:
                # The label and text for the table
                table_label = re.sub(
                    r"[/\\]",
                    "_",
                    script_prefix,
                ) + "_problems"
                text = (
                    r"Table~\ref{" + table_label + "} summarizes the gender "
                    "problems encountered during the analysis."
                )
                print >>o_file, latex_template.wrap_lines(text)

                # Getting the template
                longtable_template = latex_template.jinja2_env.get_template(
                    "longtable_template.tex",
                )

                # Rendering
                print >>o_file, longtable_template.render(
                    table_caption="Summarization of the gender problems "
                                  "encountered during Plink's analysis. "
                                  "HET is the heterozygosity rate on the X "
                                  r"chromosome. \%NOCALL is the percentage of "
                                  "no calls on the Y chromosome.",
                    table_label=table_label,
                    nb_col=len(table[1]),
                    col_alignments="llrrlrrrr",
                    text_size="scriptsize",
                    header_data=zip(table[0], [1 for i in table[0]]),
                    tabular_data=sorted(table[1:], key=lambda item: item[1]),
                )

            # Getting the templates
            graphic_template = latex_template.jinja2_env.get_template(
                "graphics_template.tex",
            )

            # If there is a figure, we add it here
            if os.path.isfile(script_prefix + ".png"):
                # Adding the figure
                figure_label = re.sub(r"[/\\]", "_", script_prefix)
                text = (
                    r"Figure~\ref{" + figure_label + r"} shows the $\bar{y}$ "
                    r"intensities versus the $\bar{x}$ intensities for each "
                    "samples. Problematic samples are shown using triangles."
                )
                print >>o_file, latex_template.wrap_lines(text)

                # Getting the paths
                graphics_path, path = os.path.split(script_prefix + ".png")
                graphics_path = os.path.relpath(graphics_path, base_dir)

                print >>o_file, float_template.render(
                    float_type="figure",
                    float_placement="H",
                    float_caption="Gender check using Plink. Mean $x$ and $y$ "
                                  "intensities are shown for each sample. "
                                  "Males are shown in blue, and females in "
                                  "red. Triangles show problematic samples "
                                  "(green for males, mauve for females). "
                                  "Unknown gender are shown in gray.",
                    float_label=figure_label,
                    float_content=graphic_template.render(
                        width=r"0.8\textwidth",
                        path=latex_template.sanitize_fig_name(path),
                    ),
                )
                # Adding the path where the graphic is
                graphics_paths.add(graphics_path)

            # If there is a 'sexcheck.LRR_BAF' directory, then there are LRR
            # and BAF plots.
            if os.path.isdir(script_prefix + ".LRR_BAF"):
                figures = glob(
                    os.path.join(script_prefix + ".LRR_BAF", "*.png"),
                )

                if len(figures) > 0:
                    # Getting the sample IDs
                    sample_ids = [
                        re.search(
                            "^baf_lrr_(\S+)_lrr_baf.png$",
                            os.path.basename(figure),
                        ) for figure in figures
                    ]
                    sample_ids = [
                        "unknown sample" if not sample else sample.group(1)
                        for sample in sample_ids
                    ]

                    # Sorting according to sample IDs
                    sorted_indexes = sorted(range(len(figures)),
                                            key=figures.__getitem__)
                    figures = [figures[i] for i in sorted_indexes]
                    sample_ids = [sample_ids[i] for i in sorted_indexes]

                    # Getting the labels
                    labels = [
                        re.sub(
                            r"[/\\]",
                            "_",
                            script_prefix + "_baf_lrr_" +
                            os.path.splitext(sample)[0],
                        ) for sample in sample_ids
                    ]

                    fig_1 = labels[0]
                    fig_2 = ""
                    if len(figures) > 1:
                        fig_2 = labels[-1]
                    text = (
                        "Figure" + ("s" if len(figures) > 1 else "") +
                        r"~\ref{" + fig_1 + "} " +
                        (r"to \ref{" + fig_2 + "} " if fig_2 else "") +
                        "show" + (" " if len(figures) > 1 else "s ") + "the "
                        "log R ratio and the B allele frequency versus the "
                        "position on chromosome X and Y for the problematic "
                        "sample{}.".format("s" if len(figures) > 1 else "")
                    )
                    print >>o_file, latex_template.wrap_lines(text)

                    zipped = zip(figures, sample_ids, labels)
                    for figure, sample_id, label in zipped:
                        sample_id = latex_template.sanitize_tex(sample_id)

                        # Getting the paths
                        graphics_path, path = os.path.split(figure)
                        graphics_path = os.path.relpath(graphics_path,
                                                        base_dir)

                        caption = (
                            "Plots showing the log R ratio and the B allele "
                            "frequency for chromosome X and Y (on the left "
                            "and right, respectively) for sample "
                            "{}.".format(sample_id)
                        )
                        print >>o_file, float_template.render(
                            float_type="figure",
                            float_placement="H",
                            float_caption=caption,
                            float_label=label,
                            float_content=graphic_template.render(
                                width=r"\textwidth",
                                path=latex_template.sanitize_fig_name(path),
                            ),
                        )
                # Adding the path where the graphic is
                graphics_paths.add(graphics_path)

    except IOError:
        msg = "{}: cannot write LaTeX summary".format(latex_file)
        raise ProgramError(msg)

    # Writing the summary results
    with open(os.path.join(base_dir, "results_summary.txt"), "a") as o_file:
        print >>o_file, "# {}".format(script_prefix)
        print >>o_file, "Number of samples with gender problem"
        print >>o_file, "  - no genetic gender\t{:,d}".format(nb_no_genetic)
        print >>o_file, "  - discordant gender\t{:,d}".format(nb_discordant)
        print >>o_file, "---"

    # We know this step does not produce a new data set, so we return the
    # original one
    return _StepResult(
        next_file=in_prefix,
        next_file_type=required_type,
        latex_summary=latex_file,
        description=sex_check.desc,
        long_description=sex_check.long_desc,
        graph_path=graphics_paths,
    )


def run_plate_bias(in_prefix, in_type, out_prefix, base_dir, options):
    """Runs step7 (plate bias).

    :param in_prefix: the prefix of the input files.
    :param in_type: the type of the input files.
    :param out_prefix: the output prefix.
    :param base_dir: the output directory.
    :param options: the options needed.

    :type in_prefix: str
    :type in_type: str
    :type out_prefix: str
    :type base_dir: str
    :type options: list

    :returns: a tuple containing the prefix of the output files (the input
              prefix for the next script) and the type of the output files
              (``bfile``).

    This function calls the :py:mod:`pyGenClean.PlateBias.plate_bias` module.
    The required file type for this module is ``bfile``, hence the need to use
    the :py:func:`check_input_files` to check if the file input file type is
    the good one, or to create it if needed.

    .. note::
        The :py:mod:`pyGenClean.PlateBias.plate_bias` module doesn't return
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
    script_prefix = os.path.join(out_prefix, "plate_bias")
    options += ["--{}".format(required_type), in_prefix,
                "--out", script_prefix]

    # We run the script
    try:
        plate_bias.main(options)
    except plate_bias.ProgramError as e:
        msg = "plate_bias: {}".format(e)
        raise ProgramError(msg)

    # The name of the summary file
    filename = script_prefix + ".significant_SNPs.summary"
    if not os.path.isfile(filename):
        raise ProgramError("{}: no such file".format(filename))

    # Getting the number of each plate
    plate_counter = None
    with open(filename, "r") as i_file:
        header = {
            name: i for i, name in
            enumerate(i_file.readline().rstrip("\r\n").split("\t"))
        }
        if "plate" not in header:
            msg = "{}: missing column plate".format(filename)
            raise ProgramError(msg)

        # Counting the number of markers for each plate
        plate_counter = Counter(
            line.rstrip("\r\n").split("\t")[header["plate"]] for line in i_file
        )

    # Creating the table
    table = [["plate name", "number of markers"]]
    for plate_name, number in plate_counter.most_common():
        table.append([
            latex_template.sanitize_tex(plate_name),
            "{:,d}".format(number),
        ])

    # The number of unique markers
    filename = script_prefix + ".significant_SNPs.txt"
    nb_markers = None
    with open(filename, "r") as i_file:
        nb_markers = len({line.rstrip("\r\n") for line in i_file})

    # Getting the p value threshold
    p_threshold = str(plate_bias.parser.get_default("pfilter"))
    if "--pfilter" in options:
        p_threshold = str(options[options.index("--pfilter") + 1])

    # We write a LaTeX summary
    latex_file = os.path.join(script_prefix + ".summary.tex")
    try:
        with open(latex_file, "w") as o_file:
            print >>o_file, latex_template.subsection(plate_bias.pretty_name)
            text = (
                "After performing the plate bias analysis using Plink, a "
                "total of {:,d} unique marker{} had a significant result "
                "({} a value less than {}).".format(
                    nb_markers,
                    "s" if nb_markers > 1 else "",
                    r"\textit{i.e.}",
                    latex_template.format_numbers(p_threshold),
                )
            )
            print >>o_file, latex_template.wrap_lines(text)

            if nb_markers > 0:
                table_label = re.sub(
                    r"[/\\]",
                    "_",
                    script_prefix,
                ) + "_plate_bias"
                text = (
                    r"Table~\ref{" + table_label + "} summarizes the plate "
                    "bias results."
                )
                print >>o_file, latex_template.wrap_lines(text)

                # Getting the template
                longtable_template = latex_template.jinja2_env.get_template(
                    "longtable_template.tex",
                )

                # The table caption
                table_caption = (
                    "Summary of the plate bias analysis performed by Plink. "
                    "For each plate, the number of significant marker{} is "
                    "shown (threshold of {}). The plates are sorted according "
                    "to the total number of significant results.".format(
                        "s" if nb_markers > 1 else "",
                        latex_template.format_numbers(p_threshold),
                    )
                )
                print >>o_file, longtable_template.render(
                    table_caption=table_caption,
                    table_label=table_label,
                    nb_col=len(table[1]),
                    col_alignments="lr",
                    text_size="normalsize",
                    header_data=zip(table[0], [1 for i in table[0]]),
                    tabular_data=table[1:],
                )

    except IOError:
        msg = "{}: cannot write LaTeX summary".format(latex_file)
        raise ProgramError(msg)

    # Writing the summary results
    with open(os.path.join(base_dir, "results_summary.txt"), "a") as o_file:
        print >>o_file, "# {}".format(script_prefix)
        print >>o_file, ("Number of markers with plate bias (p<{})\t"
                         "{:,d}".format(p_threshold, nb_markers))
        print >>o_file, "---"

    # We know this step doesn't produce an new data set, so we return the old
    # prefix and the old in_type
    return _StepResult(
        next_file=in_prefix,
        next_file_type=required_type,
        latex_summary=latex_file,
        description=plate_bias.desc,
        long_description=plate_bias.long_desc,
        graph_path=None,
    )


def run_remove_heterozygous_haploid(in_prefix, in_type, out_prefix, base_dir,
                                    options):
    """Runs step8 (remove heterozygous haploid).

    :param in_prefix: the prefix of the input files.
    :param in_type: the type of the input files.
    :param out_prefix: the output prefix.
    :param base_dir: the output directory.
    :param options: the options needed.

    :type in_prefix: str
    :type in_type: str
    :type out_prefix: str
    :type base_dir: str
    :type options: list

    :returns: a tuple containing the prefix of the output files (the input
              prefix for the next script) and the type of the output files
              (``bfile``).

    This function calls the
    :py:mod:`pyGenClean.HeteroHap.remove_heterozygous_haploid` module. The
    required file type for this module is ``bfile``, hence the need to use the
    :py:func:`check_input_files` to check if the file input file type is the
    good one, or to create it if needed.

    """
    # Creating the output directory
    os.mkdir(out_prefix)

    # We know we need bfile
    required_type = "bfile"
    check_input_files(in_prefix, in_type, required_type)

    # We need to inject the name of the input file and the name of the output
    # prefix
    script_prefix = os.path.join(out_prefix, "without_hh_genotypes")
    options += ["--{}".format(required_type), in_prefix,
                "--out", script_prefix]

    # We run the script
    try:
        remove_heterozygous_haploid.main(options)
    except remove_heterozygous_haploid.ProgramError as e:
        msg = "remove_heterozygous_haploid: {}".format(e)
        raise ProgramError(msg)

    # We get the number of genotypes that were set to missing
    nb_hh_missing = None
    with open(script_prefix + ".log", "r") as i_file:
        nb_hh_missing = re.search(
            r"(\d+) heterozygous haploid genotypes; set to missing",
            i_file.read(),
        )
    if nb_hh_missing:
        nb_hh_missing = int(nb_hh_missing.group(1))
    else:
        nb_hh_missing = 0

    # We write a LaTeX summary
    latex_file = os.path.join(script_prefix + ".summary.tex")
    try:
        with open(latex_file, "w") as o_file:
            print >>o_file, latex_template.subsection(
                remove_heterozygous_haploid.pretty_name
            )
            text = (
                "After Plink's heterozygous haploid analysis, a total of "
                "{:,d} genotype{} were set to missing.".format(
                    nb_hh_missing,
                    "s" if nb_hh_missing > 1 else "",
                )
            )
            print >>o_file, latex_template.wrap_lines(text)

    except IOError:
        msg = "{}: cannot write LaTeX summary".format(latex_file)
        raise ProgramError(msg)

    # Writing the summary results
    with open(os.path.join(base_dir, "results_summary.txt"), "a") as o_file:
        print >>o_file, "# {}".format(script_prefix)
        print >>o_file, ("Number of heterozygous haploid genotypes set to "
                         "missing\t{:,d}".format(nb_hh_missing))
        print >>o_file, "---"

    # We know this step produces an new data set (bfile), so we return it
    return _StepResult(
        next_file=os.path.join(out_prefix, "without_hh_genotypes"),
        next_file_type="bfile",
        latex_summary=latex_file,
        description=remove_heterozygous_haploid.desc,
        long_description=remove_heterozygous_haploid.long_desc,
        graph_path=None,
    )


def run_find_related_samples(in_prefix, in_type, out_prefix, base_dir,
                             options):
    """Runs step9 (find related samples).

    :param in_prefix: the prefix of the input files.
    :param in_type: the type of the input files.
    :param out_prefix: the output prefix.
    :param base_dir: the output directory.
    :param options: the options needed.

    :type in_prefix: str
    :type in_type: str
    :type out_prefix: str
    :type base_dir: str
    :type options: list

    :returns: a tuple containing the prefix of the output files (the input
              prefix for the next script) and the type of the output files
              (``bfile``).

    This function calls the
    :py:mod:`pyGenClean.RelatedSamples.find_related_samples` module. The
    required file type for this module is ``bfile``, hence the need to use the
    :py:func:`check_input_files` to check if the file input file type is the
    good one, or to create it if needed.

    .. note::
        The :py:mod:`pyGenClean.RelatedSamples.find_related_samples` module
        doesn't return usable output files. Hence, this function returns the
        input file prefix and its type.

    """
    # Creating the output directory
    os.mkdir(out_prefix)

    # We know we need bfile
    required_type = "bfile"
    check_input_files(in_prefix, in_type, required_type)

    # We need to inject the name of the input file and the name of the output
    # prefix
    script_prefix = os.path.join(out_prefix, "ibs")
    options += ["--{}".format(required_type), in_prefix,
                "--out", script_prefix]

    # The IBS2 ratio
    ibs2_ratio = find_related_samples._ibs2_ratio_default
    if "--ibs2-ratio" in options:
        ibs2_ratio = options[options.index("--ibs2-ratio") + 1]

    # The indep pairwiase
    r2_value = find_related_samples._indep_pairwise_r2_default
    if "--indep-pairwise" in options:
        r2_value = options[options.index("--indep-pairwise") + 3]

    # We run the script
    try:
        find_related_samples.main(options)
    except find_related_samples.ProgramError as e:
        msg = "find_related_samples: {}".format(e)
        raise ProgramError(msg)

    # Reading the file containing all samples that are related
    #   - ibs.related_individuals
    related_samples = set()
    with open(script_prefix + ".related_individuals", "r") as i_file:
        header = {
            name: i for i, name in
            enumerate(createRowFromPlinkSpacedOutput(i_file.readline()))
        }
        for required_col in ["FID1", "IID1", "FID2", "IID2"]:
            if required_col not in header:
                msg = "{}: missing column {}".format(
                    script_prefix + ".related_individuals",
                    required_col,
                )
                raise ProgramError(msg)

        # Reading the rest of the data
        for line in i_file:
            row = createRowFromPlinkSpacedOutput(line)
            related_samples.add((row[header["FID1"]], row[header["IID1"]]))
            related_samples.add((row[header["FID2"]], row[header["IID2"]]))

    # Reading file containing samples that should be discarded
    #   - ibs.discarded_related_individuals
    discarded_samples = None
    with open(script_prefix + ".discarded_related_individuals", "r") as i_file:
        discarded_samples = {
            tuple(i.rstrip("\r\n").split("\t")) for i in i_file
        }

    # Counting the number of markers used for computing IBS values
    nb_markers_ibs = 0
    with open(script_prefix + ".pruned_data.bim", "r") as i_file:
        for line in i_file:
            nb_markers_ibs += 1

    # Reading the merged related individuals file
    table = []
    with open(script_prefix + ".merged_related_individuals", "r") as i_file:
        for line in i_file:
            table.append([
                latex_template.sanitize_tex(item)
                for item in line.rstrip("\r\n").split("\t")
            ])

    # We write a LaTeX summary
    latex_file = os.path.join(script_prefix + ".summary.tex")
    graphics_paths = set()
    try:
        with open(latex_file, "w") as o_file:
            print >>o_file, latex_template.subsection(
                find_related_samples.pretty_name,
            )
            text = (
                "According to Plink relatedness analysis (using {:,d} "
                "marker{}), {:,d} unique sample{} {} related to at least one "
                "other sample. A total of {:,d} sample{} {} randomly selected "
                "for downstream exclusion from the dataset.".format(
                    nb_markers_ibs,
                    "s" if nb_markers_ibs > 1 else "",
                    len(related_samples),
                    "s" if len(related_samples) > 1 else "",
                    "were" if len(related_samples) > 1 else "was",
                    len(discarded_samples),
                    "s" if len(discarded_samples) > 1 else "",
                    "were" if len(discarded_samples) > 1 else "was",
                )
            )
            print >>o_file, latex_template.wrap_lines(text)

            # Adding the first figure (if present)
            fig_1 = script_prefix + ".related_individuals_z1.png"
            fig_1_label = re.sub(r"[/\\]", "_", script_prefix) + "_z1"
            if os.path.isfile(fig_1):
                text = (
                    r"Figure~\ref{" + fig_1_label + "} shows $Z_1$ versus "
                    r"$IBS2_{ratio}^\ast$ for all related samples found by "
                    r"Plink."
                )
                print >>o_file, latex_template.wrap_lines(text)

            # Adding the second figure (if present)
            fig_2 = script_prefix + ".related_individuals_z2.png"
            fig_2_label = re.sub(r"[/\\]", "_", script_prefix) + "_z2"
            if os.path.isfile(fig_2):
                text = (
                    r"Figure~\ref{" + fig_2_label + "} shows $Z_2$ versus "
                    r"$IBS2_{ratio}^\ast$ for all related samples found by "
                    r"Plink."
                )
                print >>o_file, latex_template.wrap_lines(text)

            if len(table) > 1:
                # There is data, so we add
                table_label = re.sub(r"[/\\]", "_", script_prefix) + "_related"
                text = (
                    r"Table~\ref{" + table_label + "} lists the related "
                    "sample pairs with estimated relationship."
                )
                print >>o_file, latex_template.wrap_lines(text)

            # Getting the required template
            float_template = latex_template.jinja2_env.get_template(
                "float_template.tex",
            )
            graphic_template = latex_template.jinja2_env.get_template(
                "graphics_template.tex",
            )

            figures = (fig_1, fig_2)
            labels = (fig_1_label, fig_2_label)
            graph_types = ("$Z_1$", "$Z_2$")
            for fig, label, graph_type in zip(figures, labels, graph_types):
                if os.path.isfile(fig):
                    # Getting the paths
                    graphics_path, path = os.path.split(fig)
                    graphics_path = os.path.relpath(graphics_path, base_dir)

                    # Printing
                    caption = (
                        graph_type + r" versus $IBS2_{ratio}^\ast$ for all "
                        "related samples found by Plink in the IBS analysis."
                    )
                    print >>o_file, float_template.render(
                        float_type="figure",
                        float_placement="H",
                        float_caption=caption,
                        float_label=label,
                        float_content=graphic_template.render(
                            width=r"0.8\textwidth",
                            path=latex_template.sanitize_fig_name(path),
                        ),
                    )
                    # Adding the path where the graphic is
                    graphics_paths.add(graphics_path)

            # Adding the table
            if len(table) > 1:
                # Getting the template
                longtable_template = latex_template.jinja2_env.get_template(
                    "longtable_template.tex",
                )

                # The table caption
                table_caption = (
                    "List of all related samples with estimated relationship. "
                    "Sample pairs are grouped according to their estimated "
                    "family (the index column)."
                )
                print >>o_file, longtable_template.render(
                    table_caption=table_caption,
                    table_label=table_label,
                    nb_col=len(table[1]),
                    col_alignments="rlllll",
                    text_size="scriptsize",
                    header_data=zip(table[0], [1 for i in table[0]]),
                    tabular_data=table[1:],
                )

    except IOError:
        msg = "{}: cannot write LaTeX summary".format(latex_file)
        raise ProgramError(msg)

    # Writing the summary results
    with open(os.path.join(base_dir, "results_summary.txt"), "a") as o_file:
        print >>o_file, "# {}".format(script_prefix)
        print >>o_file, ("Number of markers used for IBS analysis\t"
                         "{:,d}".format(nb_markers_ibs))
        print >>o_file, ("Number of unique related samples\t"
                         "{:,d}".format(len(related_samples)))
        print >>o_file, "---"

    # The long description
    long_description = find_related_samples.long_desc.format(
        ratio="{ratio}",
        ratio_value=ibs2_ratio,
        r_squared=r2_value,
    )

    # We know this step doesn't produce an new data set, so we return the old
    # prefix and the old in_type
    return _StepResult(
        next_file=in_prefix,
        next_file_type=required_type,
        latex_summary=latex_file,
        description=find_related_samples.desc,
        long_description=long_description,
        graph_path=graphics_paths,
    )


def run_check_ethnicity(in_prefix, in_type, out_prefix, base_dir, options):
    """Runs step10 (check ethnicity).

    :param in_prefix: the prefix of the input files.
    :param in_type: the type of the input files.
    :param out_prefix: the output prefix.
    :param base_dir: the output directory.
    :param options: the options needed.

    :type in_prefix: str
    :type in_type: str
    :type out_prefix: str
    :type base_dir: str
    :type options: list

    :returns: a tuple containing the prefix of the output files (the input
              prefix for the next script) and the type of the output files
              (``bfile``).

    This function calls the :py:mod:`pyGenClean.Ethnicity.check_ethnicity`
    module. The required file type for this module is ``bfile``, hence the need
    to use the :py:func:`check_input_files` to check if the file input file
    type is the good one, or to create it if needed.

    .. note::
        The :py:mod:`pyGenClean.Ethnicity.check_ethnicity` module doesn't
        return usable output files. Hence, this function returns the input file
        prefix and its type.

    """
    # Creating the output directory
    os.mkdir(out_prefix)

    # We know we need bfile
    required_type = "bfile"
    check_input_files(in_prefix, in_type, required_type)

    # We need to inject the name of the input file and the name of the output
    # prefix
    script_prefix = os.path.join(out_prefix, "ethnicity")
    options += ["--{}".format(required_type), in_prefix,
                "--out", script_prefix]

    # We run the script
    try:
        check_ethnicity.main(options)
    except check_ethnicity.ProgramError as e:
        msg = "check_ethnicity: {}".format(e)
        raise ProgramError(msg)

    # Getting the multiplier value
    multiplier = check_ethnicity.parser.get_default("multiplier")
    if "--multiplier" in options:
        multiplier = options[options.index("--multiplier") + 1]

    # Getting the population of which the outliers were computed from
    outliers_of = check_ethnicity.parser.get_default("outliers_of")
    if "--outliers-of" in options:
        outliers_of = options[options.index("--outliers-of") + 1]

    # Was the reference populations required?
    skip_ref_pops = "--skip-ref-pops" in options

    # Computing the number of outliers
    outliers = None
    if not skip_ref_pops:
        with open(script_prefix + ".outliers", "r") as i_file:
            outliers = {
                tuple(line.rstrip("\r\n").split("\t")) for line in i_file
            }

    # Computing the number of markers used
    nb_markers_mds = 0
    with open(script_prefix + ".ibs.pruned_data.bim", "r") as i_file:
        for line in i_file:
            nb_markers_mds += 1

    # We write a LaTeX summary
    latex_file = os.path.join(script_prefix + ".summary.tex")
    graphics_paths = set()
    try:
        # TODO: IF THIS CHANGE, code in find_outliers needs to change to...
        with open(latex_file, "w") as o_file:
            print >>o_file, latex_template.subsection(
                check_ethnicity.pretty_name,
            )
            text = None
            if skip_ref_pops:
                text = (
                    "Principal components analysis was performed using {:,d} "
                    "marker{} on the study dataset only.".format(
                        nb_markers_mds,
                        "s" if nb_markers_mds > 1 else "",
                    )
                )

            else:
                text = (
                    "Using {:,d} marker{} and a multiplier of {}, there was a "
                    "total of {:,d} outlier{} of the {} population.".format(
                        nb_markers_mds,
                        "s" if nb_markers_mds > 1 else "",
                        multiplier,
                        len(outliers),
                        "s" if len(outliers) > 1 else "",
                        outliers_of,
                    )
                )
            print >>o_file, latex_template.wrap_lines(text)

            # Adding the figure if it exists
            fig = script_prefix + ".outliers.png"
            if os.path.isfile(fig):
                # Getting the paths
                graphics_path, path = os.path.split(fig)
                graphics_path = os.path.relpath(graphics_path, base_dir)

                # Getting the required template
                float_template = latex_template.jinja2_env.get_template(
                    "float_template.tex",
                )
                graphic_template = latex_template.jinja2_env.get_template(
                    "graphics_template.tex",
                )

                # The label
                label = re.sub(r"[/\\]", "_", script_prefix) + "_outliers"

                text = (
                    r"Figure~\ref{" + label + "} shows the first two "
                    "principal components of the MDS analysis, where outliers "
                    "of the " + outliers_of + " population are shown in grey."
                )
                print >>o_file, latex_template.wrap_lines(text)

                # Printing
                caption = (
                    "MDS plots showing the first two principal components of "
                    "the source dataset with the reference panels. The "
                    "outliers of the {} population are shown in grey, while "
                    "samples of the source dataset that resemble the {} "
                    "population are shown in orange. A multiplier of {} was "
                    "used to find the {:,d} outlier{}.".format(
                        outliers_of,
                        outliers_of,
                        multiplier,
                        len(outliers),
                        "s" if len(outliers) > 1 else "",
                    )
                )
                print >>o_file, float_template.render(
                    float_type="figure",
                    float_placement="H",
                    float_caption=latex_template.wrap_lines(caption),
                    float_label=label,
                    float_content=graphic_template.render(
                        width=r"0.8\textwidth",
                        path=latex_template.sanitize_fig_name(path),
                    ),
                )
                # Adding the path where the graphic is
                graphics_paths.add(graphics_path)

            # Adding the screeplot if it exists
            fig = script_prefix + ".smartpca.scree_plot.png"
            if os.path.isfile(fig):
                # Getting the paths
                graphics_path, path = os.path.split(fig)
                graphics_path = os.path.relpath(graphics_path, base_dir)

                # Getting the required template
                float_template = latex_template.jinja2_env.get_template(
                    "float_template.tex",
                )
                graphic_template = latex_template.jinja2_env.get_template(
                    "graphics_template.tex",
                )

                # The label
                label = re.sub(r"[/\\]", "_", script_prefix) + "_screeplot"

                text = (
                    r"Figure~\ref{" + label + "} shows the scree plot for the "
                    "principal components of the MDS analysis."
                )
                print >>o_file, latex_template.wrap_lines(text)

                # Printing
                caption = (
                    "Scree plot for the principal components of the MDS "
                    "analysis."
                )
                print >>o_file, float_template.render(
                    float_type="figure",
                    float_placement="H",
                    float_caption=caption,
                    float_label=label,
                    float_content=graphic_template.render(
                        height=r"0.95\textheight",
                        path=latex_template.sanitize_fig_name(path),
                    ),
                )
                # Adding the path where the graphic is
                graphics_paths.add(graphics_path)

    except IOError:
        msg = "{}: cannot write LaTeX summary".format(latex_file)
        raise ProgramError(msg)

    # Writing the summary results
    with open(os.path.join(base_dir, "results_summary.txt"), "a") as o_file:
        print >>o_file, "# {}".format(script_prefix)
        print >>o_file, ("Number of markers used for MDS analysis\t"
                         "{:,d}".format(nb_markers_mds))
        if not skip_ref_pops:
            print >>o_file, ("Number of {} outliers\t"
                             "{:,d}".format(outliers_of, len(outliers)))
        print >>o_file, "---"

    # We know this step doesn't produce an new data set, so we return the old
    # prefix and the old in_type
    return _StepResult(
        next_file=in_prefix,
        next_file_type=required_type,
        latex_summary=latex_file,
        description=check_ethnicity.desc,
        long_description=check_ethnicity.long_desc,
        graph_path=graphics_paths,
    )


def run_flag_maf_zero(in_prefix, in_type, out_prefix, base_dir, options):
    """Runs step11 (flag MAF zero).

    :param in_prefix: the prefix of the input files.
    :param in_type: the type of the input files.
    :param out_prefix: the output prefix.
    :param base_dir: the output directory.
    :param options: the options needed.

    :type in_prefix: str
    :type in_type: str
    :type out_prefix: str
    :type base_dir: str
    :type options: list

    :returns: a tuple containing the prefix of the output files (the input
              prefix for the next script) and the type of the output files
              (``bfile``).

    This function calls the :py:mod:`pyGenClean.FlagMAF.flag_maf_zero` module.
    The required file type for this module is ``bfile``, hence the need to use
    the :py:func:`check_input_files` to check if the file input file type is
    the good one, or to create it if needed.

    .. note::
        The :py:mod:`pyGenClean.FlagMAF.flag_maf_zero` module doesn't return
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
    script_prefix = os.path.join(out_prefix, "flag_maf_0")
    options += ["--{}".format(required_type), in_prefix,
                "--out", script_prefix]

    # We run the script
    try:
        flag_maf_zero.main(options)
    except flag_maf_zero.ProgramError as e:
        msg = "flag_maf_zero: {}".format(e)
        raise ProgramError(msg)

    # Reading the file to compute the number of flagged markers
    nb_flagged = None
    flagged_fn = script_prefix + ".list"
    with open(flagged_fn, "r") as i_file:
        nb_flagged = len(i_file.read().splitlines())

    # We write a LaTeX summary
    latex_file = os.path.join(script_prefix + ".summary.tex")
    try:
        with open(latex_file, "w") as o_file:
            print >>o_file, latex_template.subsection(
                flag_maf_zero.pretty_name
            )
            safe_fn = latex_template.sanitize_tex(os.path.basename(flagged_fn))
            text = (
                "After computing minor allele frequencies (MAF) of all "
                "markers using Plink, a total of {:,d} marker{} had a MAF "
                "of zero and were flagged ({}).".format(
                    nb_flagged,
                    "s" if nb_flagged - 1 > 1 else "",
                    "see file " + latex_template.texttt(safe_fn) +
                    " for more information"
                )
            )
            print >>o_file, latex_template.wrap_lines(text)

    except IOError:
        msg = "{}: cannot write LaTeX summary".format(latex_file)
        raise ProgramError(msg)

    # Writing the summary results
    with open(os.path.join(base_dir, "results_summary.txt"), "a") as o_file:
        print >>o_file, "# {}".format(script_prefix)
        print >>o_file, ("Number of markers flagged for MAF of 0\t"
                         "{:,d}".format(nb_flagged))
        print >>o_file, "---"

    # We know this step doesn't produce an new data set, so we return the old
    # prefix and the old in_type
    return _StepResult(
        next_file=in_prefix,
        next_file_type=required_type,
        latex_summary=latex_file,
        description=flag_maf_zero.desc,
        long_description=flag_maf_zero.long_desc,
        graph_path=None,
    )


def run_flag_hw(in_prefix, in_type, out_prefix, base_dir, options):
    """Runs step12 (flag HW).

    :param in_prefix: the prefix of the input files.
    :param in_type: the type of the input files.
    :param out_prefix: the output prefix.
    :param base_dir: the output directory.
    :param options: the options needed.

    :type in_prefix: str
    :type in_type: str
    :type out_prefix: str
    :type base_dir: str
    :type options: list

    :returns: a tuple containing the prefix of the output files (the input
              prefix for the next script) and the type of the output files
              (``bfile``).

    This function calls the :py:mod:`pyGenClean.FlagHW.flag_hw` module. The
    required file type for this module is ``bfile``, hence the need to use the
    :py:func:`check_input_files` to check if the file input file type is the
    good one, or to create it if needed.

    .. note::
        The :py:mod:`pyGenClean.FlagHW.flag_hw` module doesn't return usable
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
    script_prefix = os.path.join(out_prefix, "flag_hw")
    options += ["--{}".format(required_type), in_prefix,
                "--out", script_prefix]

    # We run the script
    try:
        flag_hw.main(options)
    except flag_hw.ProgramError as e:
        msg = "flag_hw: {}".format(e)
        raise ProgramError(msg)

    # Finding the two files containing the list of flagged markers
    filenames = glob(script_prefix + ".snp_flag_threshold_[0-9]*")
    thresholds = {}
    for filename in filenames:
        # Finding the threshold of the file
        threshold = re.sub(
            r"^flag_hw.snp_flag_threshold_",
            "",
            os.path.basename(filename),
        )

        # Counting the number of markers in the file
        nb_markers = None
        with open(filename, "r") as i_file:
            nb_markers = len(i_file.read().splitlines())

        # Saving the values
        thresholds[threshold] = (nb_markers, filename)

    # We create the LaTeX summary
    latex_file = os.path.join(script_prefix + ".summary.tex")
    try:
        with open(latex_file, "w") as o_file:
            print >>o_file, latex_template.subsection(
                flag_hw.pretty_name
            )

            # Data to write
            sorted_keys = sorted(thresholds.keys(), key=float)

            text = (
                "Markers which failed Hardy-Weinberg equilibrium test (using "
                "Plink) were flagged. A total of {:,d} marker{} failed with a "
                "threshold of {}. A total of {:,d} marker{} failed with a "
                "threshold of {}. For a total list, check the files {} and "
                "{}, respectively.".format(
                    thresholds[sorted_keys[0]][0],
                    "s" if thresholds[sorted_keys[0]][0] - 1 > 1 else "",
                    latex_template.format_numbers(sorted_keys[0]),
                    thresholds[sorted_keys[1]][0],
                    "s" if thresholds[sorted_keys[1]][0] - 1 > 1 else "",
                    latex_template.format_numbers(sorted_keys[1]),
                    latex_template.texttt(
                        latex_template.sanitize_tex(os.path.basename(
                            thresholds[sorted_keys[0]][1],
                        )),
                    ),
                    latex_template.texttt(
                        latex_template.sanitize_tex(os.path.basename(
                            thresholds[sorted_keys[1]][1],
                        )),
                    ),
                )
            )
            print >>o_file, latex_template.wrap_lines(text)

    except IOError:
        msg = "{}: cannot write LaTeX summary".format(latex_file)
        raise ProgramError(msg)

    # Writing the summary results
    with open(os.path.join(base_dir, "results_summary.txt"), "a") as o_file:
        print >>o_file, "# {}".format(script_prefix)
        print >>o_file, "Number of markers flagged for HW"
        print >>o_file, "  - {}\t{:,d}".format(
                            sorted_keys[0],
                            thresholds[sorted_keys[0]][0],
                        )
        print >>o_file, "  - {}\t{:,d}".format(
                            sorted_keys[1],
                            thresholds[sorted_keys[1]][0],
                        )
        print >>o_file, "---"

    # We know this step doesn't produce an new data set, so we return the old
    # prefix and the old in_type
    return _StepResult(
        next_file=in_prefix,
        next_file_type=required_type,
        latex_summary=latex_file,
        description=flag_hw.desc,
        long_description=flag_hw.long_desc,
        graph_path=None,
    )


def run_compare_gold_standard(in_prefix, in_type, out_prefix, base_dir,
                              options):
    """Compares with a gold standard data set (compare_gold_standard.

    :param in_prefix: the prefix of the input files.
    :param in_type: the type of the input files.
    :param out_prefix: the output prefix.
    :param base_dir: the output directory.
    :param options: the options needed.

    :type in_prefix: str
    :type in_type: str
    :type out_prefix: str
    :type base_dir: str
    :type options: list

    :returns: a tuple containing the prefix of the output files (the input
              prefix for the next script) and the type of the output files
              (``bfile``).

    This function calls the :py:mod:`pyGenClean.Misc.compare_gold_standard`
    module. The required file type for this module is ``bfile``, hence the need
    to use the :py:func:`check_input_files` to check if the file input file
    type is the good one, or to create it if needed.

    .. note::
        The :py:mod:`pyGenClean.Misc.compare_gold_standard` module doesn't
        return usable output files. Hence, this function returns the input file
        prefix and its type.

    """
    # Creating the output directory
    os.mkdir(out_prefix)

    # We know we need bfile
    required_type = "bfile"
    check_input_files(in_prefix, in_type, required_type)

    # We need to inject the name of the input file and the name of the output
    # prefix
    script_prefix = os.path.join(out_prefix, "compare_with_gold")
    options += ["--{}".format(required_type), in_prefix,
                "--out", script_prefix]

    # We run the script
    try:
        compare_gold_standard.main(options)
    except compare_gold_standard.ProgramError as e:
        msg = "compare_gold_standard: {}".format(e)
        raise ProgramError(msg)

    # We create the LaTeX summary
    latex_file = os.path.join(script_prefix + ".summary.tex")
    try:
        with open(latex_file, "w") as o_file:
            print >>o_file, latex_template.subsection(
                compare_gold_standard.pretty_name
            )

    except IOError:
        msg = "{}: cannot write LaTeX summary".format(latex_file)
        raise ProgramError(msg)

    # We know this step doesn't produce an new data set, so we return the old
    # prefix and the old in_type
    return _StepResult(
        next_file=in_prefix,
        next_file_type=required_type,
        latex_summary=latex_file,
        description=compare_gold_standard.desc,
        long_description=compare_gold_standard.long_desc,
        graph_path=None,
    )


def run_subset_data(in_prefix, in_type, out_prefix, base_dir, options):
    """Subsets the data.

    :param in_prefix: the prefix of the input files.
    :param in_type: the type of the input files.
    :param out_prefix: the output prefix.
    :param base_dir: the output directory.
    :param options: the options needed.

    :type in_prefix: str
    :type in_type: str
    :type out_prefix: str
    :type base_dir: str
    :type options: list

    :returns: a tuple containing the prefix of the output files (the input
              prefix for the next script) and the type of the output files
              (``bfile``).

    This function calls the
    :py:mod:`pyGenClean.pyGenClean.PlinkUtils.subset_data` module. The required
    file type for this module is ``bfile``, hence the need to use the
    :py:func:`check_input_files` to check if the file input file type is the
    good one, or to create it if needed.

    .. note::
        The output file type is the same as the input file type.

    """
    # Creating the output directory
    os.mkdir(out_prefix)

    # The prefix of the script
    script_prefix = os.path.join(out_prefix, "subset")

    # The extension of the marker and sample file
    markers_ext = None
    samples_ext = None

    # Looking at what we have
    required_type = None
    if in_type == "bfile":
        required_type = "bfile"
        markers_ext = ".bim"
        samples_ext = ".fam"

    elif in_type == "tfile":
        required_type = "tfile"
        markers_ext = ".tped"
        samples_ext = ".tfam"

    else:
        required_type = "file"
        markers_ext = ".map"
        samples_ext = ".ped"

    if "--is-bfile" in set(options):
        required_type = "bfile"

    # Checking the input file
    check_input_files(in_prefix, in_type, required_type)

    # What is going to be done
    is_extract = "--extract" in options
    is_exclude = "--exclude" in options
    is_keep = "--keep" in options
    is_remove = "--remove" in options

    # Checking if we have a reason for markers
    m_subset_reason = None
    if "--reason-marker" in options:
        m_subset_reason = latex_template.sanitize_tex(
            options.pop(options.index("--reason-marker")+1),
        )
        options.pop(options.index("--reason-marker"))

        if (not is_extract) and (not is_exclude):
            # There was a reason for marker exclusion, but no exclusion will be
            # performed, so we skip in the report
            m_subset_reason = None

    # Checking if we have a reason for samples
    s_subset_reason = None
    if "--reason-sample" in options:
        s_subset_reason = latex_template.sanitize_tex(
            options.pop(options.index("--reason-sample")+1),
        )
        options.pop(options.index("--reason-sample"))

        if (not is_keep) and (not is_remove):
            # There was a reason for sample exclusion, but no exclusion will be
            # performed, so we skip in the report
            s_subset_reason = None

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

    # The files before the subset
    samples_before_fn = in_prefix
    markers_before_fn = in_prefix

    # The files after the subset
    samples_after_fn = script_prefix
    markers_after_fn = script_prefix

    # The name of the subset file for markers and samples
    sample_subset_fn = None
    marker_subset_fn = None

    # The samples and markers that were removed
    removed_samples = {}
    removed_markers = {}

    # The set of samples and markers before subset
    samples_before = None
    markers_before = None

    # The set of remaining samples and markers after subset
    samples_after = None
    markers_after = None

    # Markers were extracted
    nb_extract = None
    if is_extract:
        # Counting the number of markers that were extracted
        marker_subset_fn = options[options.index("--extract") + 1]
        with open(marker_subset_fn, "r") as i_file:
            nb_extract = sum(1 for line in i_file)

    # Markers were excluded
    nb_exclude = None
    if is_exclude:
        marker_subset_fn = options[options.index("--exclude") + 1]
        with open(marker_subset_fn, "r") as i_file:
            nb_exclude = sum(1 for line in i_file)

    # Checking the difference between both (for markers)
    if is_extract or is_exclude:
        markers_before_fn += markers_ext
        markers_after_fn += ".bim" if required_type == "bfile" else markers_ext

        markers_before = None
        with open(markers_before_fn, "r") as i_file:
            markers_before = {
                createRowFromPlinkSpacedOutput(line)[1] for line in i_file
            }

        with open(markers_after_fn, "r") as i_file:
            markers_after = {
                createRowFromPlinkSpacedOutput(line)[1] for line in i_file
            }

        reason = os.path.relpath(marker_subset_fn, base_dir)
        if m_subset_reason is not None:
            reason = m_subset_reason
        removed_markers = {
            name: "subset {}".format(reason)
            for name in markers_before - markers_after
        }

    # Samples were kept
    nb_keep = None
    if is_keep:
        sample_subset_fn = options[options.index("--keep") + 1]
        with open(sample_subset_fn, "r") as i_file:
            nb_keep = sum(1 for line in i_file)

    # Samples were removed
    nb_remove = None
    if is_remove:
        sample_subset_fn = options[options.index("--remove") + 1]
        with open(sample_subset_fn, "r") as i_file:
            nb_remove = sum(1 for line in i_file)

    if is_keep or is_remove:
        samples_before_fn += samples_ext
        samples_after_fn += ".fam" if required_type == "bfile" else samples_ext

        samples_before = None
        with open(samples_before_fn, "r") as i_file:
            samples_before = {
                tuple(createRowFromPlinkSpacedOutput(line)[:2])
                for line in i_file
            }

        with open(samples_after_fn, "r") as i_file:
            samples_after = {
                tuple(createRowFromPlinkSpacedOutput(line)[:2])
                for line in i_file
            }

        reason = os.path.relpath(sample_subset_fn, base_dir)
        if s_subset_reason is not None:
            reason = s_subset_reason
        removed_samples = {
            "\t".join(name): "subset {}".format(reason)
            for name in samples_before - samples_after
        }

    # Reading the log file to gather what is left
    log_file = None
    with open(script_prefix + ".log", "r") as i_file:
        log_file = i_file.read()

    # Checking the number of markers at the beginning
    nb_marker_start = re.search(
        r"(\d+) (\(of \d+\) )?markers to be included "
        "from \[ {}".format(in_prefix),
        log_file,
    )
    if nb_marker_start:
        nb_marker_start = int(nb_marker_start.group(1))

    # Checking the number of markers at the end
    nb_marker_end = re.search(
        r"Before frequency and genotyping pruning, there are (\d+) SNPs",
        log_file,
    )
    if nb_marker_end:
        nb_marker_end = int(nb_marker_end.group(1))

    # Checking the number of samples
    if (markers_after is not None) and (nb_marker_end != len(markers_after)):
        raise ProgramError("Something went wrong with Plink's subset (numbers "
                           "are different from data and log file)")
    if nb_marker_start - nb_marker_end != len(removed_markers):
        raise ProgramError("Something went wrong with Plink's subset (numbers "
                           "are different from data and log file)")

    # Checking the number of samples at the beginning
    nb_sample_start = re.search(
        r"(\d+) individuals read from \[ {}".format(in_prefix),
        log_file,
    )
    if nb_sample_start:
        nb_sample_start = int(nb_sample_start.group(1))

    # Checking the number of samples at the end
    nb_sample_end = re.search(
        r"(\d+) founders and (\d+) non-founders found",
        log_file,
    )
    if nb_sample_end:
        nb_sample_end = int(nb_sample_end.group(1)) + \
                        int(nb_sample_end.group(2))

    # Checking the number of samples
    if (samples_after is not None) and (nb_sample_end != len(samples_after)):
        raise ProgramError("Something went wrong with Plink's subset (numbers "
                           "are different from data and log file)")
    if nb_sample_start - nb_sample_end != len(removed_samples):
        raise ProgramError("Something went wrong with Plink's subset (numbers "
                           "are different from data and log file)")

    # Creating the reasons
    reasons = []
    if m_subset_reason is not None:
        reasons.append(m_subset_reason)
    if s_subset_reason is not None:
        reasons.append(s_subset_reason)

    # We write a LaTeX summary
    latex_file = os.path.join(script_prefix + ".summary.tex")
    try:
        with open(latex_file, "w") as o_file:
            section_name = subset_data.pretty_name
            if len(reasons) > 0:
                section_name += " ({})".format(", ".join(reasons))
            print >>o_file, latex_template.subsection(section_name)
            text = ""
            if is_extract:
                text += (
                    "The file for marker extraction contained {:,d} marker{}. "
                    "Out of a total of {:,d} marker{}, {:,d} "
                    "remained ({:,d} excluded).".format(
                        nb_extract,
                        "s" if nb_extract > 1 else "",
                        nb_marker_start,
                        "s" if nb_marker_start > 1 else "",
                        nb_marker_end,
                        nb_marker_start - nb_marker_end,
                    )
                )
            if is_exclude:
                text += (
                    "The file for marker exclusion contained {:,d} marker{}. "
                    "Out of a total of {:,d} marker{}, {:,d} "
                    "remained ({:,d} excluded). ".format(
                        nb_exclude,
                        "s" if nb_exclude > 1 else "",
                        nb_marker_start,
                        "s" if nb_marker_start > 1 else "",
                        nb_marker_end,
                        nb_marker_start - nb_marker_end,
                    )
                )
            if is_keep:
                text += (
                    "The file containing samples to keep contained {:,d} "
                    "sample{}. Out of a total of {:,d} sample{}, {:,d} "
                    "remained ({:,d} removed). ".format(
                        nb_keep,
                        "s" if nb_keep > 1 else "",
                        nb_sample_start,
                        "s" if nb_sample_start > 1 else "",
                        nb_sample_end,
                        nb_sample_start - nb_sample_end,
                    )
                )
            if is_remove:
                text += (
                    "The file containing samples to remove contained {:,d} "
                    "sample{}. Out of a total of {:,d} sample{}, {:,d} "
                    "remained ({:,d} removed). ".format(
                        nb_remove,
                        "s" if nb_remove > 1 else "",
                        nb_sample_start,
                        "s" if nb_sample_start > 1 else "",
                        nb_sample_end,
                        nb_sample_start - nb_sample_end,
                    )
                )
            print >>o_file, latex_template.wrap_lines(text)

    except IOError:
        msg = "{}: cannot write LaTeX summary".format(latex_file)
        raise ProgramError(msg)

    # Writing the excluded samples
    o_filename = os.path.join(base_dir, "excluded_samples.txt")
    with open(o_filename, "a") as o_file:
        for value in removed_samples.items():
            print >>o_file, "\t".join(value)

    # Writing the excluded markers to file
    o_filename = os.path.join(base_dir, "excluded_markers.txt")
    with open(o_filename, "a") as o_file:
        for value in removed_markers.items():
            print >>o_file, "\t".join(value)

    # Writing the summary results
    with open(os.path.join(base_dir, "results_summary.txt"), "a") as o_file:
        print >>o_file, "# {}".format(script_prefix)
        if nb_marker_end != nb_marker_start:
            reason = "_file_path:" + marker_subset_fn
            if m_subset_reason is not None:
                reason = m_subset_reason

            print >>o_file, "Number of markers excluded"
            print >>o_file, ("  - {reason}\t{nb:,d}\t-{nb:,d}".format(
                reason=reason,
                nb=nb_marker_start - nb_marker_end,
            ))
            print >>o_file, "---"
        if nb_sample_end != nb_sample_start:
            reason = "_file_path:" + sample_subset_fn
            if s_subset_reason is not None:
                reason = s_subset_reason

            print >>o_file, "Number of samples removed"
            print >>o_file, ("  - {reason}\t{nb:,d}\t\t-{nb:,d}".format(
                reason=reason,
                nb=nb_sample_start - nb_sample_end,
            ))
            print >>o_file, "---"

    # Modifying the description
    description = subset_data.desc
    if len(reasons) > 0:
        description += " ({})".format(", ".join(reasons))

    # We know this step does produce a new data set (bfile), so we return it
    return _StepResult(
        next_file=os.path.join(out_prefix, "subset"),
        next_file_type=required_type,
        latex_summary=latex_file,
        description=description,
        long_description=subset_data.long_desc,
        graph_path=None,
    )


def run_command(command):
    """Run a command using subprocesses.

    :param command: the command to run.

    :type command: list

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


def count_markers_samples(prefix, file_type):
    """Counts the number of markers and samples in plink file.

    :param prefix: the prefix of the files.
    :param file_type: the file type.

    :type prefix: str
    :type file_type: str

    :returns: the number of markers and samples (in a tuple).

    """
    # The files that will need counting
    sample_file = None
    marker_file = None

    if file_type == "bfile":
        # Binary files (.bed, .bim and .fam)
        sample_file = prefix + ".fam"
        marker_file = prefix + ".bim"

    elif file_type == "file":
        # Pedfile (.ped and .map)
        sample_file = prefix + ".ped"
        marker_file = prefix + ".map"

    elif file_type == "tfile":
        # Transposed pedfile (.tped and .tfam)
        sample_file = prefix + ".tfam"
        marker_file = prefix + ".tped"

    # Counting (this may take some time)
    nb_samples = 0
    with open(sample_file, "r") as f:
        for line in f:
            nb_samples += 1

    nb_markers = 0
    with open(marker_file, "r") as f:
        for line in f:
            nb_markers += 1

    return nb_markers, nb_samples


def check_input_files(prefix, the_type, required_type):
    """Check that the file is of a certain file type.

    :param prefix: the prefix of the input files.
    :param the_type: the type of the input files (bfile, tfile or file).
    :param required_type: the required type of the input files (bfile, tfile or
                          file).

    :type prefix: str
    :type the_type: str
    :type required_type: str

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
        logger.info("Converting {} from {} to {}".format(
            prefix,
            the_type,
            required_type,
        ))
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
        logger.info("Converting {} from {} to {}".format(
            prefix,
            the_type,
            required_type,
        ))
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
        logger.info("Converting {} from {} to {}".format(
            prefix,
            the_type,
            required_type,
        ))
        run_command(plink_command)

        # Everything is now fine
        return True

    else:
        msg = "{}: unknown file format".format(required_type)
        raise ProgramError(msg)


def all_files_exist(file_list):
    """Check if all files exist.

    :param file_list: the names of files to check.

    :type file_list: list

    :returns: ``True`` if all files exist, ``False`` otherwise.

    """
    all_exist = True
    for filename in file_list:
        all_exist = all_exist and os.path.isfile(filename)
    return all_exist


def read_config_file(filename):
    """Reads the configuration file.

    :param filename: the name of the file containing the configuration.

    :type filename: str

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
    * ``contamination`` (:py:func:`run_contamination`)
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
            msg = ("{}: section {}: no variable called 'script'".format(
                filename,
                section,
            ))
            raise ProgramError(msg)
        if script_name not in available_modules:
            msg = ("{}: section {}: script {}: invalid script name".format(
                filename,
                section,
                script_name,
            ))
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

    ======================= =======  ==========================================
             Options         Type                      Description
    ======================= =======  ==========================================
    ``--bfile``             String   The input binary file prefix from Plink.
    ``--tfile``             String   The input transposed file prefix from
                                     Plink.
    ``--file``              String   The input file prefix from Plink.
    ``--conf``              String   The parameter file for the data clean up.
    ``--report-author``     String   The current project number.
    ``--report-number``     String   The current project author.
    ``--report-background`` String   Text of file containing the background
                                     section of the report.
    ======================= =======  ==========================================

    .. note::
        No option check is done here (except for the one automatically done by
        :py:mod:`argparse`). Those need to be done elsewhere (see
        :py:func:`checkArgs`).

    """
    return parser.parse_args()


# The parser object
desc = "Runs the data clean up (pyGenClean version {}).".format(__version__)
parser = argparse.ArgumentParser(description=desc)
parser.add_argument("-v", "--version", action="version",
                    version="pyGenClean version {}".format(__version__))

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

# The options
group = parser.add_argument_group("Report Options")
report_options(group)

group = parser.add_argument_group("Configuration File")
group.add_argument("--conf", type=str, metavar="FILE", required=True,
                   help="The parameter file for the data clean up.")


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
    "contamination": contamination,
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
    "contamination": run_contamination,
    "flag_maf_zero": run_flag_maf_zero,
    "flag_hw": run_flag_hw,
    "subset": run_subset_data,
    "compare_gold_standard": run_compare_gold_standard,
}


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
