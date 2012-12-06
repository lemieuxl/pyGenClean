#!/usr/bin/env python2.7

import os
import sys
import shutil
import datetime
import argparse
import subprocess
import ConfigParser

import StatGenDataCleanUp.Step1.duplicated_samples as duplicated_samples
import StatGenDataCleanUp.Step2.duplicated_snps as duplicated_snps
import StatGenDataCleanUp.Step3.clean_noCall_hetero_snps as noCall_hetero_snps
import StatGenDataCleanUp.Step4.sample_missingness as sample_missingness
import StatGenDataCleanUp.Step5.snp_missingness as snp_missingness
import StatGenDataCleanUp.Step6.sex_check as sex_check
import StatGenDataCleanUp.Step7.plate_bias as plate_bias
import StatGenDataCleanUp.Step8.remove_heterozygous_haploid as remove_heterozygous_haploid
import StatGenDataCleanUp.Step9.remove_IBS as remove_IBS
import StatGenDataCleanUp.Step10.check_ethnicity as check_ethnicity
import StatGenDataCleanUp.Step11.flag_maf_zero as flag_maf_zero
import StatGenDataCleanUp.Step12.flag_hw as flag_hw
import StatGenDataCleanUp.Misc.compare_gold_standard as compare_gold_standard
import PlinkUtils.subset_data as SubsetData

def main():
    # Getting and checking the options
    args = parse_args()
    check_args(args)

    # Reading the configuration file
    order, conf = read_config_file(args.conf)

    # The directory name
    dirname = "data_clean_up."
    dirname += datetime.datetime.today().strftime("%Y-%m-%d_%H.%M.%S")
    if os.path.isdir(dirname):
        answer = "N"
        if not args.overwrite:
            # The directory already exists...
            print >>sys.stderr, ("WARNING: {}: directory already "
                                 "exists".format(dirname))
            print >>sys.stderr, "Overwrite [Y/N]? ",
            answer = raw_input()
        if args.overwrite or answer.upper() == "Y":
            # Delete everything with the directory
            shutil.rmtree(dirname)
        elif answer.upper() == "N":
            print >>sys.stderr, "STOPING NOW"
            sys.exit(0)
        else:
            msg = "{}: not a valid answer (Y or N)".format(answer)
            raise ProgramError(msg)

    # Creating the output directory
    os.mkdir(dirname)

    # Executing the data clean up
    current_input_file = None
    current_input_type = None
    if args.tfile is not None:
        current_input_file = args.tfile
        current_input_type = "tfile"
    elif args.bfile is not None:
        current_input_file = args.bfile
        current_input_type = "bfile"
    else:
        current_input_file = args.file
        current_input_type = "file"

    for number in order:
        script_name, options = conf[number]
        output_prefix = os.path.join(dirname,
                                     "{}_{}".format(number, script_name))
        function_to_use = available_functions[script_name]
        print "\nRunning {} {}".format(number, script_name)
        print ("   - Using {} as prefix for input "
               "files".format(current_input_file))
        print "   - Results will be in [ {} ]".format(output_prefix)
        current_input_file, current_input_type = function_to_use(current_input_file,
                                                                 current_input_type,
                                                                 output_prefix,
                                                                 options)


def run_duplicated_samples(in_prefix, in_type, out_prefix, options):
    """Runs step1 (duplicated samples)."""
    # Creating the output directory
    os.mkdir(out_prefix)

    # We know we need tfile
    required_type = "tfile"
    check_input_files(in_prefix, in_type, required_type,
                      os.path.basename(out_prefix))

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
    """Runs step2 (duplicated snps)."""
    # Creating the output directory
    os.mkdir(out_prefix)

    # We know we need a tfile
    required_type = "tfile"
    check_input_files(in_prefix, in_type, required_type,
                      os.path.basename(out_prefix))

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


def run_noCall_hetero_snps(in_prefix, in_type, out_prefix, options):
    """Runs step 3 (clean no call and hetero)."""
    # Creating the output directory
    os.mkdir(out_prefix)

    # We know we need a tfile
    required_type = "tfile"
    check_input_files(in_prefix, in_type, required_type,
                      os.path.basename(out_prefix))

    # We need to inject the name of the input file and the name of the output
    # prefix
    options += ["--{}".format(required_type), in_prefix,
                "--out", os.path.join(out_prefix, "clean_noCall_hetero")]

    # We run the script
    try:
        noCall_hetero_snps.main(options)
    except noCall_hetero_snps.ProgramError as e:
        msg = "noCall_hetero_snps: {}".format(e)
        raise ProgramError(msg)

    # We know this step does produce a new data set (tfile), so we return it
    return os.path.join(out_prefix, "clean_noCall_hetero"), "tfile"


def run_sample_missingness(in_prefix, in_type, out_prefix, options):
    """Runs step4 (clean mind)."""
    # Creating the output directory
    os.mkdir(out_prefix)

    # We are looking at what we have
    required_type = "tfile"
    if in_type == "bfile":
        required_type = "bfile"
    check_input_files(in_prefix, in_type, required_type,
                      os.path.basename(out_prefix))

    # We need to inject the name of the input file and the name of the output
    # prefix
    options += ["--ifile", in_prefix,
                "--out", os.path.join(out_prefix, "clean_mind")]
    if required_type == "bfile":
        options.append("--is-bfile")

    # We run the script
    try:
        sample_missingness.main(options)
    except sample_missingness.ProgramError as e:
        msg = "sample_missingness: {}".format(e)
        raise ProgramError(msg)

    # We know this step does produce a new data set (bfile), so we return it
    return os.path.join(out_prefix, "clean_mind"), "bfile"


def run_snp_missingness(in_prefix, in_type, out_prefix, options):
    """Run step5 (clean geno)."""
    # Creating the output directory
    os.mkdir(out_prefix)

    # We know we need a bfile
    required_type = "bfile"
    check_input_files(in_prefix, in_type, required_type,
                      os.path.basename(out_prefix))

    # We need to inject the name of the input file and the name of the output
    # prefix
    options += ["--{}".format(required_type), in_prefix,
                "--out", os.path.join(out_prefix, "clean_geno")]

    # We run the script
    try:
        snp_missingness.main(options)
    except snp_missingness.ProgramError as e:
        msg = "snp_missingness: {}".format(e)
        raise ProgramError(msg)

    # We know this step does produce a new data set (bfile), so we return it
    return os.path.join(out_prefix, "clean_geno"), "bfile"


def run_sex_check(in_prefix, in_type, out_prefix, options):
    """Runs step6 (sexcheck)."""
    # Creating the output directory
    os.mkdir(out_prefix)

    # We know we need a bfile
    required_type = "bfile"
    check_input_files(in_prefix, in_type, required_type,
                      os.path.basename(out_prefix))

    # We need to inject the name of the input file and the name of the output
    # prefix
    options += ["--{}".format(required_type), in_prefix,
                "--out", os.path.join(out_prefix, "sexcheck")]

    # We run the script
    try:
        sex_check.main(options)
    except sex_check.ProgramError as e:
        msg = "sex_check {}".format(e)
        raise ProgramError(msg)

    # We know this step does not produce a new data set, so we return the
    # original one
    return in_prefix, required_type


def run_plate_bias(in_prefix, in_type, out_prefix, options):
    """Runs step7 (plate bias)."""
    # Creating the output directory
    os.mkdir(out_prefix)

    # We know we need bfile
    required_type = "bfile"
    check_input_files(in_prefix, in_type, required_type,
                      os.path.basename(out_prefix))

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
    """Runs step8 (remove heterozygous haploid)."""
    # Creating the output directory
    os.mkdir(out_prefix)

    # We know we need bfile
    required_type = "bfile"
    check_input_files(in_prefix, in_type, required_type,
                      os.path.basename(out_prefix))

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


def run_remove_IBS(in_prefix, in_type, out_prefix, options):
    """Runs step9 (remove IBS)."""
    # Creating the output directory
    os.mkdir(out_prefix)

    # We know we need bfile
    required_type = "bfile"
    check_input_files(in_prefix, in_type, required_type,
                      os.path.basename(out_prefix))

    # We need to inject the name of the input file and the name of the output
    # prefix
    options += ["--{}".format(required_type), in_prefix,
                "--out", os.path.join(out_prefix, "ibs")]

    # We run the script
    try:
        remove_IBS.main(options)
    except remove_IBS.ProgramError as e:
        msg = "remove_IBS: {}".format(e)
        raise ProgramError(msg)

    # We know this step doesn't produce an new data set, so we return the old
    # prefix and the old in_type
    return in_prefix, required_type


def run_check_ethnicity(in_prefix, in_type, out_prefix, options):
    """Runs step10 (check ethnicity)."""
    # Creating the output directory
    os.mkdir(out_prefix)

    # We know we need bfile
    required_type = "bfile"
    check_input_files(in_prefix, in_type, required_type,
                      os.path.basename(out_prefix))

    # We need to inject the name of the input file and the name of the output
    # prefix
    options += ["--{}".format(required_type), in_prefix,
                "--out", os.path.join(out_prefix, "ethnic")]

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
    """Runs step11 (flag MAF zero)."""
    # Creating the output directory
    os.mkdir(out_prefix)

    # We know we need bfile
    required_type = "bfile"
    check_input_files(in_prefix, in_type, required_type,
                      os.path.basename(out_prefix))

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
    """Runs step12 (flag HW)."""
    # Creating the output directory
    os.mkdir(out_prefix)

    # We know we need bfile
    required_type = "bfile"
    check_input_files(in_prefix, in_type, required_type,
                      os.path.basename(out_prefix))

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
    """Compares with a gold standard data set (compare_gold_standard."""
    # Creating the output directory
    os.mkdir(out_prefix)

    # We know we need bfile
    required_type = "bfile"
    check_input_files(in_prefix, in_type, required_type,
                      os.path.basename(out_prefix))

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
    """Subsets the data."""
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
    check_input_files(in_prefix, in_type, required_type,
                      os.path.basename(out_prefix))

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
        SubsetData.main(options)
    except SubsetData.ProgramError as e:
        msg = "subset_data: {}".format(e)
        raise ProgramError(msg)

    # We know this step does produce a new data set (bfile), so we return it
    return os.path.join(out_prefix, "subset"), required_type


def run_command(command):
    """Run a command using subprocesses."""
    output = None
    try:
        output = subprocess.check_output(command, stderr=subprocess.STDOUT,
                                         shell=False)
    except subprocess.CalledProcessError:
        msg = "couldn't run command\n{}".format(command)
        raise ProgramError(msg)


def check_input_files(prefix, the_type, required_type, name):
    """Check that the file is of a certain file type."""
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
    """Check if all files exist."""
    all_exist = True
    for filename in file_list:
        all_exist = all_exist and os.path.isfile(filename)
    return all_exist


def read_config_file(filename):
    """Reads the configuration file."""
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
        if script_name not in available_steps:
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
                    if variable_name == "indep-pairwise":
                        # This is a special option
                        options.extend(variable_value.split(" "))
                    else:
                        options.append(variable_value)

        # Saving the configuration
        configuration[section] = (script_name, options)

    return sections, configuration


def check_args(args):
    """Checks the arguments and options."""
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
    """Parses the command line options and arguments."""
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
prog = "run_data_clean_up"
desc = """Runs the data clean up."""
parser = argparse.ArgumentParser(description=desc, prog=prog)
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

# The maximal step number (from Step1 to StepN)
available_steps = {"duplicated_samples", "duplicated_snps",
                   "noCall_hetero_snps", "sample_missingness",
                   "snp_missingness", "sex_check", "plate_bias",
                   "remove_heterozygous_haploid", "remove_IBS",
                   "check_ethnicity", "flag_maf_zero", "flag_hw", "subset",
                   "compare_gold_standard"}
available_functions = {"duplicated_samples": run_duplicated_samples,
                       "duplicated_snps": run_duplicated_snps,
                       "noCall_hetero_snps": run_noCall_hetero_snps,
                       "sample_missingness": run_sample_missingness,
                       "snp_missingness": run_snp_missingness,
                       "sex_check": run_sex_check,
                       "plate_bias": run_plate_bias,
                       "remove_heterozygous_haploid": run_remove_heterozygous_haploid,
                       "remove_IBS": run_remove_IBS,
                       "check_ethnicity": run_check_ethnicity,
                       "flag_maf_zero": run_flag_maf_zero,
                       "flag_hw": run_flag_hw,
                       "subset": run_subset_data,
                       "compare_gold_standard": run_compare_gold_standard}

# Calling the main, if necessery
if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print >>sys.stderr, "Cancelled by user"
        sys.exit(0)
    except ProgramError as e:
        parser.error(e.message)
