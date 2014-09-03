#!/usr/bin/env python

# How to build source distribution (might add --plat-name win32)
# python setup.py sdist --format bztar
# python setup.py sdist --format gztar
# python setup.py sdist --format zip
# python setup.py bdist --format msi

import os
import sys
import glob
import shutil
from setuptools import setup

# Check the Python version
major, minor, micro, s, tmp = sys.version_info
if major==2 and minor<7 or major<2:
    raise SystemExit("""pyGenClean requires Python 2.7 or later.""")
if major==3:
    raise SystemExit("""pyGenClean doesn't work on Python 3...""")

# Creating the "scripts" directory
if os.path.isdir("scripts"):
    shutil.rmtree("scripts")
os.mkdir("scripts")

# Automatically copying the scripts
shutil.copyfile(os.path.join("pyGenClean", "run_data_clean_up.py"),
                os.path.join("scripts", "run_pyGenClean"))
shutil.copyfile(os.path.join("pyGenClean", "DupSamples", "duplicated_samples.py"),
                os.path.join("scripts", "pyGenClean_duplicated_samples"))
shutil.copyfile(os.path.join("pyGenClean", "DupSNPs", "duplicated_snps.py"),
                os.path.join("scripts", "pyGenClean_duplicated_snps"))
shutil.copyfile(os.path.join("pyGenClean", "NoCallHetero", "clean_noCall_hetero_snps.py"),
                os.path.join("scripts", "pyGenClean_clean_noCall_hetero_snps"))
shutil.copyfile(os.path.join("pyGenClean", "NoCallHetero", "heterozygosity_plot.py"),
                os.path.join("scripts", "pyGenClean_heterozygosity_plot"))
shutil.copyfile(os.path.join("pyGenClean", "SampleMissingness", "sample_missingness.py"),
                os.path.join("scripts", "pyGenClean_sample_missingness"))
shutil.copyfile(os.path.join("pyGenClean", "MarkerMissingness", "snp_missingness.py"),
                os.path.join("scripts", "pyGenClean_snp_missingness"))
shutil.copyfile(os.path.join("pyGenClean", "SexCheck", "sex_check.py"),
                os.path.join("scripts", "pyGenClean_sex_check"))
shutil.copyfile(os.path.join("pyGenClean", "SexCheck", "gender_plot.py"),
                os.path.join("scripts", "pyGenClean_gender_plot"))
shutil.copyfile(os.path.join("pyGenClean", "SexCheck", "baf_lrr_plot.py"),
                os.path.join("scripts", "pyGenClean_baf_lrr_plot"))
shutil.copyfile(os.path.join("pyGenClean", "PlateBias", "plate_bias.py"),
                os.path.join("scripts", "pyGenClean_plate_bias"))
shutil.copyfile(os.path.join("pyGenClean", "HeteroHap", "remove_heterozygous_haploid.py"),
                os.path.join("scripts", "pyGenClean_remove_heterozygous_haploid"))
shutil.copyfile(os.path.join("pyGenClean", "RelatedSamples", "find_related_samples.py"),
                os.path.join("scripts", "pyGenClean_find_related_samples"))
shutil.copyfile(os.path.join("pyGenClean", "RelatedSamples", "merge_related_samples.py"),
                os.path.join("scripts", "pyGenClean_merge_related_samples"))
shutil.copyfile(os.path.join("pyGenClean", "Ethnicity", "check_ethnicity.py"),
                os.path.join("scripts", "pyGenClean_check_ethnicity"))
shutil.copyfile(os.path.join("pyGenClean", "Ethnicity", "find_outliers.py"),
                os.path.join("scripts", "pyGenClean_find_outliers"))
shutil.copyfile(os.path.join("pyGenClean", "FlagMAF", "flag_maf_zero.py"),
                os.path.join("scripts", "pyGenClean_flag_maf_zero"))
shutil.copyfile(os.path.join("pyGenClean", "FlagHW", "flag_hw.py"),
                os.path.join("scripts", "pyGenClean_flag_hw"))
shutil.copyfile(os.path.join("pyGenClean", "Misc", "compare_gold_standard.py"),
                os.path.join("scripts", "pyGenClean_compare_gold_standard"))
shutil.copyfile(os.path.join("PlinkUtils", "compare_bim.py"),
                os.path.join("scripts", "pyGenClean_compare_bim"))
shutil.copyfile(os.path.join("PlinkUtils", "plot_MDS_standalone.py"),
                os.path.join("scripts", "pyGenClean_plot_MDS"))
shutil.copyfile(os.path.join("PlinkUtils", "subset_data.py"),
                os.path.join("scripts", "pyGenClean_subset_data"))
shutil.copyfile(os.path.join("pyGenClean", "Ethnicity", "plot_eigenvalues.py"),
                os.path.join("scripts", "pyGenClean_plot_eigenvalues"))

# Changing the mod of the files
if not sys.platform.startswith("win"):
    os.chmod("scripts", 0755)
    for script_name in glob.glob(os.path.join("scripts", "*")):
        os.chmod(script_name, 0755)

setup(name="pyGenClean",
      version="2.0b",
      description="Automated data clean up pipeline",
      author="Louis-Philippe Lemieux Perreault",
      author_email="louis-philippe.lemieux.perreault@umontreal.ca",
      url="http://www.statgen.org",
      license="GPL",
      scripts=[os.path.join("scripts", "{}".format(i))
                    for i in os.listdir("scripts")],
      install_requires=["matplotlib >= 1.2.0", "numpy >= 1.6.2",
                        "scipy >= 0.11.0", "scikit-learn >= 0.12.1",
                        "drmaa >= 0.5"],
      packages=["pyGenClean", "PlinkUtils"] +
               ["pyGenClean.{}".format(i)
                            for i in os.listdir("pyGenClean")
                            if os.path.isdir(os.path.join("pyGenClean", i)) and
                               not i.endswith("git")],
      classifiers=['Operating System :: Linux',
                   'Programming Language :: Python',
                   'Programming Language :: Python :: 2.7'])
