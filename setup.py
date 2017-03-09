#!/usr/bin/env python

# How to build source distribution
#   - python setup.py sdist --format bztar
#   - python setup.py sdist --format gztar
#   - python setup.py sdist --format zip
#   - python setup.py bdist --format msi

# How to build for conda
#   - python setup.py bdist_conda
#   - conda convert -p all /PATH/TO/FILE -o conda_dist
#   - cd conda_dist && conda index *


import os
import sys
from setuptools import setup


MAJOR = 1
MINOR = 8
MICRO = 3
VERSION = "{0}.{1}.{2}".format(MAJOR, MINOR, MICRO)


def check_python_version():
    """Checks the python version, exists if != 2.7."""
    python_major, python_minor = sys.version_info[:2]

    if python_major != 2 or python_minor != 7:
        sys.stderr.write("pyGenClean requires python 2.7")
        sys.exit(1)


def write_version_file(fn=None):
    if fn is None:
        fn = os.path.join(
            os.path.dirname(os.path.abspath(__file__)),
            os.path.join("pyGenClean", "version.py"),
        )

    content = ("\n# THIS FILE WAS GENERATED AUTOMATICALLY BY PYGENCLEAN\n"
               'pygenclean_version = "{version}"\n')

    a = open(fn, "w")
    try:
        a.write(content.format(version=VERSION))
    finally:
        a.close()


def setup_package():
    # Checking the python version
    check_python_version()

    # Saving the version into a file
    write_version_file()

    setup(
        name="pyGenClean",
        version=VERSION,
        description="Automated data clean up pipeline for genetic data",
        long_description=("This package provides tools to automatically "
                          "perform genetic data clean up (QC steps) prior to "
                          "a genome-wide association study."),
        author="Louis-Philippe Lemieux Perreault",
        author_email="louis-philippe.lemieux.perreault@umontreal.ca",
        url="https://github.com/lemieuxl/pyGenClean",
        license="GPL",
        entry_points={
            "console_scripts": [
                "run_pyGenClean=pyGenClean.run_data_clean_up:safe_main",
                "pyGenClean_duplicated_samples=pyGenClean.DupSamples.duplicated_samples:safe_main",
                "pyGenClean_duplicated_snps=pyGenClean.DupSNPs.duplicated_snps:safe_main",
                "pyGenClean_clean_noCall_hetero_snps=pyGenClean.NoCallHetero.clean_noCall_hetero_snps:safe_main",
                "pyGenClean_heterozygosity_plot=pyGenClean.NoCallHetero.heterozygosity_plot:safe_main",
                "pyGenClean_sample_missingness=pyGenClean.SampleMissingness.sample_missingness:safe_main",
                "pyGenClean_snp_missingness=pyGenClean.MarkerMissingness.snp_missingness:safe_main",
                "pyGenClean_sex_check=pyGenClean.SexCheck.sex_check:safe_main",
                "pyGenClean_gender_plot=pyGenClean.SexCheck.gender_plot:safe_main",
                "pyGenClean_baf_lrr_plot=pyGenClean.SexCheck.baf_lrr_plot:safe_main",
                "pyGenClean_plate_bias=pyGenClean.PlateBias.plate_bias:safe_main",
                "pyGenClean_remove_heterozygous_haploid=pyGenClean.HeteroHap.remove_heterozygous_haploid:safe_main",
                "pyGenClean_find_related_samples=pyGenClean.RelatedSamples.find_related_samples:safe_main",
                "pyGenClean_merge_related_samples=pyGenClean.RelatedSamples.merge_related_samples:safe_main",
                "pyGenClean_check_ethnicity=pyGenClean.Ethnicity.check_ethnicity:safe_main",
                "pyGenClean_find_outliers=pyGenClean.Ethnicity.find_outliers:safe_main",
                "pyGenClean_flag_maf_zero=pyGenClean.FlagMAF.flag_maf_zero:safe_main",
                "pyGenClean_flag_hw=pyGenClean.FlagHW.flag_hw:safe_main",
                "pyGenClean_compare_gold_standard=pyGenClean.Misc.compare_gold_standard:safe_main",
                "pyGenClean_compare_bim=pyGenClean.PlinkUtils.compare_bim:safe_main",
                "pyGenClean_plot_MDS=pyGenClean.PlinkUtils.plot_MDS_standalone:safe_main",
                "pyGenClean_subset_data=pyGenClean.PlinkUtils.subset_data:safe_main",
                "pyGenClean_plot_eigenvalues=pyGenClean.Ethnicity.plot_eigenvalues:safe_main",
                "pyGenClean_merge_reports=pyGenClean.LaTeX.merge_reports:safe_main",
                "pyGenClean_check_contamination=pyGenClean.Contamination.contamination:safe_main",
            ],
        },
        install_requires=["matplotlib >= 1.2.0", "numpy >= 1.6.2",
                          "scipy >= 0.11.0", "scikit-learn >= 0.12.1",
                          "Jinja2 >= 2.7.3"],
        packages=["pyGenClean", "pyGenClean.Ethnicity", "pyGenClean.PlateBias",
                  "pyGenClean.DupSamples", "pyGenClean.SexCheck",
                  "pyGenClean.MarkerMissingness", "pyGenClean.FlagMAF",
                  "pyGenClean.FlagHW", "pyGenClean.RelatedSamples",
                  "pyGenClean.DupSNPs", "pyGenClean.Misc", "pyGenClean.LaTeX",
                  "pyGenClean.HeteroHap", "pyGenClean.SampleMissingness",
                  "pyGenClean.NoCallHetero", "pyGenClean.PlinkUtils",
                  "pyGenClean.Contamination"],
        package_data={"pyGenClean.LaTeX": ["templates/*.tex"]},
        classifiers=[
            "Development Status :: 5 - Production/Stable",
            "Intended Audience :: Science/Research",
            "License :: OSI Approved :: GNU General Public License (GPL)",
            "Operating System :: Unix",
            "Operating System :: POSIX :: Linux",
            "Operating System :: MacOS :: MacOS X",
            "Operating System :: Microsoft :: Windows",
            "Programming Language :: Python",
            "Programming Language :: Python :: 2.7",
            "Topic :: Scientific/Engineering :: Bio-Informatics"
        ],
        keywords="bioinformatics quality control genetic",
    )
    return


if __name__ == "__main__":
    setup_package()
