#!/usr/bin/env python
"""Setup script for pyGenClean."""

# How to build source distribution
#   - python setup.py sdist --format bztar
#   - python setup.py sdist --format gztar
#   - python setup.py sdist --format zip
#   - python setup.py bdist --format msi

# How to build for conda
#   - python setup.py bdist_conda
#   - conda convert -p all /PATH/TO/FILE -o conda_dist
#   - cd conda_dist && conda index *


import sys
from os import path

from setuptools import setup, find_packages


MAJOR = 2
MINOR = 0
MICRO = "0b1"
VERSION = "{0}.{1}.{2}".format(MAJOR, MINOR, MICRO)


def check_python_version():
    """Checks the python version, error if < 3.5."""
    python_major, python_minor = sys.version_info[:2]

    if python_major != 3 or python_minor < 5:
        sys.stderr.write("pyGenClean requires python 3.5")
        sys.exit(1)


def write_version_file():
    """Writes the version to file."""
    fn = path.join(
        path.dirname(path.abspath(__file__)), "pygenclean", "version.py",
    )

    content = ("\n# THIS FILE WAS GENERATED AUTOMATICALLY BY PYGENCLEAN\n"
               'pygenclean_version = "{version}"\n')

    with open(fn, "w") as f:
        f.write(content.format(version=VERSION))


def setup_package():
    """Setup the package."""
    # Checking the python version
    check_python_version()

    # Saving the version into a file
    write_version_file()

    setup(
        name="pyGenClean",
        version=VERSION,
        description="Automated data clean up pipeline for genetic data",
        long_description="This package provides tools to automatically "
                         "perform genetic data clean up (QC steps) prior to "
                         "a genome-wide association study.",
        author="Louis-Philippe Lemieux Perreault",
        author_email="louis-philippe.lemieux.perreault@umontreal.ca",
        url="https://github.com/lemieuxl/pyGenClean",
        license="GPL",
        entry_points={
            "console_scripts": [
                "pyGenClean=pygenclean.cli:main",
            ],
        },
        install_requires=[],
        packages=find_packages(),
        include_package_data=True,
        classifiers=[
            "Development Status :: 5 - Production/Stable",
            "Intended Audience :: Science/Research",
            "License :: OSI Approved :: GNU General Public License (GPL)",
            "Operating System :: Unix",
            "Operating System :: POSIX :: Linux",
            "Operating System :: MacOS :: MacOS X",
            "Operating System :: Microsoft :: Windows",
            "Programming Language :: Python",
            "Programming Language :: Python :: 3",
            "Topic :: Scientific/Engineering :: Bio-Informatics"
        ],
        keywords="bioinformatics quality control genetic",
    )


if __name__ == "__main__":
    setup_package()
