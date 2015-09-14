
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


import re


# The module version
__version__ = "0.9.2"


def get_version():
    """Returns the version of the module.

    :returns: (major, minor, micro)

    """

    return tuple(__version__.split("."))


def createRowFromPlinkSpacedOutput(line):
    """Remove leading spaces and change spaces to tabs.

    :param: line: a line from a Plink's report file.

    :type: line: str

    :returns: an array containing each field from the input line.

    Plink's output files are usually created so that they are human readable.
    Hence, instead of separating fields using tabulation, it uses a certain
    amount of spaces to create columns. Using the :py:mod:`re` module, the
    fields are split.

    .. testsetup::

        from pyGenClean.PlinkUtils import createRowFromPlinkSpacedOutput

    .. doctest::

        >>> line = " CHR               SNP         BP   A1      A2"
        >>> createRowFromPlinkSpacedOutput(line)
        ['CHR', 'SNP', 'BP', 'A1', 'A2']

    """

    return re.sub(
        "  *",
        "\t",
        re.sub("^  *", "", line.rstrip("\r\n")),
    ).split("\t")
