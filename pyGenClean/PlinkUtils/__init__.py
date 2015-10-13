
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
import tempfile
from subprocess import Popen, PIPE

from ..pipeline_error import ProgramError


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


def get_plink_version():
    """Gets the Plink version from the binary.

    :returns: the version of the Plink software
    :rtype: str

    This function uses :py:class:`subprocess.Popen` to gather the version of
    the Plink binary. Since executing the software to gather the version
    creates an output file, it is deleted.

    .. warning::
        This function only works as long as the version is returned as
        ``| PLINK! | NNN |`` (where, ``NNN`` is the version), since we use
        regular expresion to extract the version number from the standard
        output of the software.

    """
    # Running the command
    tmp_fn = None
    with tempfile.NamedTemporaryFile(delete=False) as tmpfile:
        tmp_fn = tmpfile.name + "_pyGenClean"

    # The command to run
    command = ["plink", "--noweb", "--out", tmp_fn]

    output = None
    try:
        proc = Popen(command, stdout=PIPE, stderr=PIPE)
        output = proc.communicate()[0].decode()
    except OSError:
        raise ProgramError("plink: command not found")

    # Deleting the output file automatically created by Plink
    if os.path.isfile(tmp_fn + ".log"):
        os.remove(tmp_fn + ".log")

    # Finding the version
    version = re.search(r"\|\s+PLINK!\s+\|\s+(\S+)\s+\|", output)
    if version is None:
        version = "unknown"
    else:
        version = version.group(1)

    return version
