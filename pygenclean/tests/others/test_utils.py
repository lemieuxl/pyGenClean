"""Test the utils module."""

import bz2
import gzip
import lzma
import random
import uuid
from pathlib import Path

import pytest

from ... import utils
from ...error import ProgramError


def test_decode_chrom():
    """Test the 'decode_chrom' function."""
    # Testing all chromosome by integer (we allow for 0, because of unknown
    # chromosome)
    for chrom in range(0, 27):
        assert chrom == utils.decode_chrom(str(chrom))

    # Testing the named chromosomes (X, Y, XY, M and MT)
    for chrom, str_chrom in zip(range(23, 27), ("X", "Y", "XY", "MT")):
        assert chrom == utils.decode_chrom(str_chrom)
        assert chrom == utils.decode_chrom(str_chrom.lower())

    # Testing invalid chromosomes
    with pytest.raises(ProgramError) as program_error:
        utils.decode_chrom("Z")
    assert program_error.value.message == "Z: invalid chromosome"

    # Testing invalid integers (below 0)
    with pytest.raises(ProgramError) as program_error:
        utils.decode_chrom("-1")
    assert program_error.value.message == "-1: invalid chromosome"

    # Testing invalid integers (above 26)
    with pytest.raises(ProgramError) as program_error:
        utils.decode_chrom("27")
    assert program_error.value.message == "27: invalid chromosome"


def test_decode_sex():
    """Test the 'decode_sex' function.

    There's only 3 possible values:
        1: Male
        2: Female
        Any other value: Unknown

    """
    assert "Male" == utils.decode_sex("1")
    assert "Female" == utils.decode_sex("2")
    assert "Unknown" == utils.decode_sex("0")
    assert "Unknown" == utils.decode_sex("No")


def test_is_compressed(tmp_path: Path):
    """Test the 'is_gzip', 'is_xz' and 'is_bzip2' functions.

    We will test each compression file with each function.

    """
    dummy_file = tmp_path / "test"

    # Testing GZIP file with all the functions
    with gzip.open(dummy_file, "wt") as f:
        for i in range(10):
            print(f"line_{i + 1}", file=f)

    assert utils.is_gzip(dummy_file)
    assert not utils.is_xz(dummy_file)
    assert not utils.is_bzip2(dummy_file)

    # Testing XZ file with all the functions
    with lzma.open(dummy_file, "wt") as f:
        for i in range(10):
            print(f"line_{i + 1}", file=f)

    assert not utils.is_gzip(dummy_file)
    assert utils.is_xz(dummy_file)
    assert not utils.is_bzip2(dummy_file)

    # Testing BZIP2 files with all the functions
    with bz2.open(dummy_file, "wt") as f:
        for i in range(10):
            print(f"line_{i + 1}", file=f)

    assert not utils.is_gzip(dummy_file)
    assert not utils.is_xz(dummy_file)
    assert utils.is_bzip2(dummy_file)


def test_get_open_func(tmp_path: Path):
    """Test the 'get_open_func' function."""
    dummy_file = tmp_path / "test"

    # First, trying a normal file
    content = str(uuid.uuid4())
    with open(dummy_file, "w") as f:
        print(content, file=f)

    with utils.get_open_func(dummy_file)() as f:
        assert f.read() == content + "\n"

    # Then, trying a GZIP file
    content = str(uuid.uuid4())
    with gzip.open(dummy_file, "wt") as f:
        print(content, file=f)

    with utils.get_open_func(dummy_file)() as f:
        assert f.read() == content + "\n"

    # Then, trying a XZ file
    content = str(uuid.uuid4())
    with lzma.open(dummy_file, "wt") as f:
        print(content, file=f)

    with utils.get_open_func(dummy_file)() as f:
        assert f.read() == content + "\n"

    # Then, trying a BZIP2 file
    content = str(uuid.uuid4())
    with bz2.open(dummy_file, "wt") as f:
        print(content, file=f)

    with utils.get_open_func(dummy_file)() as f:
        assert f.read() == content + "\n"


def test_split_extra_args():
    """Test the 'split_extra_args' function."""
    assert not utils.split_extra_args(None)
    assert utils.split_extra_args("--bfile 'a prefix'; other bad *") == [
        "--bfile", "'a prefix;'", "other", "bad", "'*'",
    ]


def test_flip_alleles():
    """Test the 'flip_alleles' function."""
    assert utils.flip_alleles("AC0") == set("TG0")
    assert utils.flip_alleles("TG0") == set("AC0")

    with pytest.raises(ProgramError) as program_error:
        utils.flip_alleles("ACGTZ")
    assert program_error.value.message == "ACGTZ: unknown alleles"


def test_compatible_alleles():
    """Test the 'compatible_alleles' function."""
    # Same alleles
    assert utils.compatible_alleles({"A", "C"}, {"A", "C"})
    assert utils.compatible_alleles({"A", "C"}, {"0", "C"})
    assert utils.compatible_alleles({"A", "C"}, {"A", "0"})
    assert utils.compatible_alleles({"A", "C"}, {"0"})

    # Complement alleles
    assert utils.compatible_alleles({"A", "C"}, {"T", "G"})
    assert utils.compatible_alleles({"A", "C"}, {"0", "G"})
    assert utils.compatible_alleles({"A", "C"}, {"T", "0"})
    assert utils.compatible_alleles({"A", "C"}, {"0"})

    # Incompatible alleles
    assert not utils.compatible_alleles({"A", "C"}, {"T", "C"})
    assert not utils.compatible_alleles({"C", "G"}, {"T", "A"})


def test_check_allele_status():
    """Test the 'check_allele_status' function.

    Note that it is impossible to get 'problem' when one set of alleles has 2
    items, and the other has 1 item.

    """
    # Simple cases
    assert utils.check_allele_status({"A", "C"}, {"A", "C"}) is None
    assert utils.check_allele_status({"A"}, {"T"}) == "homo_flip"
    assert utils.check_allele_status({"A", "C"}, {"T", "G"}) == "flip"
    assert utils.check_allele_status({"A", "C"}, {"A", "G"}) == "problem"
    assert utils.check_allele_status({"A"}, {"C", "A"}) == "homo_hetero"
    assert utils.check_allele_status({"A"}, {"T", "C"}) == "homo_hetero_flip"
    assert utils.check_allele_status({"A"}, {"A", "T"}) == "homo_hetero"

    # With missing values
    assert utils.check_allele_status({"A", "0"}, {"A", "0"}) is None
    assert utils.check_allele_status({"A"}, {"T", "0"}) == "homo_flip"
    assert utils.check_allele_status({"A", "0"}, {"T"}) == "homo_flip"
    assert utils.check_allele_status({"A", "0"}, {"G", "0"}) == "problem"


def test_count_lines(tmp_path: Path):
    """Test the 'count_lines' function."""
    dummy_file = tmp_path / "dummy.txt"

    # Each line has a new line character
    nb_lines = random.randint(100, 200)
    with open(dummy_file, "w") as f:
        for i in range(nb_lines):
            print(f"line {i + 1}", file=f)
    assert utils.count_lines(dummy_file) == nb_lines

    # Last line doesn't have a new line character
    nb_lines = random.randint(100, 200)
    with open(dummy_file, "w") as f:
        for i in range(nb_lines):
            f.write(f"line {i + 1}\n")
        f.write("final line without new line")
    assert utils.count_lines(dummy_file) == nb_lines + 1
