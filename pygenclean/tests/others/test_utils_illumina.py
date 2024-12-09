"""Tests the illumina sub module from the utils module."""


import gzip
import random
from pathlib import Path

from ...utils import illumina


def test_get_data_position_with_tag_uncompressed(tmp_path: Path):
    """Test the function where we find the number of lines to skip."""
    # Writing the data
    dummy_file = tmp_path / "dummy_file.txt"
    nb_skip = random.randint(1, 10)
    with open(dummy_file, "w") as f:
        # Printing between 1 and 10 lines
        for i in range(nb_skip):
            print(f"line {i + 1}", file=f)

        # Priting the data tag
        print("[Data]", file=f)

        # Printing between 10 and 100 extra lines
        for i in range(random.randint(10, 100)):
            print(f"extra line {i + 1}", file=f)

    # Checking the results
    assert illumina.get_data_position(dummy_file) == nb_skip + 1


def test_get_data_position_with_tag_compressed(tmp_path: Path):
    """Test the function where we find the number of lines to skip."""
    # Writing the data
    dummy_file = tmp_path / "dummy_file.txt.gz"
    nb_skip = random.randint(1, 10)
    with gzip.open(dummy_file, "wt") as f:
        # Printing between 1 and 10 lines
        for i in range(nb_skip):
            print(f"line {i + 1}", file=f)

        # Priting the data tag
        print("[Data]", file=f)

        # Printing between 10 and 100 extra lines
        for i in range(random.randint(10, 100)):
            print(f"extra line {i + 1}", file=f)

    # Checking the results
    assert illumina.get_data_position(dummy_file) == nb_skip + 1


def test_get_data_position_without_tag_uncompressed(tmp_path: Path):
    """Test the function where we find the number of lines to skip."""
    # Writing the data
    dummy_file = tmp_path / "dummy_file.txt"
    with open(dummy_file, "w") as f:
        # Printing between 1 and 10 lines
        for i in range(random.randint(110, 120)):
            print(f"line {i + 1}", file=f)

    # Checking the results
    assert illumina.get_data_position(dummy_file) == 0


def test_get_data_position_without_tag_compressed(tmp_path: Path):
    """Test the function where we find the number of lines to skip."""
    # Writing the data
    dummy_file = tmp_path / "dummy_file.txt.gz"
    with gzip.open(dummy_file, "wt") as f:
        # Printing between 1 and 10 lines
        for i in range(random.randint(110, 120)):
            print(f"line {i + 1}", file=f)

    # Checking the results
    assert illumina.get_data_position(dummy_file) == 0
