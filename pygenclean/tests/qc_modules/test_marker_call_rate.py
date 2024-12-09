"""Test the marker call rate QC module.

TODO: Variants on chromosome X and Y...

"""


import random
from argparse import Namespace
from pathlib import Path

import numpy as np
import pytest
from pyplink import PyPlink
from pytest_mock import MockerFixture

from ...error import ProgramError
from ...qc_modules.marker_call_rate import marker_call_rate
from .. import utils


def test_main_no_variants_removed(tmp_path: Path):
    """Test the main function (no variants were removed).

    No variants will be removed, since there is no missing value. This is to
    check what happens when no variants have a high enough missing rate.

    """
    # Generating the genotypes
    genotypes = utils.generate_genotypes(nb_variants=10, nb_samples=100)

    # Generating the input Plink files from the genotypes
    prefix = str(tmp_path / "original")
    utils.generate_plink_files(prefix, genotypes)

    # The output prefix
    out = str(tmp_path / "out")

    # Executing the script
    argv = [
        "--bfile", prefix,
        "--geno", "0.02",
        "--out", str(out),
    ]
    marker_call_rate.main(argv=argv)

    # The '.removed_snps' file should be empty
    with open(out + ".removed_snps") as f:
        assert len(f.read()) == 0

    # The BED, BIM and FAM files should be identical
    for suffix in (".bed", ".bim", ".fam"):
        with open(prefix + suffix, "rb") as f1, open(out + suffix, "rb") as f2:
            assert f1.read() == f2.read()


def test_main_variants_removed(tmp_path: Path):
    """Test the main function (variants were removed).

    There will be multiple variants with missing rate high enough, and one
    variant with missing rate equals to the 'geno' value. This is to check if
    the threshold is inclusiive or not.

    """
    # Generating the genotypes
    genotypes = utils.generate_genotypes(nb_variants=100, nb_samples=100)

    # Randomly selecting 10 variants
    indices = random.sample(range(genotypes.shape[0]), 10)

    # The first 9 variants will have a missing rate > 0.02
    for i in indices[:-1]:
        nb_missing = random.randint(3, 100)
        for j in random.sample(range(genotypes.shape[1]), nb_missing):
            genotypes[i, j] = -1

    # The last variant will have exactly 2% missing
    for j in random.sample(range(genotypes.shape[1]), 2):
        genotypes[indices[-1], j] = -1

    # Generating the input Plink files from the genotypes
    prefix = str(tmp_path / "original")
    utils.generate_plink_files(prefix, genotypes)

    # The output prefix
    out = str(tmp_path / "out")

    # Executing the script
    argv = [
        "--bfile", prefix,
        "--geno", "0.02",
        "--out", str(out),
    ]
    marker_call_rate.main(argv=argv)

    # The '.removed_snps' file should have 9 variants
    with open(out + ".removed_snps") as f:
        removed_variants = set(f.read().splitlines())
        assert removed_variants == {f"var{i + 1}" for i in indices[:-1]}

    # The FAM files should be identical
    with open(prefix + ".fam", "rb") as f1, open(out + ".fam", "rb") as f2:
        assert f1.read() == f2.read()

    # Comparing the BIM file
    with open(prefix + ".bim") as f1, open(out + ".bim") as f2:
        bim1 = [
            line for i, line in enumerate(f1.read().splitlines())
            if i not in indices[:-1]
        ]
        bim2 = f2.read().splitlines()
        assert bim1 == bim2

    # Comparing the BED file
    with PyPlink(prefix) as bed1, PyPlink(out) as bed2:
        bed1 = np.array([
            geno for i, (_, geno) in enumerate(bed1)
            if i not in indices[:-1]
        ])
        bed2 = np.array([geno for _, geno in bed2])
        assert np.all(bed1 == bed2)


def test_compare_bim(mocker: MockerFixture, tmp_path: Path):
    """Just test something."""
    removed_snps = {"snp1", "snp2", "snp3"}

    # Mocking the compare_bim from the plink module
    mocked_compare_bim = mocker.patch(
        "pygenclean.qc_modules.marker_call_rate.marker_call_rate.plink_utils"
        ".compare_bim",
        return_value=(removed_snps, None, set()),
    )

    # The final output file
    output = tmp_path / "fake_2.removed_snps"

    # The name of the dummy file
    dummy_bim = tmp_path / "fake_2.bim"
    args = Namespace(bfile="fake_1", out=str(dummy_bim.with_suffix("")))
    marker_call_rate.compare_bim(args)

    # Checking the arguments used (names of the BIM files) are valid
    mocked_compare_bim.assert_called_once_with(
        bim_a="fake_1.bim", bim_b=str(dummy_bim),
    )

    # Making sure the file exists
    assert output.is_file()

    # Checing its content
    with open(output) as f:
        assert set(f.read().splitlines()) == removed_snps


def test_compare_bim_fails(mocker: MockerFixture):
    """Just test something."""
    just_in_1 = {"snp1", "snp2", "snp3"}
    just_in_2 = {"snp4"}

    # Mocking the compare_bim from the plink module
    mocker.patch(
        "pygenclean.utils.plink.compare_bim",
        return_value=(just_in_1, None, just_in_2),
    )

    args = Namespace(bfile="fake_1", out="fake_2")
    with pytest.raises(ProgramError):
        marker_call_rate.compare_bim(args)


def test_run_plink(mocker: MockerFixture):
    """Test the run Plink command."""
    # Mocking the 'execute_external_command'
    mocked_execute = mocker.patch(
        "pygenclean.qc_modules.marker_call_rate.marker_call_rate"
        ".execute_external_command",
    )

    # Executing the command
    args = Namespace(
        plink_107=False,
        bfile="dummy_bfile",
        geno=0.2,
        out="dummy_out",
    )
    marker_call_rate.run_plink(args)

    # Checking the command which was created
    mocked_execute.assert_called_once_with(
        command=[
            "plink1.9",
            "--noweb",
            "--bfile", "dummy_bfile",
            "--geno", "0.2",
            "--make-bed",
            "--out", "dummy_out",
        ],
    )


def test_run_plink_old(mocker: MockerFixture):
    """Test the run Plink command (with old Plink)."""
    # Mocking the 'execute_external_command'
    mocked_execute = mocker.patch(
        "pygenclean.qc_modules.marker_call_rate.marker_call_rate"
        ".execute_external_command",
    )

    # Executing the command
    args = Namespace(
        plink_107=True,
        bfile="dummy_bfile",
        geno=0.3,
        out="dummy_out",
    )
    marker_call_rate.run_plink(args)

    # Checking the command which was created
    mocked_execute.assert_called_once_with(
        command=[
            "plink",
            "--noweb",
            "--bfile", "dummy_bfile",
            "--geno", "0.3",
            "--make-bed",
            "--out", "dummy_out",
        ],
    )


def test_check_args_ok(mocker: MockerFixture):
    """Test the argument checker (all good)."""
    mocked_check_files = mocker.patch(
        "pygenclean.qc_modules.marker_call_rate.marker_call_rate.plink_utils"
        ".check_files",
        return_value=True,
    )

    args = Namespace(bfile="dummy_prefix", geno=0.1)
    marker_call_rate.check_args(args)
    mocked_check_files.assert_called_once_with(args.bfile)


def test_check_args_fail_bfile(mocker: MockerFixture):
    """Test the argument checker (bfile is not good)."""
    mocker.patch(
        "pygenclean.qc_modules.marker_call_rate.marker_call_rate.plink_utils"
        ".check_files",
        return_value=False,
    )

    args = Namespace(bfile="dummy_prefix", geno=0.1)
    with pytest.raises(ProgramError):
        marker_call_rate.check_args(args)


def test_check_args_fail_geno(mocker: MockerFixture):
    """Test the argument checker (geno is not good)."""
    mocker.patch(
        "pygenclean.qc_modules.marker_call_rate.marker_call_rate.plink_utils"
        ".check_files",
        return_value=True,
    )

    # Below 0
    args = Namespace(bfile="dummy_prefix", geno=-0.1)
    with pytest.raises(ProgramError):
        marker_call_rate.check_args(args)

    # Above 1
    args = Namespace(bfile="dummy_prefix", geno=1.1)
    with pytest.raises(ProgramError):
        marker_call_rate.check_args(args)
