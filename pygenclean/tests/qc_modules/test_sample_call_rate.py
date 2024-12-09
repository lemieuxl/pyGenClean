"""Test the sample call rate QC module."""


import random
from argparse import Namespace
from pathlib import Path

import numpy as np
import pytest
from pyplink import PyPlink
from pytest_mock import MockerFixture

from ...error import ProgramError
from ...qc_modules.sample_call_rate import sample_call_rate
from .. import utils


def test_main_no_samples_removed(tmp_path: Path):
    """Test the main function (no samples were removed).

    No samples will be removed, since there is no missing value. This is to
    check what happens when no samples have a high enough missing rate.

    """
    # Generating the genotypes
    genotypes = utils.generate_genotypes(nb_variants=100, nb_samples=10)

    # Generating the input Plink files from the genotypes
    prefix = str(tmp_path / "original")
    utils.generate_plink_files(prefix, genotypes)

    # The output prefix
    out = str(tmp_path / "out")

    # Executing the script
    argv = [
        "--bfile", prefix,
        "--mind", "0.02",
        "--out", str(out),
    ]
    sample_call_rate.main(argv=argv)

    # There should not be any '.irem' file
    assert not (tmp_path / (out + ".irem")).is_file()

    # The BED, BIM and FAM files should be identical
    for suffix in (".bed", ".bim", ".fam"):
        with open(prefix + suffix, "rb") as f1, open(out + suffix, "rb") as f2:
            assert f1.read() == f2.read()


def test_main_samples_removed(tmp_path: Path):
    """Test the main function (samples were removed).

    There will be multiple samples with missing rate high enough, and one
    sample with missing rate equals to the 'mind' value. This is to check if
    the threshold is inclusiive or not.

    """
    # Generating the genotypes
    genotypes = utils.generate_genotypes(nb_variants=100, nb_samples=100)

    # Randomly selecting 10 variants
    indices = sorted(random.sample(range(genotypes.shape[1]), 10))

    # The first 9 sammples will have a missing rate > 0.02
    for j in indices[:-1]:
        nb_missing = random.randint(3, 100)
        for i in random.sample(range(genotypes.shape[0]), nb_missing):
            genotypes[i, j] = -1

    # The last sample will have exactly 2% missing
    for i in random.sample(range(genotypes.shape[1]), 2):
        genotypes[j, indices[-1]] = -1

    # Generating the input Plink files from the genotypes
    prefix = str(tmp_path / "original")
    utils.generate_plink_files(prefix, genotypes)

    # The output prefix
    out = str(tmp_path / "out")

    # Executing the script
    argv = [
        "--bfile", prefix,
        "--mind", "0.02",
        "--out", str(out),
    ]
    sample_call_rate.main(argv=argv)

    # The '.irem' file should have 9 samples
    with open(out + ".irem") as f:
        removed_samples = set(
            tuple(line.split()) for line in f.read().splitlines()
        )
        assert removed_samples == {(f"s{i + 1}", ) * 2 for i in indices[:-1]}

    # Checking the imiss file
    with open(out + ".imiss") as f:
        for i, line in enumerate(f):
            # header
            if i == 0:
                continue

            fid, iid, miss_pheno, n_miss, n_geno, f_miss = line.split()

            # The index of the current samples
            sample_index = indices[i - 1]

            # Finding the number of missing genotypes
            nb_missing = np.count_nonzero(genotypes[:, sample_index] == -1)
            missing_rate = nb_missing / genotypes.shape[0]

            # Checking the values
            assert fid == iid
            assert iid == f"s{sample_index + 1}"
            assert miss_pheno == "Y"
            assert int(n_miss) == nb_missing
            assert int(n_geno) == genotypes.shape[0]
            assert float(f_miss) == pytest.approx(missing_rate)

    # Comparing the FAM file
    with open(prefix + ".fam") as f1, open(out + ".fam") as f2:
        fam1 = [
            line for i, line in enumerate(f1.read().splitlines())
            if i not in indices[:-1]
        ]
        fam2 = f2.read().splitlines()
        assert fam1 == fam2

    # The BIM files should be the same, except for allele frequencies (because
    # samples were removed). We're excluding the samples
    genotypes = np.delete(genotypes, indices[:-1], 1)
    to_flip = []
    with open(prefix + ".bim") as f1, open(out + ".bim") as f2:
        content1 = [line.split() for line in f1.read().splitlines()]
        content2 = [line.split() for line in f2.read().splitlines()]
        assert len(content1) == len(content2)

        for i, (row1, row2) in enumerate(zip(content1, content2)):
            # Comparing the first part (chrom, pos, cm and pos)
            assert row1[:4] == row2[:4]

            # Generating the MAF for the new genotypes
            var_geno = genotypes[i, :]
            var_geno = var_geno[var_geno != -1]
            maf = np.mean(var_geno) / 2

            if maf > 0.5:
                assert row1[4] == row2[5]
                assert row1[5] == row2[4]
                to_flip.append(i)

            else:
                assert row1[4:] == row2[4:]

    # Comparing the BED files
    with PyPlink(prefix) as bed1, PyPlink(out) as bed2:
        bed1 = np.array([geno for _, geno in bed1])
        bed2 = np.array([geno for _, geno in bed2])

        # Removing the samples in the first file
        bed1 = np.delete(bed1, indices[:-1], 1)

        # Flipping markers to flip
        if len(to_flip) > 0:
            to_flip = np.array(to_flip)
            bed1[to_flip] = 2 - bed1[to_flip]

        assert np.all(bed1 == bed2)


def test_run_plink(mocker: MockerFixture):
    """Test the run Plink command."""
    # Mocking the 'execute_external_command'
    mocked_execute = mocker.patch(
        "pygenclean.qc_modules.sample_call_rate.sample_call_rate"
        ".execute_external_command",
    )

    # Executing the command
    args = Namespace(
        plink_107=False,
        bfile="dummy_bfile",
        mind=0.1,
        out="dummy_out",
    )
    sample_call_rate.run_plink(args)

    # Checking the command which was created
    mocked_execute.assert_called_once_with(
        command=[
            "plink1.9",
            "--noweb",
            "--bfile", "dummy_bfile",
            "--mind", "0.1",
            "--make-bed",
            "--out", "dummy_out",
        ],
    )


def test_run_plink_old(mocker: MockerFixture):
    """Test the run Plink command (with old Plink)."""
    # Mocking the 'execute_external_command'
    mocked_execute = mocker.patch(
        "pygenclean.qc_modules.sample_call_rate.sample_call_rate"
        ".execute_external_command",
    )

    # Executing the command
    args = Namespace(
        plink_107=True,
        bfile="dummy_bfile",
        mind=0.2,
        out="dummy_out",
    )
    sample_call_rate.run_plink(args)

    # Checking the command which was created
    mocked_execute.assert_called_once_with(
        command=[
            "plink",
            "--noweb",
            "--bfile", "dummy_bfile",
            "--mind", "0.2",
            "--make-bed",
            "--out", "dummy_out",
        ],
    )


def test_check_args_ok(mocker: MockerFixture):
    """Test the argument checker (all good)."""
    mocked_check_files = mocker.patch(
        "pygenclean.qc_modules.sample_call_rate.sample_call_rate.plink_utils"
        ".check_files",
        return_value=True,
    )

    args = Namespace(bfile="dummy_prefix", mind=0.1)
    sample_call_rate.check_args(args)
    mocked_check_files.assert_called_once_with(args.bfile)


def test_check_args_fail_bfile(mocker: MockerFixture):
    """Test the argument checker (bfile is not good)."""
    mocker.patch(
        "pygenclean.qc_modules.sample_call_rate.sample_call_rate.plink_utils"
        ".check_files",
        return_value=False,
    )

    args = Namespace(bfile="dummy_prefix", mind=0.1)
    with pytest.raises(ProgramError):
        sample_call_rate.check_args(args)


def test_check_args_fail_geno(mocker: MockerFixture):
    """Test the argument checker (mind is not good)."""
    mocker.patch(
        "pygenclean.qc_modules.sample_call_rate.sample_call_rate.plink_utils"
        ".check_files",
        return_value=True,
    )

    # Below 0
    args = Namespace(bfile="dummy_prefix", mind=-0.1)
    with pytest.raises(ProgramError):
        sample_call_rate.check_args(args)

    # Above 1
    args = Namespace(bfile="dummy_prefix", mind=1.1)
    with pytest.raises(ProgramError):
        sample_call_rate.check_args(args)
