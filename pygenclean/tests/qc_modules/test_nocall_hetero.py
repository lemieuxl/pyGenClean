"""Test the no call / heterozygous only QC module.

TODO: Variants on chromosome MT.

"""


import os
import random
from argparse import Namespace
from pathlib import Path

import numpy as np
import pytest
from pyplink import PyPlink
from pytest_mock import MockerFixture

from ...error import ProgramError
from ...qc_modules.nocall_hetero import nocall_hetero
from .. import utils


def test_main_no_variants_removed(tmp_path: Path):
    """Test the main function (no variants were removed).

    No variants will be removed, since there is no missing value. This is to
    check what happens when no variants have a high enough missing rate.

    """
    # Generating the genotypes
    genotypes = utils.generate_genotypes(nb_variants=10, nb_samples=100)

    # The last variant will be on the MT chromosome, and it'll be 100% hetero
    chromosomes = np.random.randint(1, 22, size=genotypes.shape[0])
    chromosomes[-1] = 26
    genotypes[-1, :] = 1

    # The last sample is always no call
    genotypes[:, -1] = -1

    # Generating the input Plink files from the genotypes
    prefix = str(tmp_path / "original")
    utils.generate_plink_files(prefix, genotypes, chromosomes=chromosomes)

    # The output prefix
    out = str(tmp_path / "out")

    # Executing the script
    argv = [
        "--bfile", prefix,
        "--out", str(out),
    ]
    nocall_hetero.main(argv=argv)

    # The '.all_failed' file should be empty
    with open(out + ".all_failed") as f:
        assert len(f.read()) == 0

    # The '.all_hetero' file should be empty
    with open(out + ".all_hetero") as f:
        assert len(f.read()) == 0

    # There should not be a '.exclude' file
    assert not os.path.isfile(out + ".exclude")

    # The BED, BIM and FAM files should be identical
    for suffix in (".bed", ".bim", ".fam"):
        with open(prefix + suffix, "rb") as f1, open(out + suffix, "rb") as f2:
            assert f1.read() == f2.read()


def test_main_variants_no_call_all_hetero_removed(tmp_path: Path):
    """Test the main function (variants were removed).

    There will be 10 variants is 100% missing genotypes, and 10 variants with
    100% heterozygous genotypes. The last variant will be 100% nocall / hetero,
    but will be located on chromosome MT, so it should remain in the dataset.

    """
    # Generating the genotypes
    genotypes = utils.generate_genotypes(nb_variants=100, nb_samples=100)

    # Randomly selecting 20 variants
    indices = random.sample(range(genotypes.shape[0] - 1), 20)

    # The first 10 will be 100% no call
    genotypes[indices[:10], :] = -1

    # The last 10 wiil be 100% heterozygous (with some no calls)
    for i in indices[10:]:
        genotypes[i, :] = random.choices([-1, 1], k=genotypes.shape[1])

    # The last variants will be on the MT chromosome
    chromosomes = np.random.randint(1, 22, size=genotypes.shape[0])
    chromosomes[-1] = 26
    genotypes[-1, :] = random.choices([-1, 1], k=genotypes.shape[1])

    # Generating the input Plink files from the genotypes
    prefix = str(tmp_path / "original")
    utils.generate_plink_files(prefix, genotypes, chromosomes=chromosomes)

    # The output prefix
    out = str(tmp_path / "out")

    # Executing the script
    argv = [
        "--bfile", prefix,
        "--out", str(out),
    ]
    nocall_hetero.main(argv=argv)

    # The '.all_failed' should contain the first 10 selected variants
    with open(out + ".all_failed") as f:
        expected = {f"var{i + 1}" for i in indices[:10]}
        assert set(f.read().splitlines()) == expected

    # The '.all_hetero' should contain the last 10 selected variants
    with open(out + ".all_hetero") as f:
        expected = {f"var{i + 1}" for i in indices[10:]}
        assert set(f.read().splitlines()) == expected

    # The '.exclude' file should contain the 20 selected variants
    with open(out + ".exclude") as f:
        assert set(f.read().splitlines()) == {f"var{i + 1}" for i in indices}

    # The FAM files should be identical
    with open(prefix + ".fam", "rb") as f1, open(out + ".fam", "rb") as f2:
        assert f1.read() == f2.read()

    # Comparing the BIM file
    with open(prefix + ".bim") as f1, open(out + ".bim") as f2:
        bim1 = [
            line for i, line in enumerate(f1.read().splitlines())
            if i not in indices
        ]
        bim2 = f2.read().splitlines()
        assert bim1 == bim2

    # Comparing the BED file
    with PyPlink(prefix) as bed1, PyPlink(out) as bed2:
        bed1 = np.array([
            geno for i, (_, geno) in enumerate(bed1)
            if i not in indices
        ])
        bed2 = np.array([geno for _, geno in bed2])
        assert np.all(bed1 == bed2)


def test_main_variants_no_call_removed(tmp_path: Path):
    """Test the main function (variants were removed).

    There will be 10 variants is 100% missing genotypes, and 10 variants with
    100% heterozygous genotypes. Thos variants will stay in the dataset because
    of the '--keep-het' flag. The last variant will be 100% nocall / hetero,
    but will be located on chromosome MT, so it should remain in the dataset.

    """
    # Generating the genotypes
    genotypes = utils.generate_genotypes(nb_variants=100, nb_samples=100)

    # Randomly selecting 20 variants
    indices = random.sample(range(genotypes.shape[0] - 1), 20)

    # The first 10 will be 100% no call
    genotypes[indices[:10], :] = -1

    # The last 10 wiil be 100% heterozygous (with some no calls)
    for i in indices[10:]:
        genotypes[i, :] = random.choices([-1, 1], k=genotypes.shape[1])

    # The last variants will be on the MT chromosome
    chromosomes = np.random.randint(1, 22, size=genotypes.shape[0])
    chromosomes[-1] = 26
    genotypes[-1, :] = random.choices([-1, 1], k=genotypes.shape[1])

    # Generating the input Plink files from the genotypes
    prefix = str(tmp_path / "original")
    utils.generate_plink_files(prefix, genotypes, chromosomes=chromosomes)

    # The output prefix
    out = str(tmp_path / "out")

    # Executing the script
    argv = [
        "--bfile", prefix,
        "--keep-het",
        "--out", str(out),
    ]
    nocall_hetero.main(argv=argv)

    # The '.all_failed' should contain the first 10 selected variants
    with open(out + ".all_failed") as f:
        expected = {f"var{i + 1}" for i in indices[:10]}
        assert set(f.read().splitlines()) == expected

    # The '.all_hetero' should contain the last 10 selected variants, even
    # though they were not removed
    with open(out + ".all_hetero") as f:
        expected = {f"var{i + 1}" for i in indices[10:]}
        assert set(f.read().splitlines()) == expected

    # The '.exclude' file should contain the first 10 selected variants
    with open(out + ".exclude") as f:
        expected = {f"var{i + 1}" for i in indices[:10]}
        assert set(f.read().splitlines()) == expected

    # The FAM files should be identical
    with open(prefix + ".fam", "rb") as f1, open(out + ".fam", "rb") as f2:
        assert f1.read() == f2.read()

    # Comparing the BIM file
    with open(prefix + ".bim") as f1, open(out + ".bim") as f2:
        bim1 = [
            line for i, line in enumerate(f1.read().splitlines())
            if i not in indices[:10]
        ]
        bim2 = f2.read().splitlines()
        assert bim1 == bim2

    # Comparing the BED file
    with PyPlink(prefix) as bed1, PyPlink(out) as bed2:
        bed1 = np.array([
            geno for i, (_, geno) in enumerate(bed1)
            if i not in indices[:10]
        ])
        bed2 = np.array([geno for _, geno in bed2])
        assert np.all(bed1 == bed2)


def test_check_args_ok(mocker: MockerFixture):
    """Test the argument checker (all good)."""
    mocked_check_files = mocker.patch(
        "pygenclean.qc_modules.nocall_hetero.nocall_hetero.plink_utils"
        ".check_files",
        return_value=True,
    )

    args = Namespace(bfile="dummy_prefix")
    nocall_hetero.check_args(args)
    mocked_check_files.assert_called_once_with(args.bfile)


def test_check_args_fail_bfile(mocker: MockerFixture):
    """Test the argument checker (bfile is not good)."""
    mocker.patch(
        "pygenclean.qc_modules.nocall_hetero.nocall_hetero.plink_utils"
        ".check_files",
        return_value=False,
    )

    args = Namespace(bfile="dummy_prefix")
    with pytest.raises(ProgramError) as program_error:
        nocall_hetero.check_args(args)

    assert str(program_error.value) == "dummy_prefix: no such binary files"
