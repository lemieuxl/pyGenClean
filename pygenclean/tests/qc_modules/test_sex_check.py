"""Test the sex check QC module."""


import random
from argparse import Namespace
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from pyplink import PyPlink
from pytest_mock import MockerFixture

from ...error import ProgramError
from ...qc_modules.sex_check import sex_check
from ..utils import generate_genotypes, generate_plink_files


def test_check_args_ok(mocker: MockerFixture):
    """Test the argument checker (all good)."""
    mocked_check_files = mocker.patch(
        "pygenclean.qc_modules.sex_check.sex_check.plink_utils.check_files",
        return_value=True,
    )

    args = Namespace(
        bfile="dummy_prefix",
        nb_chr_23=1000,
        intensity_plot=False,
        baf_lrr=False,
    )
    sex_check.check_args(args)
    mocked_check_files.assert_called_once_with(args.bfile)


def test_check_args_fail_bfile(mocker: MockerFixture):
    """Test the argument checker (bfile is not good)."""
    mocker.patch(
        "pygenclean.qc_modules.sex_check.sex_check.plink_utils.check_files",
        return_value=False,
    )

    args = Namespace(
        bfile="dummy_prefix",
        nb_chr_23=1000,
        intensity_plot=False,
        baf_lrr=False,
    )
    with pytest.raises(ProgramError) as program_error:
        sex_check.check_args(args)

    assert str(program_error.value) == "dummy_prefix: no such binary files"


def test_check_args_nb_chr_23(mocker: MockerFixture):
    """Test the argument checker (negative number of markers on chr X)."""
    mocker.patch(
        "pygenclean.qc_modules.sex_check.sex_check.plink_utils.check_files",
        return_value=True,
    )

    args = Namespace(
        bfile="dummy_prefix",
        nb_chr_23=-1,
        intensity_plot=False,
        baf_lrr=False,
    )
    with pytest.raises(ProgramError) as program_error:
        sex_check.check_args(args)

    assert str(program_error.value) == (
        "-1: number of markers on chr23 must be positive"
    )


def test_check_args_with_intensities(mocker: MockerFixture):
    """Test the argument checker (itensities provided)."""
    mocker.patch(
        "pygenclean.qc_modules.sex_check.sex_check.plink_utils.check_files",
        return_value=True,
    )

    args = Namespace(
        bfile="dummy_prefix",
        nb_chr_23=1000,
        intensity_plot=True,
        sex_intensities="dummy_file",
        intensity_plot_extra_args="--dummy 1; --multiple '3 4 5'",
        baf_lrr=False,
    )

    sex_check.check_args(args)
    assert args.intensity_plot_extra_args == [
        "--dummy", "'1;'", "--multiple", "'3 4 5'",
    ]


def test_check_args_no_intensities(mocker: MockerFixture):
    """Test the argument checker (no itensities provided)."""
    mocker.patch(
        "pygenclean.qc_modules.sex_check.sex_check.plink_utils.check_files",
        return_value=True,
    )

    args = Namespace(
        bfile="dummy_prefix",
        nb_chr_23=1000,
        intensity_plot=True,
        sex_intensities=None,
        baf_lrr=False,
    )
    with pytest.raises(ProgramError) as program_error:
        sex_check.check_args(args)

    assert str(program_error.value) == (
        "Asking for intensity plot, but no '--intensities' provided"
    )


def test_check_args_with_baf_lrr(mocker: MockerFixture):
    """Test the argument checker (with intensity directory)."""
    mocker.patch(
        "pygenclean.qc_modules.sex_check.sex_check.plink_utils.check_files",
        return_value=True,
    )

    args = Namespace(
        bfile="dummy_prefix",
        nb_chr_23=1000,
        intensity_plot=False,
        baf_lrr=True,
        per_sample_baf_lrr_dir="dummy_directory",
        baf_lrr_extra_args="--dummy 1; --multiple '3 4 5'",
    )

    sex_check.check_args(args)
    assert args.baf_lrr_extra_args == [
        "--dummy", "'1;'", "--multiple", "'3 4 5'",
    ]


def test_check_args_no_baf_lrr(mocker: MockerFixture):
    """Test the argument checker (no itensities provided)."""
    mocker.patch(
        "pygenclean.qc_modules.sex_check.sex_check.plink_utils.check_files",
        return_value=True,
    )

    args = Namespace(
        bfile="dummy_prefix",
        nb_chr_23=1000,
        intensity_plot=False,
        baf_lrr=True,
        per_sample_baf_lrr_dir=None,
    )
    with pytest.raises(ProgramError) as program_error:
        sex_check.check_args(args)

    assert str(program_error.value) == (
        "Asking for BAF & LRR plot, but no '--intensity-dir' provided"
    )


def test_get_sex_mismatch(tmp_path: Path):
    """Test the 'get_sex_mismatch' function.

    The default values for `female_f` and `male_f` are 0.3 and 0.7,
    respectively.

    There are "problems" (STATUS != OK) and mismatches (PEDSEX != SNPSEX).

    We consider females (PEDSEX = 2) to be "OK" if SNPSEX is unknown but F is
    below the threshold.

    We consider males (PEDSEX = 1) to be "OK" if SNPSEX is unknown but F is
    above the threshold.

    """
    dummy_sexcheck = tmp_path / "dummy.sexcheck"
    list_problem_sex = tmp_path / "dummy.list_problem_sex"
    list_problem_sex_ids = tmp_path / "dummy.list_problem_sex_ids"

    # Generating content (3 samples for each category)
    nb_samples = 3
    content = []
    expected_problems = []

    # OK males
    for i in range(nb_samples):
        f_value = random.uniform(0.9, 1.0)
        content.append((f"F{i}", f"I{i}", "1", "1", "OK", f"{f_value:.4f}"))

    # OK females
    for i in range(i, i + nb_samples):
        f_value = random.uniform(0.0, 0.2)
        content.append((f"F{i}", f"I{i}", "2", "2", "OK", f"{f_value:.4f}"))

    # Problematic males
    for i in range(i, i + nb_samples):
        f_value = random.uniform(0.0, 0.2)
        content.append(
            (f"F{i}", f"I{i}", "1", "2", "PROBLEM", f"{f_value:.4f}"),
        )
        expected_problems.append(content[-1])

    # Problematic females
    for i in range(i, i + nb_samples):
        f_value = random.uniform(0.9, 1.0)
        content.append(
            (f"F{i}", f"I{i}", "2", "1", "PROBLEM", f"{f_value:.4f}"),
        )
        expected_problems.append(content[-1])

    # Problematic females, but F < 0.3 (rescued)
    for i in range(i, i + nb_samples):
        f_value = random.uniform(0.0, 0.299)
        content.append(
            (f"F{i}", f"I{i}", "2", "0", "PROBLEM", f"{f_value:.4f}"),
        )

    # Problematic males, but F > 0.7 (rescued)
    for i in range(i, i + nb_samples):
        f_value = random.uniform(0.701, 1.0)
        content.append(
            (f"F{i}", f"I{i}", "1", "0", "PROBLEM", f"{f_value:.4f}"),
        )

    # Unknown with male SNPSEX and F > 0.7 (hence males)
    for i in range(i, i + nb_samples):
        f_value = random.uniform(0.75, 1.0)
        content.append(
            (f"F{i}", f"I{i}", "0", "1", "PROBLEM", f"{f_value:.4f}"),
        )
        expected_problems.append(content[-1])

    # Unknown with male SNPSEX and F < 0.3 (hence females)
    for i in range(i, i + nb_samples):
        f_value = random.uniform(0.0, 0.27)
        content.append(
            (f"F{i}", f"I{i}", "0", "2", "PROBLEM", f"{f_value:.4f}"),
        )
        expected_problems.append(content[-1])

    with open(dummy_sexcheck, "w") as f:
        fmt_string = "%4s %4s %7s %7s %8s %7s\n"

        header = ("FID", "IID", "PEDSEX", "SNPSEX", "STATUS", "F")
        f.write(fmt_string % header)

        for row in content:
            f.write(fmt_string % row)

    # Executing the function
    observed_mismatches = sex_check.get_sex_mismatch(
        filename=dummy_sexcheck,
        female_f=0.3,
        male_f=0.7,
        prefix=str(dummy_sexcheck.with_suffix("")),
    )

    # Checking the returned values
    assert observed_mismatches == {tuple(row[:2]) for row in expected_problems}

    # Checking the '.list_problem_sex' file
    assert list_problem_sex.is_file()
    with open(list_problem_sex) as f:
        # Reading the file
        observed = [
            tuple(line.split("\t")) for line in f.read().splitlines()
        ]

    # Checking the header
    assert observed[0] == ("FID", "IID", "PEDSEX", "SNPSEX", "STATUS", "F")

    # Checking the content
    for observed_row, expected_row in zip(observed[1:], expected_problems):
        # First 5 values are the same
        assert observed_row[:5] == expected_row[:5]

        # The last value (index 5) might have rounding issues (removing last 0)
        assert observed_row[5] == expected_row[5].rstrip("0")

    # Checking the '.list_problem_sex_ids' file
    assert list_problem_sex_ids.is_file()
    with open(list_problem_sex_ids) as f:
        observed = [tuple(line.split("\t")) for line in f.read().splitlines()]
    assert observed == [tuple(row[:2]) for row in expected_problems]


def test_get_sex_mismatch_no_input():
    """Test the 'get_sex_mismatch' function when no input file exists."""
    # Executing the function
    with pytest.raises(ProgramError) as program_error:
        sex_check.get_sex_mismatch(
            filename="dummy.sexcheck",
            female_f=0.3,
            male_f=0.7,
            prefix="dummy",
        )
    assert program_error.value.message == (
        "dummy.sexcheck: something went wrong with PLink"
    )


def test_get_genotypes(tmp_path: Path):
    """Test the 'get_genotypes' function."""
    # Generating genotypes
    genotypes = generate_genotypes(nb_samples=10, nb_variants=100)

    dummy_bfile = tmp_path / "dummy.bed"
    generate_plink_files(str(dummy_bfile.with_suffix("")), genotypes=genotypes)

    # We'll randomly select some variants and samples
    variants_indices = np.array(random.sample(
        range(genotypes.shape[0]), k=random.randint(40, 60),
    ))
    sample_indices = np.sort(np.array(random.sample(
        range(genotypes.shape[1]), k=random.randint(4, 7),
    )))

    # The selected samples (as booleans)
    selected_samples = np.zeros(genotypes.shape[1], dtype=bool)
    selected_samples[sample_indices] = True

    # The selected variants, as a list of string
    selected_variants = [f"var{i}" for i in variants_indices]

    # Executing the function
    with PyPlink(str(dummy_bfile.with_suffix(""))) as bed:
        observed = sex_check.get_genotypes(
            samples=selected_samples,
            markers=selected_variants,
            bed=bed,
        )

    # Checking the results
    assert np.all(observed == genotypes[variants_indices][:, sample_indices])


def test_write_no_call(tmp_path: Path):
    """Test the 'write_no_call' function."""
    dummy_file = tmp_path / "dummy.txt"

    # Generating genotypes for 30 samples
    nb_samples = 30
    genotypes = generate_genotypes(nb_samples=nb_samples, nb_variants=100)

    # Adding missing values
    all_nb_missing = {}
    for i in range(genotypes.shape[1]):
        nb_missing = random.randint(10, 30)
        genotypes[:, i][random.sample(range(100), k=nb_missing)] = -1
        all_nb_missing[(f"f{i}", f"s{i}")] = nb_missing

    # Generating a FAM file
    fam = pd.DataFrame({
        "fid": [f"f{i}" for i in range(nb_samples)],
        "iid": [f"s{i}" for i in range(nb_samples)],
        "father": [f"father_{i}" for i in range(nb_samples)],
        "mother": [f"mother_{i}" for i in range(nb_samples)],
        "gender": [random.randint(1, 2) for _ in range(nb_samples)],
        "status": [-9] * nb_samples,
    }).set_index(["fid", "iid"])

    # Executing the function
    sex_check.write_no_call(dummy_file, fam=fam, genotypes=genotypes)

    # Checking the file
    assert dummy_file.is_file()
    with open(dummy_file) as f:
        for i, line in enumerate(f):
            row = line.rstrip().split("\t")
            assert len(row) == 6

            # Checking the header
            if i == 0:
                assert row == [
                    "FID", "IID", "PEDSEX", "N_GENO", "N_MISS", "F_MISS",
                ]
                continue

            sample = tuple(row[:2])
            assert int(row[2]) == fam.loc[sample, "gender"]
            assert int(row[3]) == genotypes.shape[0]
            assert int(row[4]) == all_nb_missing[sample]
            assert float(row[5]) == all_nb_missing[sample] / genotypes.shape[0]


def test_write_no_call_only_one(tmp_path: Path):
    """Test the 'write_no_call' function (only one variant).

    This is for edge cases detection...

    """
    dummy_file = tmp_path / "dummy.txt"

    # Generating genotypes for 30 samples
    nb_samples = 30
    genotypes = generate_genotypes(nb_samples=nb_samples, nb_variants=1)

    # Adding missing values
    all_nb_missing = {}
    for i in range(genotypes.shape[1]):
        nb_missing = random.randint(0, 1)
        all_nb_missing[(f"f{i}", f"s{i}")] = nb_missing

        if nb_missing:
            genotypes[0, i] = -1

    # Generating a FAM file
    fam = pd.DataFrame({
        "fid": [f"f{i}" for i in range(nb_samples)],
        "iid": [f"s{i}" for i in range(nb_samples)],
        "father": [f"father_{i}" for i in range(nb_samples)],
        "mother": [f"mother_{i}" for i in range(nb_samples)],
        "gender": [random.randint(1, 2) for _ in range(nb_samples)],
        "status": [-9] * nb_samples,
    }).set_index(["fid", "iid"])

    # Executing the function
    sex_check.write_no_call(dummy_file, fam=fam, genotypes=genotypes)

    # Checking the file
    assert dummy_file.is_file()
    with open(dummy_file) as f:
        for i, line in enumerate(f):
            row = line.rstrip().split("\t")
            assert len(row) == 6

            # Checking the header
            if i == 0:
                assert row == [
                    "FID", "IID", "PEDSEX", "N_GENO", "N_MISS", "F_MISS",
                ]
                continue

            sample = tuple(row[:2])
            assert int(row[2]) == fam.loc[sample, "gender"]
            assert int(row[3]) == genotypes.shape[0]
            assert int(row[4]) == all_nb_missing[sample]
            assert float(row[5]) == all_nb_missing[sample] / genotypes.shape[0]


def test_write_no_call_no_variants(tmp_path: Path):
    """Test the 'write_no_call' function (no variants).

    This is for edge cases detection...

    """
    dummy_file = tmp_path / "dummy.txt"

    # Generating genotypes for 30 samples
    nb_samples = 30
    genotypes = generate_genotypes(nb_samples=nb_samples, nb_variants=10)

    # Removing variants from the genotypes
    genotypes = genotypes[np.zeros(genotypes.shape[0], dtype=bool), :]

    # Generating a FAM file
    fam = pd.DataFrame({
        "fid": [f"f{i}" for i in range(nb_samples)],
        "iid": [f"s{i}" for i in range(nb_samples)],
        "father": [f"father_{i}" for i in range(nb_samples)],
        "mother": [f"mother_{i}" for i in range(nb_samples)],
        "gender": [random.randint(1, 2) for _ in range(nb_samples)],
        "status": [-9] * nb_samples,
    }).set_index(["fid", "iid"])

    # Executing the function
    sex_check.write_no_call(dummy_file, fam=fam, genotypes=genotypes)

    # Checking the file
    assert dummy_file.is_file()
    with open(dummy_file) as f:
        for i, line in enumerate(f):
            row = line.rstrip().split("\t")
            assert len(row) == 6

            # Checking the header
            if i == 0:
                assert row == [
                    "FID", "IID", "PEDSEX", "N_GENO", "N_MISS", "F_MISS",
                ]
                continue

            sample = tuple(row[:2])
            assert int(row[2]) == fam.loc[sample, "gender"]
            assert row[3] == "0"
            assert row[4] == "0"
            assert row[5] == "-9"


def test_write_heterozygosity(tmp_path: Path):
    """Test the 'write_heterozygosity' function."""
    dummy_file = tmp_path / "dummy.txt"

    # Generating genotypes (only homozygous, i.e. 0 and 2)
    nb_samples = 30
    genotypes = 2 * np.random.randint(0, 2, size=(100, nb_samples),
                                      dtype=np.int8)

    # Adding heterozygous values
    all_nb_hetero = {}
    for i in range(genotypes.shape[1]):
        nb_hetero = random.randint(40, 60)
        genotypes[:, i][random.sample(range(100), k=nb_hetero)] = 1
        all_nb_hetero[(f"f{i}", f"s{i}")] = nb_hetero

    # Generating a FAM file
    fam = pd.DataFrame({
        "fid": [f"f{i}" for i in range(nb_samples)],
        "iid": [f"s{i}" for i in range(nb_samples)],
        "father": [f"father_{i}" for i in range(nb_samples)],
        "mother": [f"mother_{i}" for i in range(nb_samples)],
        "gender": [random.randint(1, 2) for _ in range(nb_samples)],
        "status": [-9] * nb_samples,
    }).set_index(["fid", "iid"])

    # Executing the function
    sex_check.write_heterozygosity(dummy_file, fam=fam, genotypes=genotypes)

    # Checking the file
    assert dummy_file.is_file()
    with open(dummy_file) as f:
        for i, line in enumerate(f):
            row = line.rstrip().split("\t")
            assert len(row) == 4

            # Checking the header
            if i == 0:
                assert row == ["FID", "IID", "PEDSEX", "HETERO"]
                continue

            sample = tuple(row[:2])
            assert int(row[2]) == fam.loc[sample, "gender"]
            assert float(row[3]) == all_nb_hetero[sample] / genotypes.shape[0]


def test_write_heterozygosity_only_one(tmp_path: Path):
    """Test the 'write_heterozygosity' function (only one variant).

    This is for edge cases detection...

    """
    dummy_file = tmp_path / "dummy.txt"

    # Generating genotypes (only homozygous, i.e. 0 and 2)
    nb_samples = 30
    genotypes = 2 * np.random.randint(0, 1, size=(100, nb_samples),
                                      dtype=np.int8)

    # Adding missing values
    all_nb_hetero = {}
    for i in range(genotypes.shape[1]):
        nb_hetero = random.randint(0, 1)
        all_nb_hetero[(f"f{i}", f"s{i}")] = nb_hetero

        if nb_hetero:
            genotypes[0, i] = 1

    # Generating a FAM file
    fam = pd.DataFrame({
        "fid": [f"f{i}" for i in range(nb_samples)],
        "iid": [f"s{i}" for i in range(nb_samples)],
        "father": [f"father_{i}" for i in range(nb_samples)],
        "mother": [f"mother_{i}" for i in range(nb_samples)],
        "gender": [random.randint(1, 2) for _ in range(nb_samples)],
        "status": [-9] * nb_samples,
    }).set_index(["fid", "iid"])

    # Executing the function
    sex_check.write_heterozygosity(dummy_file, fam=fam, genotypes=genotypes)

    # Checking the file
    assert dummy_file.is_file()
    with open(dummy_file) as f:
        for i, line in enumerate(f):
            row = line.rstrip().split("\t")
            assert len(row) == 4

            # Checking the header
            if i == 0:
                assert row == ["FID", "IID", "PEDSEX", "HETERO"]
                continue

            sample = tuple(row[:2])
            assert int(row[2]) == fam.loc[sample, "gender"]
            assert float(row[3]) == all_nb_hetero[sample] / genotypes.shape[0]


def test_write_heterozygosity_no_variants(tmp_path: Path):
    """Test the 'write_heterozygosity' function (no variants).

    This is for edge cases detection...

    """
    dummy_file = tmp_path / "dummy.txt"

    # Generating genotypes for 30 samples
    nb_samples = 30
    genotypes = generate_genotypes(nb_samples=nb_samples, nb_variants=10)

    # Removing variants from the genotypes
    genotypes = genotypes[np.zeros(genotypes.shape[0], dtype=bool), :]

    # Generating a FAM file
    fam = pd.DataFrame({
        "fid": [f"f{i}" for i in range(nb_samples)],
        "iid": [f"s{i}" for i in range(nb_samples)],
        "father": [f"father_{i}" for i in range(nb_samples)],
        "mother": [f"mother_{i}" for i in range(nb_samples)],
        "gender": [random.randint(1, 2) for _ in range(nb_samples)],
        "status": [-9] * nb_samples,
    }).set_index(["fid", "iid"])

    # Executing the function
    sex_check.write_heterozygosity(dummy_file, fam=fam, genotypes=genotypes)

    # Checking the file
    assert dummy_file.is_file()
    with open(dummy_file) as f:
        for i, line in enumerate(f):
            row = line.rstrip().split("\t")
            assert len(row) == 4

            # Checking the header
            if i == 0:
                assert row == ["FID", "IID", "PEDSEX", "HETERO"]
                continue

            sample = tuple(row[:2])
            assert int(row[2]) == fam.loc[sample, "gender"]
            assert row[3] == "-9"


def test_compute_statistics(tmp_path: Path):
    """Test the 'compute_statistics' function.

    We'll generate 100 markers on the X chromosome, and 100 markers on the Y
    chromosome. We'll make sure to only select a subset of sample, to see if
    the subset is well performed by the tool.

    We don't really need to make sure females have low call rate on the Y, or
    that males are homozygous on the X, since we only want to check call rate
    and heterozygosity.

    We'll add no call everywhere, to make sur they are not counted in the
    total.

    """
    # Generating the genotypes (30 samples, 200 variants)
    genotypes = generate_genotypes(nb_variants=200, nb_samples=30)

    # Genotypes on chromosome X should only be homozygous
    genotypes[:100] = 2 * np.random.randint(
        0, 2, size=(100, 30), dtype=np.int8,
    )

    # Generating missing and heterozygous values
    all_nb_missing = {}
    all_hetero_rate = {}
    for i in range(genotypes.shape[1]):
        # Missing values for all variants, but counting only for chromosome Y
        missing_values = random.sample(range(200), k=random.randint(5, 20))
        genotypes[:, i][missing_values] = -1
        all_nb_missing[(f"s{i}", f"s{i}")] = len([
            index for index in missing_values if index >= 100
        ])

        # Heterozygous values (only for chromosome X)
        hetero_values = [
            index
            for index in random.sample(range(100), k=random.randint(5, 20))
            if index not in missing_values
        ]
        genotypes[:100, i][hetero_values] = 1
        all_hetero_rate[(f"s{i}", f"s{i}")] = (
            len(hetero_values) /
            (100 - len([index for index in missing_values if index < 100]))
        )

    # Adding autosomes with some missing values
    autosomes = generate_genotypes(nb_variants=100, nb_samples=30)
    for i in range(autosomes.shape[1]):
        missing_values = random.sample(range(100), k=random.randint(5, 20))
        autosomes[:, i][missing_values] = -1

    # Writing the Plink files
    dummy_bfile = str(tmp_path / "dummy")
    generate_plink_files(
        dummy_bfile,
        genotypes=np.vstack([autosomes, genotypes]),
        chromosomes=np.array(
            sorted(random.choices(range(1, 23), k=100))
            + ([23] * 100)
            + ([24] * 100),
        ),
    )

    # Executing the function without samples (files should be empty)
    sex_check.compute_statistics(
        bfile=dummy_bfile,
        samples=set(),
        prefix=str(tmp_path / "output"),
    )

    # Comparing the hetero file
    with open(tmp_path / "output.chr23.hetero.tsv") as f:
        assert f.read() == "FID\tIID\tPEDSEX\tHETERO\n"

    # Comparing the missing file
    with open(tmp_path / "output.chr24.no_call.tsv") as f:
        assert f.read() == "FID\tIID\tPEDSEX\tN_GENO\tN_MISS\tF_MISS\n"

    # We'll select a random number of samples
    selected_samples = {
        (f"s{i}", f"s{i}")
        for i in random.sample(
            range(genotypes.shape[1]), k=random.randint(15, 20),
        )
    }

    # Executing the function without samples (files should be empty)
    sex_check.compute_statistics(
        bfile=dummy_bfile,
        samples=selected_samples,
        prefix=str(tmp_path / "output"),
    )

    # Checking the hetero file
    with open(tmp_path / "output.chr23.hetero.tsv") as f:
        for i, line in enumerate(f):
            row = line.rstrip().split("\t")

            # The header
            if i == 0:
                assert row == ["FID", "IID", "PEDSEX", "HETERO"]
                continue

            fid, iid, _, hetero = row

            assert (fid, iid) in selected_samples
            assert float(hetero) == all_hetero_rate[(fid, iid)]

    # Checking the missing file
    with open(tmp_path / "output.chr24.no_call.tsv") as f:
        for i, line in enumerate(f):
            row = line.rstrip().split("\t")

            # The header
            if i == 0:
                assert row == ["FID", "IID", "PEDSEX", "N_GENO", "N_MISS",
                               "F_MISS"]
                continue

            fid, iid, _, n_geno, n_miss, f_miss = row

            assert (fid, iid) in selected_samples
            assert int(n_geno) == 100
            assert int(n_miss) == all_nb_missing[(fid, iid)]
            assert float(f_miss) == all_nb_missing[(fid, iid)] / 100
