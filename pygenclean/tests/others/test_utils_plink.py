"""Tests the plink sub module from the utils module."""


import random
from datetime import datetime
from pathlib import Path

import pytest
from pytest_mock import MockerFixture

from ...utils import decode_sex, plink


def test_get_markers_on_chrom(tmp_path: Path):
    """Test the 'get_markers_on_chrom' funciton."""
    # We want a big variaty of chromosomes
    chromosomes = random.choices(range(1, 27), k=100)
    chromosomes.sort()

    # The choosen chromosomes
    choosen_chromosomes = set(random.sample(list(range(1, 27)), k=3))

    # Generating the BIM file
    bimfile = tmp_path / "temp.bim"
    expected_variants = set()
    with open(bimfile, "w") as f:
        for i, chrom in enumerate(chromosomes):
            var_name = f"var{i + 1}"
            print(chrom, var_name, 0, i + 101, "A", "B", sep="\t", file=f)

            if chrom in choosen_chromosomes:
                expected_variants.add(var_name)

    # Testing the function
    assert expected_variants == plink.get_markers_on_chrom(
        bimfile, choosen_chromosomes
    )


def test_get_markers_on_chrom_none(tmp_path: Path):
    """Test the 'get_markers_on_chrom' funciton when there are none."""
    # The choosen chromosomes
    choosen_chromosomes = set(random.sample(list(range(1, 27)), k=3))

    # We want a big variaty of chromosomes, except for choosen_chromosomes
    chromosomes = [
        chrom
        for chrom in random.choices(range(1, 27), k=100)
        if chrom not in choosen_chromosomes
    ]
    chromosomes.sort()

    # Generating the BIM file
    bimfile = tmp_path / "temp.bim"
    with open(bimfile, "w") as f:
        for i, chrom in enumerate(chromosomes):
            var_name = f"var{i + 1}"
            print(chrom, var_name, 0, i + 101, "A", "B", sep="\t", file=f)

    # Testing the function
    assert set() == plink.get_markers_on_chrom(
        bimfile, choosen_chromosomes
    )


def test_get_sample_sexes(tmp_path: Path):
    """Test the 'get_sample sexes' function."""
    expected_sex = {
        (f"fid{i + 1}", f"iid{i + 1}"): random.randint(0, 2)
        for i in range(100)
    }

    famfile = tmp_path / "test.fam"
    with open(famfile, "w") as f:
        for i in range(100):
            print(
                f"fid{i + 1}",
                f"iid{i + 1}",
                0,
                0,
                expected_sex[(f"fid{i + 1}", f"iid{i + 1}")],
                -9,
                file=f,
            )

    # Reformatting (the decode_sex function should be tested elsewhere)
    expected_sex = {
        iid: decode_sex(str(sex)) for iid, sex in expected_sex.items()
    }

    assert expected_sex == plink.get_sample_sexes(famfile)


def test_get_sample_sexes_iid(tmp_path: Path):
    """Test the 'get_sample sexes' function (iid only)."""
    expected_sex = {f"iid{i + 1}": random.randint(0, 2) for i in range(100)}

    famfile = tmp_path / "test.fam"
    with open(famfile, "w") as f:
        for i in range(100):
            print(f"fid{i + 1}", f"iid{i + 1}", 0, 0,
                  expected_sex[f"iid{i + 1}"], -9, file=f)

    # Reformatting (the decode_sex function should be tested elsewhere)
    expected_sex = {
        iid: decode_sex(str(sex)) for iid, sex in expected_sex.items()
    }

    assert expected_sex == plink.get_sample_sexes(famfile, only_iid=True)


def test_parse_bim(tmp_path: Path):
    """Test the 'parse_bim' function."""
    # The data for the BIM file
    bim_content = []
    for i in range(100):
        bim_content.append((
            random.randint(1, 26),
            f"var{i + 1}",
            random.random() * 100,
            random.randint(1, 1000000),
            *random.sample("ACGT", 2),
        ))

    # Writing the content
    bim_file = tmp_path / "test.bim"
    with open(bim_file, "w") as f:
        for row in bim_content:
            print(*row, sep="\t", file=f)

    # Parsing using the function
    for expected, observed in zip(bim_content, plink.parse_bim(bim_file)):
        chrom, name, cm, pos, a1, a2 = expected

        assert chrom == observed.chrom
        assert name == observed.name
        assert cm == pytest.approx(observed.cm, 0.0000001)
        assert pos == observed.pos
        assert a1 == observed.a1
        assert a2 == observed.a2


def test_parse_fam(tmp_path: Path):
    """Test the 'parse_fam' function."""
    # The data for the FAM file
    fam_content = []
    for i in range(100):
        fam_content.append((
            f"family_{i + 1}",
            f"sample_{i + 1}",
            f"father_{i + 1}",
            f"mother_{i + 1}",
            random.randint(0, 2),
            str(random.randint(1, 100)),
        ))

    # Writing the content
    fam_file = tmp_path / "test.fam"
    with open(fam_file, "w") as f:
        for row in fam_content:
            print(*row, sep=" ", file=f)

    # Parsing using the function
    for expected, observed in zip(fam_content, plink.parse_fam(fam_file)):
        fid, iid, father, mother, sex, status = expected

        # Encoding sex
        if sex == 1:
            sex = "Male"
        elif sex == 2:
            sex = "Female"
        else:
            sex = "Unknown"

        assert fid == observed.fid
        assert iid == observed.iid
        assert father == observed.father
        assert mother == observed.mother
        assert sex == observed.sex
        assert status == observed.status


def test_check_files_ok(tmp_path: Path):
    """Test the 'parse_fam' file (ok)."""
    for extension in ("bed", "bim", "fam"):
        with open(tmp_path / f"dummy.{extension}", "w") as _:
            pass

    assert plink.check_files(tmp_path / "dummy")


def test_check_files_not_ok(tmp_path: Path):
    """Test the 'parse_fam' file (not ok)."""
    for extension in ("bed", "bim", "fam"):
        with open(tmp_path / f"dummy.{extension}", "w") as _:
            pass

    # Deleting one file at a time, and checking the function returns False
    for extension in ("bed", "bim", "fam"):
        deleted_file = tmp_path / f"dummy.{extension}"
        deleted_file.unlink()
        assert not plink.check_files(tmp_path / "dummy")
        with open(deleted_file, "w") as _:
            pass


def test_check_unique_iid_true(tmp_path: Path):
    """Test the 'check_unique_iid' function (only unique ID)."""
    fam_file = tmp_path / "dummy.fam"
    with open(fam_file, "w") as f:
        for i in range(100):
            print(f"s{i + 1}", f"s{i + 1}", f"f{i + 1}", f"m{i + 1}",
                  random.randint(1, 2), -9, file=f)

    assert plink.check_unique_iid(str(fam_file.with_suffix("")))


def test_check_unique_iid_false(tmp_path: Path):
    """Test the 'check_unique_iid' function (not unique ID)."""
    fam_file = tmp_path / "dummy.fam"
    with open(fam_file, "w") as f:
        for i in range(10):
            for j in range(30):
                print(f"f{i + 1}", f"s{j + 1}", f"f{j + 1}", f"m{j + 1}",
                      random.randint(1, 2), -9, file=f)

    assert not plink.check_unique_iid(str(fam_file.with_suffix("")))


def test_split_lines():
    """Test the 'split_lines' function."""
    line = "   123 123 121  325 2\t\t\t1 23  1   3124345\t32413 4532\t\t 23  4"
    assert ["123", "123", "121", "325", "2", "1", "23", "1", "3124345",
            "32413", "4532", "23", "4"] == plink.split_line(line)


# Calling the function with all possible combination of subset
#   - keep alone
#   - remove alone
#   - exclude alone
#   - extract alone
#   - keep with exclude
#   - keep with extract
#   - remove with exclude
#   - remove with extract


def test_subset_keep_alone(mocker: MockerFixture):
    """Test the 'subset' function (keep alone)."""
    mocked = mocker.patch("pygenclean.utils.plink.execute_external_command")

    plink.subset(
        bfile="dummy",
        out="out",
        samples="samples_keep.txt",
        sample_subset_type="keep",
    )

    # Checking
    mocked.assert_called_once_with([
        "plink1.9",
        "--noweb",
        "--bfile", "dummy",
        "--make-bed",
        "--out", "out",
        "--keep", "samples_keep.txt",
    ])


def test_subset_remove_alone(mocker: MockerFixture):
    """Test the 'subset' function (remove alone)."""
    mocked = mocker.patch("pygenclean.utils.plink.execute_external_command")

    plink.subset(
        bfile="dummy",
        out="out",
        samples="samples_remove.txt",
        sample_subset_type="remove",
    )

    # Checking
    mocked.assert_called_once_with([
        "plink1.9",
        "--noweb",
        "--bfile", "dummy",
        "--make-bed",
        "--out", "out",
        "--remove", "samples_remove.txt",
    ])


def test_subset_exclude_alone(mocker: MockerFixture):
    """Test the 'subset' function (exclude alone)."""
    mocked = mocker.patch("pygenclean.utils.plink.execute_external_command")

    plink.subset(
        bfile="dummy",
        out="out",
        markers="markers_exclude.txt",
        marker_subset_type="exclude",
    )

    # Checking
    mocked.assert_called_once_with([
        "plink1.9",
        "--noweb",
        "--bfile", "dummy",
        "--make-bed",
        "--out", "out",
        "--exclude", "markers_exclude.txt",
    ])


def test_subset_extract_alone(mocker: MockerFixture):
    """Test the 'subset' function (extract alone)."""
    mocked = mocker.patch("pygenclean.utils.plink.execute_external_command")

    plink.subset(
        bfile="dummy",
        out="out",
        markers="markers_extract.txt",
        marker_subset_type="extract",
    )

    # Checking
    mocked.assert_called_once_with([
        "plink1.9",
        "--noweb",
        "--bfile", "dummy",
        "--make-bed",
        "--out", "out",
        "--extract", "markers_extract.txt",
    ])


def test_subset_keep_exclude(mocker: MockerFixture):
    """Test the 'subset' function (keep and exclude)."""
    mocked = mocker.patch("pygenclean.utils.plink.execute_external_command")

    plink.subset(
        bfile="dummy",
        out="out",
        samples="samples_keep.txt",
        sample_subset_type="keep",
        markers="markers_exclude.txt",
        marker_subset_type="exclude",
    )

    # Checking
    mocked.assert_called_once_with([
        "plink1.9",
        "--noweb",
        "--bfile", "dummy",
        "--make-bed",
        "--out", "out",
        "--keep", "samples_keep.txt",
        "--exclude", "markers_exclude.txt",
    ])


def test_subset_keep_extract(mocker: MockerFixture):
    """Test the 'subset' function (keep and extract)."""
    mocked = mocker.patch("pygenclean.utils.plink.execute_external_command")

    plink.subset(
        bfile="dummy",
        out="out",
        samples="samples_keep.txt",
        sample_subset_type="keep",
        markers="markers_extract.txt",
        marker_subset_type="extract",
    )

    # Checking
    mocked.assert_called_once_with([
        "plink1.9",
        "--noweb",
        "--bfile", "dummy",
        "--make-bed",
        "--out", "out",
        "--keep", "samples_keep.txt",
        "--extract", "markers_extract.txt",
    ])


def test_subset_remove_exclude(mocker: MockerFixture):
    """Test the 'subset' function (remove and exclude)."""
    mocked = mocker.patch("pygenclean.utils.plink.execute_external_command")

    plink.subset(
        bfile="dummy",
        out="out",
        samples="samples_remove.txt",
        sample_subset_type="remove",
        markers="markers_exclude.txt",
        marker_subset_type="exclude",
    )

    # Checking
    mocked.assert_called_once_with([
        "plink1.9",
        "--noweb",
        "--bfile", "dummy",
        "--make-bed",
        "--out", "out",
        "--remove", "samples_remove.txt",
        "--exclude", "markers_exclude.txt",
    ])


def test_subset_remove_extract(mocker: MockerFixture):
    """Test the 'subset' function (remove and extract)."""
    mocked = mocker.patch("pygenclean.utils.plink.execute_external_command")

    plink.subset(
        bfile="dummy",
        out="out",
        samples="samples_remove.txt",
        sample_subset_type="remove",
        markers="markers_extract.txt",
        marker_subset_type="extract",
    )

    # Checking
    mocked.assert_called_once_with([
        "plink1.9",
        "--noweb",
        "--bfile", "dummy",
        "--make-bed",
        "--out", "out",
        "--remove", "samples_remove.txt",
        "--extract", "markers_extract.txt",
    ])


def test_subset_invalid_sample():
    """Test the 'subset' function (invalid sample subset type)."""
    # Remove alone
    with pytest.raises(ValueError) as value_error:
        plink.subset(
            bfile="dummy",
            out="out",
            samples="samples_remove.txt",
            sample_subset_type="removed",
        )

    assert str(value_error.value) == "removed: invalid sample subset"


def test_subset_invalid_marker():
    """Test the 'subset' function (invalid marker subset type)."""
    # Remove alone
    with pytest.raises(ValueError) as value_error:
        plink.subset(
            bfile="dummy",
            out="out",
            markers="markers_extract.txt",
            marker_subset_type="extracted",
        )

    assert str(value_error.value) == "extracted: invalid marker subset"


def test_subset_invalid():
    """Test the 'subset' function (invalid marker subset type)."""
    # Remove alone
    with pytest.raises(ValueError) as value_error:
        plink.subset(
            bfile="dummy",
            out="out",
        )

    assert str(value_error.value) == "no subset was selected"


def test_subset_keep_alone_old_plink(mocker: MockerFixture):
    """Test the 'subset' function (keep alone, old Plink)."""
    mocked = mocker.patch("pygenclean.utils.plink.execute_external_command")

    plink.subset(
        bfile="dummy",
        out="out",
        samples="samples_keep.txt",
        sample_subset_type="keep",
        use_original_plink=True
    )

    # Checking
    mocked.assert_called_once_with([
        "plink",
        "--noweb",
        "--bfile", "dummy",
        "--make-bed",
        "--out", "out",
        "--keep", "samples_keep.txt",
    ])


def test_subset_remove_alone_old_plink(mocker: MockerFixture):
    """Test the 'subset' function (remove alone, old Plink)."""
    mocked = mocker.patch("pygenclean.utils.plink.execute_external_command")

    plink.subset(
        bfile="dummy",
        out="out",
        samples="samples_remove.txt",
        sample_subset_type="remove",
        use_original_plink=True,
    )

    # Checking
    mocked.assert_called_once_with([
        "plink",
        "--noweb",
        "--bfile", "dummy",
        "--make-bed",
        "--out", "out",
        "--remove", "samples_remove.txt",
    ])


def test_subset_exclude_alone_old_plink(mocker: MockerFixture):
    """Test the 'subset' function (exclude alone, old Plink)."""
    mocked = mocker.patch("pygenclean.utils.plink.execute_external_command")

    plink.subset(
        bfile="dummy",
        out="out",
        markers="markers_exclude.txt",
        marker_subset_type="exclude",
        use_original_plink=True,
    )

    # Checking
    mocked.assert_called_once_with([
        "plink",
        "--noweb",
        "--bfile", "dummy",
        "--make-bed",
        "--out", "out",
        "--exclude", "markers_exclude.txt",
    ])


def test_subset_extract_alone_old_plink(mocker: MockerFixture):
    """Test the 'subset' function (extract alone, old Plink)."""
    mocked = mocker.patch("pygenclean.utils.plink.execute_external_command")

    plink.subset(
        bfile="dummy",
        out="out",
        markers="markers_extract.txt",
        marker_subset_type="extract",
        use_original_plink=True,
    )

    # Checking
    mocked.assert_called_once_with([
        "plink",
        "--noweb",
        "--bfile", "dummy",
        "--make-bed",
        "--out", "out",
        "--extract", "markers_extract.txt",
    ])


def test_subset_keep_exclude_old_plink(mocker: MockerFixture):
    """Test the 'subset' function (keep and exclude, old Plink)."""
    mocked = mocker.patch("pygenclean.utils.plink.execute_external_command")

    plink.subset(
        bfile="dummy",
        out="out",
        samples="samples_keep.txt",
        sample_subset_type="keep",
        markers="markers_exclude.txt",
        marker_subset_type="exclude",
        use_original_plink=True,
    )

    # Checking
    mocked.assert_called_once_with([
        "plink",
        "--noweb",
        "--bfile", "dummy",
        "--make-bed",
        "--out", "out",
        "--keep", "samples_keep.txt",
        "--exclude", "markers_exclude.txt",
    ])


def test_subset_keep_extract_old_plink(mocker: MockerFixture):
    """Test the 'subset' function (keep and extract, old Plink)."""
    mocked = mocker.patch("pygenclean.utils.plink.execute_external_command")

    plink.subset(
        bfile="dummy",
        out="out",
        samples="samples_keep.txt",
        sample_subset_type="keep",
        markers="markers_extract.txt",
        marker_subset_type="extract",
        use_original_plink=True,
    )

    # Checking
    mocked.assert_called_once_with([
        "plink",
        "--noweb",
        "--bfile", "dummy",
        "--make-bed",
        "--out", "out",
        "--keep", "samples_keep.txt",
        "--extract", "markers_extract.txt",
    ])


def test_subset_remove_exclude_old_plink(mocker: MockerFixture):
    """Test the 'subset' function (remove and exclude, old Plink)."""
    mocked = mocker.patch("pygenclean.utils.plink.execute_external_command")

    plink.subset(
        bfile="dummy",
        out="out",
        samples="samples_remove.txt",
        sample_subset_type="remove",
        markers="markers_exclude.txt",
        marker_subset_type="exclude",
        use_original_plink=True,
    )

    # Checking
    mocked.assert_called_once_with([
        "plink",
        "--noweb",
        "--bfile", "dummy",
        "--make-bed",
        "--out", "out",
        "--remove", "samples_remove.txt",
        "--exclude", "markers_exclude.txt",
    ])


def test_subset_remove_extract_old_plink(mocker: MockerFixture):
    """Test the 'subset' function (remove and extract, old Plink)."""
    mocked = mocker.patch("pygenclean.utils.plink.execute_external_command")

    plink.subset(
        bfile="dummy",
        out="out",
        samples="samples_remove.txt",
        sample_subset_type="remove",
        markers="markers_extract.txt",
        marker_subset_type="extract",
        use_original_plink=True,
    )

    # Checking
    mocked.assert_called_once_with([
        "plink",
        "--noweb",
        "--bfile", "dummy",
        "--make-bed",
        "--out", "out",
        "--remove", "samples_remove.txt",
        "--extract", "markers_extract.txt",
    ])


def test_subset_invalid_sample_old_plink():
    """Test the 'subset' function (invalid sample subset type, old Plink)."""
    # Remove alone
    with pytest.raises(ValueError) as value_error:
        plink.subset(
            bfile="dummy",
            out="out",
            samples="samples_remove.txt",
            sample_subset_type="removed",
            use_original_plink=True
        )

    assert str(value_error.value) == "removed: invalid sample subset"


def test_subset_invalid_marker_old_plink():
    """Test the 'subset' function (invalid marker subset type, old Plink)."""
    # Remove alone
    with pytest.raises(ValueError) as value_error:
        plink.subset(
            bfile="dummy",
            out="out",
            markers="markers_extract.txt",
            marker_subset_type="extracted",
            use_original_plink=True,
        )

    assert str(value_error.value) == "extracted: invalid marker subset"


def test_subset_invalid_old_plink():
    """Test the 'subset' function (invalid marker subset type, old Plink)."""
    # Remove alone
    with pytest.raises(ValueError) as value_error:
        plink.subset(
            bfile="dummy",
            out="out",
            use_original_plink=True,
        )

    assert str(value_error.value) == "no subset was selected"


def test_rename_markers(tmp_path: Path):
    """Test the 'rename_markers' function."""
    # Creating a "rename" file and a dummy bim file
    rename = tmp_path / "rename.txt"
    bim = tmp_path / "dummy.bim"

    to_rename = {}
    with open(rename, "w") as rename_f, open(bim, "w") as bim_f:
        for i in range(100):
            var_name = f"var{i + 1}"

            # Writing the BIM (and saving a copy)
            print(
                random.choice(range(1, 27)),
                var_name,
                random.random() * 100,
                random.randint(1, 10000),
                *random.sample("ACGT", 2),
                sep="\t",
                file=bim_f,
            )

            # Renaming 50% of the variants
            if random.random() < 0.5:
                to_rename[var_name] = "new_" + var_name
                print(var_name, to_rename[var_name], sep="\t", file=rename_f)

    # The content of the original BIM
    with open(bim, "rb") as f:
        ori_content = f.read()

    # Executing the function
    plink.rename_markers(str(bim.with_suffix("")), rename)

    # Finding the current file and comparing it with the original content
    now = datetime.now().strftime("%Y-%m-%d_%H.%M")
    backup_files = list(tmp_path.glob(f"{bim.name}.{now}.*"))
    assert len(backup_files) == 1
    with open(backup_files[0], "rb") as f:
        assert f.read() == ori_content

    # Making sure the names changed
    with open(bim) as new_f, open(backup_files[0]) as ori_f:
        new_bim = new_f.read().splitlines()
        ori_bim = ori_f.read().splitlines()

        # Checking length
        assert len(new_bim) == len(ori_bim)

        # checking content
        for new_row, ori_row in zip(new_bim, ori_bim):
            ori_row = ori_row.split()
            ori_row[1] = to_rename.get(ori_row[1], ori_row[1])
            assert new_row.split() == ori_row


def test_compute_freq(mocker: MockerFixture):
    """Test the 'compute_freq' function."""
    mocked = mocker.patch("pygenclean.utils.plink.execute_external_command")

    plink.compute_freq(
        bfile="dummy_prefix",
        out="out",
    )

    mocked.assert_called_once_with([
        "plink1.9",
        "--noweb",
        "--bfile", "dummy_prefix",
        "--freq",
        "--out", "out",
    ])


def test_compute_freq_with_extract(mocker: MockerFixture):
    """Test the 'compute_freq' function (with variant extraction)."""
    mocked = mocker.patch("pygenclean.utils.plink.execute_external_command")

    plink.compute_freq(
        bfile="dummy_prefix",
        out="out",
        extract="extract_file",
    )

    mocked.assert_called_once_with([
        "plink1.9",
        "--noweb",
        "--bfile", "dummy_prefix",
        "--freq",
        "--out", "out",
        "--extract", "extract_file",
    ])


def test_compute_freq_old_plink(mocker: MockerFixture):
    """Test the 'compute_freq' function (old Plink)."""
    mocked = mocker.patch("pygenclean.utils.plink.execute_external_command")

    plink.compute_freq(
        bfile="dummy_prefix",
        out="out",
        use_original_plink=True
    )

    mocked.assert_called_once_with([
        "plink",
        "--noweb",
        "--bfile", "dummy_prefix",
        "--freq",
        "--out", "out",
    ])


def test_compute_freq_with_extract_old_plink(mocker: MockerFixture):
    """Test the 'compute_freq' function (variant extraction and old Plink)."""
    mocked = mocker.patch("pygenclean.utils.plink.execute_external_command")

    plink.compute_freq(
        bfile="dummy_prefix",
        out="out",
        extract="extract_file",
        use_original_plink=True,
    )

    mocked.assert_called_once_with([
        "plink",
        "--noweb",
        "--bfile", "dummy_prefix",
        "--freq",
        "--out", "out",
        "--extract", "extract_file",
    ])


def test_compare_bim(tmp_path: Path):
    """Test the 'compare_bim' function."""
    bim1 = tmp_path / "dummy_1.bim"
    bim2 = tmp_path / "dummy_2.bim"

    expected_in_first = set()
    expected_in_second = set()
    expected_in_both = set()
    with open(bim1, "w") as f1, open(bim2, "w") as f2:
        for i in range(1000):
            row = [
                random.randint(1, 26),
                f"var{i + 1}",
                random.random() * 100,
                random.randint(100, 1000000),
                *random.sample("ACGT", 2),
            ]

            # A random number
            random_nb = random.random()

            # Only in first BIM
            if random_nb < (1/3):
                print(*row, sep="\t", file=f1)
                expected_in_first.add(row[1])

            # Only in second BIM
            elif random_nb < (2/3):
                print(*row, sep="\t", file=f2)
                expected_in_second.add(row[1])

            # In both BIMs
            else:
                print(*row, sep="\t", file=f1)
                print(*row, sep="\t", file=f2)
                expected_in_both.add(row[1])

    expected = (expected_in_first, expected_in_both, expected_in_second)
    assert plink.compare_bim(bim1, bim2) == expected


def test_compare_fam(tmp_path: Path):
    """Test the 'compare_fam' function."""
    fam1 = tmp_path / "dummy_1.fam"
    fam2 = tmp_path / "dummy_2.fam"

    expected_in_first = set()
    expected_in_second = set()
    expected_in_both = set()
    with open(fam1, "w") as f1, open(fam2, "w") as f2:
        for i in range(100):
            for j in range(100):
                row = [
                    f"f{i + 1}",
                    f"s{j + 1}",
                    f"father_{i + 1}",
                    f"mother_{i + 1}",
                    random.randint(0, 2),
                    random.randint(0, 2),
                ]

                # A random number
                random_nb = random.random()

                # Only in first BIM
                if random_nb < (1/3):
                    print(*row, sep="\t", file=f1)
                    expected_in_first.add(tuple(row[:2]))

                # Only in second BIM
                elif random_nb < (2/3):
                    print(*row, sep="\t", file=f2)
                    expected_in_second.add(tuple(row[:2]))

                # In both BIMs
                else:
                    print(*row, sep="\t", file=f1)
                    print(*row, sep="\t", file=f2)
                    expected_in_both.add(tuple(row[:2]))

    expected = (expected_in_first, expected_in_both, expected_in_second)
    assert plink.compare_fam(fam1, fam2) == expected


def test_merge_files(mocker: MockerFixture, tmp_path: Path):
    """Test the 'merge_files' function."""
    mocked = mocker.patch("pygenclean.utils.plink.execute_external_command")
    files_to_merge = tmp_path / "dummy.files_to_merge"

    plink.merge_files(
        prefixes=[f"bfile_{i + 1}" for i in range(10)],
        out=str(files_to_merge.with_suffix("")),
    )

    # Checking the file exists
    assert files_to_merge.is_file()

    # Checking the content of the file
    with open(files_to_merge) as f:
        observed_content = f.read().splitlines()
    assert observed_content == [
        f"bfile_{i + 1}.bed bfile_{i + 1}.bim bfile_{i + 1}.fam"
        for i in range(1, 10)
    ]

    # Checking the command
    mocked.assert_called_once_with([
        "plink1.9",
        "--noweb",
        "--bfile", "bfile_1",
        "--merge-list", str(files_to_merge),
        "--make-bed",
        "--out", str(files_to_merge.with_suffix("")),
    ])


def test_merge_files_old_plink(mocker: MockerFixture, tmp_path: Path):
    """Test the 'merge_files' function (old Plink)."""
    mocked = mocker.patch("pygenclean.utils.plink.execute_external_command")
    files_to_merge = tmp_path / "dummy.files_to_merge"

    plink.merge_files(
        prefixes=[f"bfile_{i + 1}" for i in range(10)],
        out=str(files_to_merge.with_suffix("")),
        use_original_plink=True,
    )

    # Checking the file exists
    assert files_to_merge.is_file()

    # Checking the content of the file
    with open(files_to_merge) as f:
        observed_content = f.read().splitlines()
    assert observed_content == [
        f"bfile_{i + 1}.bed bfile_{i + 1}.bim bfile_{i + 1}.fam"
        for i in range(1, 10)
    ]

    # Checking the command
    mocked.assert_called_once_with([
        "plink",
        "--noweb",
        "--bfile", "bfile_1",
        "--merge-list", str(files_to_merge),
        "--make-bed",
        "--out", str(files_to_merge.with_suffix("")),
    ])


def test_flip_markers(mocker: MockerFixture):
    """Test the 'flip_markers' function."""
    mocked = mocker.patch("pygenclean.utils.plink.execute_external_command")

    plink.flip_markers(
        prefix="dummy_bfile",
        out="dummy_out",
        markers="dummy_markers.txt",
    )

    mocked.assert_called_once_with([
        "plink1.9",
        "--noweb",
        "--bfile", "dummy_bfile",
        "--flip", "dummy_markers.txt",
        "--make-bed",
        "--out", "dummy_out",
    ])


def test_flip_markers_old_plink(mocker: MockerFixture):
    """Test the 'flip_markers' function (old Plink)."""
    mocked = mocker.patch("pygenclean.utils.plink.execute_external_command")

    plink.flip_markers(
        prefix="dummy_bfile",
        out="dummy_out",
        markers="dummy_markers.txt",
        use_original_plink=True,
    )

    mocked.assert_called_once_with([
        "plink",
        "--noweb",
        "--bfile", "dummy_bfile",
        "--flip", "dummy_markers.txt",
        "--make-bed",
        "--out", "dummy_out",
    ])


def test_get_acgt_geno_map():
    """Test the 'get_acgt_geno_map' function."""
    assert plink.get_acgt_geno_map("A", "B") == {
        0:  "B B",
        1:  "A B",
        2:  "A A",
        -1: "0 0"
    }
