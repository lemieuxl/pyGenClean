"""Test the tasks module (from utils)."""


import errno
import os

import pytest
from pytest_subprocess.fake_process import FakeProcess

from ...error import ProgramError
from ...utils import task


def test_execute_external_command(fp: FakeProcess):
    """Test the 'execute_external_command' function."""
    # Registering a fake command
    fp.register(
        ["unknown_command", "arg1", "arg2"],
        stdout=b"stdout",
        stderr=b"stderr",
        returncode=0,
    )

    res = task.execute_external_command(["unknown_command", "arg1", "arg2"])
    assert res == "stdout"


def test_execute_external_command_no_decode(fp: FakeProcess):
    """Test the 'execute_external_command' function (no decode)."""
    # Registering a fake command
    fp.register(
        ["unknown_command", "arg1", "arg2"],
        stdout=b"stdout",
        stderr=b"stderr",
        returncode=0,
    )

    res = task.execute_external_command(["unknown_command", "arg1", "arg2"],
                                        decode=False)
    assert res == b"stdout"


def test_execute_external_command_file_error(fp: FakeProcess):
    """Test the 'execute_external_command' function (command not found)."""
    strerror = os.strerror(errno.ENOENT)

    def _callback(*args, **kwargs):
        raise FileNotFoundError(
            errno.ENOENT, strerror, "unknown_command",
        )

    # Registering a fake command
    fp.register(
        ["unknown_command", "arg1", "arg2"],
        callback=_callback,
    )

    with pytest.raises(ProgramError) as program_error:
        task.execute_external_command(["unknown_command", "arg1", "arg2"])
    assert program_error.value.message == f"unknown_command: {strerror}"


def test_execute_external_command_error_code(fp: FakeProcess):
    """Test the 'execute_external_command' function."""
    # Registering a fake command
    fp.register(
        ["unknown_command", "arg1", "arg2"],
        stdout=b"stdout",
        stderr=b"stderr",
        returncode=1,
    )

    with pytest.raises(ProgramError) as program_error:
        task.execute_external_command(["unknown_command", "arg1", "arg2"])

    assert program_error.value.message == (
        "Something went wrong:\nstderr\nunknown_command arg1 arg2"
    )


def test_execute_external_commands(fp: FakeProcess):
    """Test the 'execute_external_commands' function."""
    # We'll use 3 different commands
    commands = []
    for i in range(3):
        commands.append([f"unknown_command_{i}", "arg1", "arg2"])
        fp.register(
            commands[-1],
            stdout=f"stdout_{i}",
        )

    # Launcing the commands
    assert task.execute_external_commands(commands) == [
        f"stdout_{i}" for i in range(3)
    ]
