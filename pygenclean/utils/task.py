"""Utility functions to execute tasks."""


import shlex
import logging
import subprocess
import multiprocessing
from typing import List

from ..error import ProgramError


__all__ = ["execute_external_command", "execute_external_commands"]


logger = logging.getLogger(__name__)


def execute_external_command(command: List[str]) -> str:
    """Executes an external command.

    Args:
        command (list): the list containing the command and arguments/options.

    """

    # Fail-proofing the command
    command = [shlex.quote(part) for part in command]
    logger.debug("Executing %s", " ".join(command))

    try:
        proc = subprocess.Popen(
            command, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        )
        outs, errs = proc.communicate()

    except FileNotFoundError as exception:
        raise ProgramError(exception.strerror) from exception

    if proc.returncode != 0:
        # Something went wrong
        command = " ".join(command)
        errs = errs.decode()
        raise ProgramError(f"Something went wrong:\n{errs}\n{command}")

    return outs


def execute_external_commands(commands: List[List[str]],
                              threads: int = 2) -> None:
    """Execute multiple commands using a worker pool.

    Args:
        commands (list): the list of command to execute.
        threads (int): the number of threads.

    """
    with multiprocessing.Pool(processes=threads) as pool:
        results = []

        # Fail-proofing the command
        for command in commands:
            command = [shlex.quote(part) for part in command]
            results.append(pool.apply_async(
                execute_external_command, (command, ),
            ))

        for res in results:
            res.get()
