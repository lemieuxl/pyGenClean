"""Utility functions to execute tasks."""


import logging
import multiprocessing
import shlex
import subprocess
from typing import List

from ..error import ProgramError


__all__ = ["execute_external_command", "execute_external_commands"]


logger = logging.getLogger(__name__)


def execute_external_command(command: List[str]) -> str:
    """Executes an external command.

    Args:
        command (list): the list containing the command and arguments/options.

    Returns:
        str: The standard output of the command.

    """

    # Fail-proofing the command
    logger.debug("Executing %s", " ".join(map(shlex.quote, command)))

    try:
        with subprocess.Popen(
            command, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        ) as proc:
            outs, errs = proc.communicate()

    except FileNotFoundError as error:
        raise ProgramError(f"{command[0]}: {error.strerror}") from error

    if proc.returncode != 0:
        # Something went wrong
        command_str = " ".join(command)
        errs_str = errs.decode()
        raise ProgramError(f"Something went wrong:\n{errs_str}\n{command_str}")

    return outs.decode()


def execute_external_commands(commands: List[List[str]],
                              threads: int = 2) -> List[str]:
    """Execute multiple commands using a worker pool.

    Args:
        commands (list): the list of command to execute.
        threads (int): the number of threads.

    Returns:
        list: A list containing the standard output of all the commands.

    """
    with multiprocessing.Pool(processes=threads) as pool:
        results = []

        # Fail-proofing the command
        for command in commands:
            results.append(pool.apply_async(
                execute_external_command, (command, ),
            ))

        return [res.get() for res in results]
