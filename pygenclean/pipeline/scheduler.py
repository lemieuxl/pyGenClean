"""The main pipeline."""


import argparse
import logging
from pprint import pprint
import time
from datetime import datetime
from os import path
from pathlib import Path
from typing import List, Optional

import tomllib

from ..error import ProgramError
from ..qc_modules import qc_modules
from ..utils import plink as plink_utils
from ..utils import timer
from ..version import pygenclean_version as __version__


DESCRIPTION = "The main pyGenClean pipeline."


logger = logging.getLogger(__name__)


@timer(logger)
def main(args: Optional[argparse.Namespace] = None,
         argv: Optional[List[str]] = None) -> None:
    """The main pipeline.

    Args:
        args (argparse.Namespace): the arguments and options.
        argv (list): the argument as list.

    """
    if args is None:
        args = parse_args(argv)
    check_args(args)

    # Creating the QC directory
    qc_dir = create_qc_dir()

    # Reading the configuration file
    conf = read_configuration(args.conf)
    pprint(conf)

    # Executing the pipeline
    usable_files = {
        "bfile": args.bfile,
    }
    for step, step_info in sorted(conf["steps"].items(),
                                  key=lambda x: int(x[0])):
        logger.info("Executing step %s: %s", step, step_info["module"])

        # The sub directory
        sub_dir = qc_dir / f"{step}_{step_info['module']}"
        sub_dir.mkdir()

        # The QC module
        qc_module = qc_modules[step_info["module"]]

        # The bfile and out file
        argv = [
            ["--bfile", usable_files["bfile"]],
            ["--out", str(sub_dir / qc_module.DEFAULT_OUT)],
        ]

        # Adding special cases
        if "from_previous" in step_info:
            argv += [
                [f"--{key}", usable_files[key]]
                for key in step_info["from_previous"]
            ]

        # Building the argument list from the configuration (excluding special
        # cases)
        argv += [
            [f"--{key}", str(value)]
            for key, value in step_info.items()
            if key not in ("module", "from_previous")
        ]

        # Executing the QC module
        usable_files = qc_module.main(
            argv=[item for subitems in argv for item in subitems],
        )


def create_qc_dir() -> Path:
    """Create the QC directory."""
    qc_dir = Path(
        "data_clean_up." + datetime.today().strftime("%Y-%m-%d_%H.%M.%S"),
    )
    while qc_dir.is_dir():
        time.sleep(1)
        qc_dir = Path(
            "data_clean_up." + datetime.today().strftime("%Y-%m-%d_%H.%M.%S"),
        )
    qc_dir.mkdir()

    return qc_dir


def read_configuration(fn: str) -> None:
    """Reads the configuration file."""
    conf = None
    with open(fn, "rb") as f:
        conf = tomllib.load(f)
    return conf


def check_args(args: argparse.Namespace) -> None:
    """Checks the arguments and options.

    Args:
        args (argparse.Namespace): the arguments and options.

    If there is a problem with an option, an exception is raised using the
    :py:class:`ProgramError` class, a message is printed to the
    :class:`sys.stderr` and the program exists with code 1.

    """
    # Check if we have the bed, bim and fam files
    if not plink_utils.check_files(args.bfile):
        raise ProgramError(f"{args.bfile}: no such binary files")

    # Checking the configuration file
    if not path.isfile(args.conf):
        raise ProgramError(f"{args.conf}: no such file")


def parse_args(argv: Optional[List[str]] = None) -> argparse.Namespace:
    """Parses the arguments and function."""
    parser = argparse.ArgumentParser(description=DESCRIPTION)

    parser.add_argument(
        "-v", "--version", action="version",
        version=f"pyGenClean {__version__}",
    )

    # Adding the arguments and options
    add_args(parser)

    return parser.parse_args(argv)


def add_args(parser: argparse.ArgumentParser) -> None:
    """Add arguments and options to the parser."""
    # The INPUT files
    group = parser.add_argument_group("Input File")
    group.add_argument(
        "--bfile", type=str, metavar="FILE", required=True,
        help="The input file prefix (will find the plink binary files by "
             "appending the prefix to the .bim, .bed and .fam files, "
             "respectively.",
    )
    group.add_argument(
        "--conf", type=str, metavar="TOML", required=True,
        help="The configuration for the data clean up (TOML format).",
    )


if __name__ == "__main__":
    main()
