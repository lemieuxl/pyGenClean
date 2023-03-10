"""The main pipeline."""


import argparse
import json
import logging
import shlex
import time
from datetime import datetime
from os import path
from pathlib import Path
from typing import Any, Dict, List, Optional

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
    qc_dir = create_qc_dir(args.qc_dir)

    # Reading the configuration file
    conf = read_configuration(args.conf)
    logger.debug("Configuration:\n%s", json.dumps(conf, indent=2))

    # Executing the pipeline
    usable_files = {
        0: {
            "usable_bfile": args.bfile,
        },
    }
    previous_step = 0
    for step, step_info in sorted(conf["steps"].items(),
                                  key=lambda x: int(x[0])):
        # The "pretty" name
        qc_module_name = step_info["module"].replace("-", "_")

        # The sub directory
        sub_dir = qc_dir / f"{step}_{qc_module_name}"
        sub_dir.mkdir()

        # The QC module
        qc_module = qc_modules[qc_module_name]

        # The bfile and out file
        argv = [
            ["--bfile", usable_files[previous_step]["usable_bfile"]],
            ["--out", str(sub_dir / qc_module.DEFAULT_OUT)],
        ]

        # Adding special cases
        if "from_step" in step_info:
            for from_step, from_step_info in step_info["from_step"].items():
                for key, value in from_step_info.items():
                    argv.append([f"--{key}", usable_files[from_step][value]])

        # Building the argument list from the configuration (excluding special
        # cases)
        argv += list(generate_options(step_info))

        # Logging the options used
        logger.info(
            "Executing step %s: %s (%s)",
            step,
            qc_module_name,
            qc_module.DESCRIPTION.rstrip("."),
        )
        for argument in argv:
            logger.info("  %s", " ".join(map(shlex.quote, argument)))

        # Executing the QC module
        usable_files[step] = qc_module.main(
            argv=[item for subitems in argv for item in subitems],
        )

        # The previous step
        previous_step = step


def generate_options(step_info: Dict[str, Any]):
    """Generate the options for ARGV."""
    for key, value in step_info.items():
        if key in ("module", "from_step"):
            continue
        if isinstance(value, bool) and value:
            yield [f"--{key}"]
        else:
            yield [f"--{key}", str(value)]


def create_qc_dir(qc_dir: Optional[str]) -> Path:
    """Create the QC directory."""
    if not qc_dir:
        qc_dir = Path(
            "data_clean_up." + datetime.today().strftime("%Y-%m-%d_%H.%M.%S"),
        )
        while qc_dir.is_dir():
            time.sleep(1)
            qc_dir = Path(
                "data_clean_up."
                + datetime.today().strftime("%Y-%m-%d_%H.%M.%S"),
            )
        qc_dir.mkdir()

    else:
        qc_dir = Path(qc_dir)
        if qc_dir.is_dir():
            raise ValueError(f"{qc_dir}: directory already exists")
        qc_dir.mkdir()

    logger.info("QC will be in '%s'", qc_dir)
    return qc_dir


def read_configuration(fn: str) -> None:
    """Reads the configuration file."""
    # Reading the configuration
    conf = None
    with open(fn, "rb") as f:
        conf = tomllib.load(f)

    # Verifying the configuration
    # TODO: make sure that minimal step is 1 and unique sequential numbers
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

    # The OUTPUT directory
    group = parser.add_argument_group("Output Options")
    group.add_argument(
        "--qc-dir", metavar="DIR", type=str,
        help="The output directory. If not is specified, the directory name "
             "will automatically be generated.",
    )


if __name__ == "__main__":
    main()
