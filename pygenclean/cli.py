"""The main command line interface for pyGenClean."""


import argparse
import logging
import sys
import time
from datetime import datetime
from os import path
from shlex import quote

import argcomplete

from .error import ProgramError
from .pipeline import scheduler as main_pipeline
from .qc_modules import qc_modules, qc_sub_modules
from .version import pygenclean_version as __version__


def main():
    """The main command line interface for pyGenClean."""
    args = parse_args()

    timestamp = configure_logging(args)
    logger = logging.getLogger("pyGenClean")

    # Logging useful information
    logger.info("This is pyGenClean version %s", __version__)
    logger.info(
        "Command used: %s %s",
        path.basename(sys.argv[0]),
        " ".join(map(quote, sys.argv[1:])),
    )

    # If we run the pipeline, we need to pass the timestamp (for the QC
    # directory creation)
    if args.command == "pipeline":
        args.qc_dir = "data_clean_up." + timestamp

    # Dispatching to the correct function
    try:
        args.func(args=args)

    except KeyboardInterrupt:
        logger.info("Cancelled by user")
        sys.exit(1)

    except ProgramError as error:
        logger.error(error.message)
        sys.exit(1)


def parse_args() -> argparse.Namespace:
    """Parses the arguments and function."""
    parser = argparse.ArgumentParser(
        description="Runs pyGenClean.",
    )

    parser.add_argument(
        "-v", "--version", action="version",
        version=f"pyGenClean vesion {__version__}",
    )

    parser.add_argument(
        "--debug", action="store_true", help="Enable debug logging.",
    )

    subparser = parser.add_subparsers(dest="command", required=True)

    # Adding the different parsers for the pipeline and all the QC modules
    add_pipeline_parser(subparser)
    add_qc_module_parsers(subparser)

    # Argument completion
    argcomplete.autocomplete(parser)

    return parser.parse_args()


def add_pipeline_parser(main_subparser: argparse._SubParsersAction) -> None:
    """Automatically add the pipeline's parser."""
    parser = main_subparser.add_parser(
        "pipeline", description=main_pipeline.DESCRIPTION,
        help=main_pipeline.DESCRIPTION,
    )

    # Version information
    parser.add_argument(
        "-v", "--version", action="version",
        version=f"pyGenClean pipeline {__version__}",
    )

    main_pipeline.add_args(parser)
    parser.set_defaults(func=main_pipeline.main)


def add_qc_module_parsers(main_subparser: argparse._SubParsersAction) -> None:
    """Automatically add all QC module's parsers."""
    for qc_module_name, qc_module in qc_modules.items():
        # The QC module parser
        parser = main_subparser.add_parser(
            qc_module.SCRIPT_NAME, description=qc_module.DESCRIPTION,
            help=qc_module.DESCRIPTION,
        )

        # Version information
        parser.add_argument(
            "-v", "--version", action="version",
            version=f"pyGenClean {qc_module.SCRIPT_NAME} {__version__}",
        )

        # The QC module sub parser
        subparsers = parser.add_subparsers(
            dest="subcommand", required=True,
            description=f"Below is a list of tools in the "
                        f"{qc_module.SCRIPT_NAME} module. Note that 'run' "
                        f"executes the main {qc_module.SCRIPT_NAME} pipeline.",
        )

        # Main qc module script
        subparser = subparsers.add_parser(
            "run", description=qc_module.DESCRIPTION,
            help=qc_module.DESCRIPTION,
        )
        qc_module.add_args(subparser)
        subparser.set_defaults(func=qc_module.main)

        # Adding sub qc module parsers, if any
        if qc_module_name in qc_sub_modules:
            for qc_sub_module in qc_sub_modules[qc_module_name].values():
                subparser = subparsers.add_parser(
                    qc_sub_module.SCRIPT_NAME,
                    description=qc_sub_module.DESCRIPTION,
                    help=qc_sub_module.DESCRIPTION,
                )
                qc_sub_module.add_args(subparser)
                subparser.set_defaults(func=qc_sub_module.main)


def configure_logging(args: argparse.Namespace) -> str:
    """Configures the logging."""
    # We want INFO by default, unless specified otherwise
    logging_level = logging.INFO
    if args.debug:
        logging_level = logging.DEBUG

    # Configuring the logging
    logging.basicConfig(
        level=logging_level,
        format="[%(asctime)s %(name)s %(levelname)s] %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )

    # Finding the name of the log file
    timestamp = datetime.today().strftime("%Y-%m-%d_%H.%M.%S")
    filename = f"pyGenClean.{timestamp}.log"
    while path.isfile(filename):
        time.sleep(1)
        timestamp = datetime.today().strftime("%Y-%m-%d_%H.%M.%S")
        filename = f"pyGenClean.{timestamp}.log"

    # Adding the file handler
    file_handler = logging.FileHandler(filename, mode="w")
    file_handler.setFormatter(logging.Formatter(
        fmt="[%(asctime)s %(name)s %(levelname)s] %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    ))
    logging.root.addHandler(file_handler)

    return timestamp


if __name__ == "__main__":
    main()
