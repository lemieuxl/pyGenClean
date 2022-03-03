"""The main command line interface for pyGenClean."""


import sys
import logging
import argparse
from os import path
from shlex import quote

from .qc_module import qc_modules, qc_sub_modules
from .error import ProgramError
from .version import pygenclean_version as __version__


def main():  # pylint: disable=missing-docstring
    args = parse_args()

    configure_logging(args)
    logger = logging.getLogger("pyGenClean")

    # Logging useful information
    logger.info("This is pyGenClean version %s", __version__)
    logger.info(
        "Command used: %s %s",
        path.basename(sys.argv[0]),
        " ".join(quote(value) for value in sys.argv[1:]),
    )

    # Dispatching to the correct function
    try:
        args.func(args=args)

    except KeyboardInterrupt:
        logger.info("Cancelled by user")
        sys.exit(0)

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

    # Adding the different parsers for all the QC modules
    add_qc_module_parsers(subparser)

    return parser.parse_args()


def add_qc_module_parsers(main_subparser: argparse._SubParsersAction) -> None:
    """Automatically add all QC module's parsers."""
    for qc_module_name, qc_module in qc_modules.items():
        # The QC module parser
        parser = main_subparser.add_parser(
            qc_module.SCRIPT_NAME, description=qc_module.DESCRIPTION,
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


def configure_logging(args: argparse.Namespace) -> None:
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

    file_handler = logging.FileHandler("pyGenClean.log", mode="w")
    file_handler.setFormatter(logging.Formatter(
        fmt="[%(asctime)s %(name)s %(levelname)s] %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    ))
    logging.root.addHandler(file_handler)


if __name__ == "__main__":
    main()
