"""The main command line interface for pyGenClean."""


import sys
import logging
import argparse
from os import path
from shlex import quote

from .error import ProgramError
from .version import pygenclean_version as __version__

from .plate_bias import plate_bias
from .sex_check import sex_check, intensity_plot, baf_lrr_plot


MODULE_MAIN = {
    sex_check.SCRIPT_NAME: {
        "main": sex_check.main,
        intensity_plot.SCRIPT_NAME: intensity_plot.main,
        baf_lrr_plot.SCRIPT_NAME: baf_lrr_plot.main
    }
}


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


def parse_args():
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

    # The sex-check module
    add_sex_check_args(subparser)
    add_plate_bias_args(subparser)

    return parser.parse_args()


def add_sex_check_args(main_subparser):
    """Parser for the sex_check utility."""
    parser = main_subparser.add_parser(
        "sex-check", description=sex_check.DESCRIPTION,
    )

    parser.add_argument(
        "-v", "--version", action="version",
        version=f"pyGenClean sex-check {__version__}",
    )

    subparsers = parser.add_subparsers(
        dest="subcommand", required=True,
        description="Below is a list of tools in the sex-check module. Note "
                    "that 'run' execute the main sex-check pipeline.",
    )

    # Main sex-check
    subparser = subparsers.add_parser(
        "run", description=sex_check.DESCRIPTION, help=sex_check.DESCRIPTION,
    )
    sex_check.add_args(subparser)
    subparser.set_defaults(func=sex_check.main)

    # Intensity plot
    subparser = subparsers.add_parser(
        intensity_plot.SCRIPT_NAME, description=intensity_plot.DESCRIPTION,
        help=intensity_plot.DESCRIPTION,
    )
    intensity_plot.add_args(subparser)
    subparser.set_defaults(func=intensity_plot.main)

    # BAF & LRR plot
    subparser = subparsers.add_parser(
        baf_lrr_plot.SCRIPT_NAME, description=baf_lrr_plot.DESCRIPTION,
        help=baf_lrr_plot.DESCRIPTION,
    )
    baf_lrr_plot.add_args(subparser)
    subparser.set_defaults(func=baf_lrr_plot.main)


def add_plate_bias_args(main_subparser):
    """Parser for the plate_bias utility."""
    parser = main_subparser.add_parser(
        "plate-bias", description=plate_bias.DESCRIPTION,
    )

    parser.add_argument(
        "-v", "--version", action="version",
        version=f"pyGenClean plate-bias {__version__}",
    )

    subparsers = parser.add_subparsers(
        dest="subcommand", required=True,
        description="Below is a list of tools in the plate-bias module. Note "
                    "that 'run' execute the main plate-bias pipeline.",
    )

    # Main plate-bias
    subparser = subparsers.add_parser(
        "run", description=plate_bias.DESCRIPTION, help=plate_bias.DESCRIPTION,
    )
    plate_bias.add_args(subparser)
    subparser.set_defaults(func=plate_bias.main)


def configure_logging(args):
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
