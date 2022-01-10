"""The main command line interface for pyGenClean."""


import sys
import logging
import argparse
from os import path
from shlex import quote

from .error import ProgramError
from .version import pygenclean_version as __version__

from .plate_bias import plate_bias
from .related_samples import related_samples
from .sex_check import sex_check, intensity_plot, baf_lrr_plot
from .flag_maf import flag_maf
from .sample_call_rate import sample_call_rate
from .marker_call_rate import marker_call_rate
from .flag_hw import flag_hw


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

    # The different modules
    add_sex_check_args(subparser)
    add_plate_bias_args(subparser)
    add_related_samples_args(subparser)
    add_flag_maf(subparser)
    add_sample_call_rate(subparser)
    add_marker_call_rate(subparser)
    add_flag_hw(subparser)

    return parser.parse_args()


def add_sex_check_args(main_subparser):
    """Parser for the sex_check utility."""
    parser = main_subparser.add_parser(
        sex_check.SCRIPT_NAME, description=sex_check.DESCRIPTION,
    )

    parser.add_argument(
        "-v", "--version", action="version",
        version=f"pyGenClean {sex_check.SCRIPT_NAME} {__version__}",
    )

    subparsers = parser.add_subparsers(
        dest="subcommand", required=True,
        description=f"Below is a list of tools in the {sex_check.SCRIPT_NAME} "
                    f"module. Note that 'run' executes the main "
                    f"{sex_check.SCRIPT_NAME} pipeline.",
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
        plate_bias.SCRIPT_NAME, description=plate_bias.DESCRIPTION,
    )

    parser.add_argument(
        "-v", "--version", action="version",
        version=f"pyGenClean {plate_bias.SCRIPT_NAME} {__version__}",
    )

    subparsers = parser.add_subparsers(
        dest="subcommand", required=True,
        description=f"Below is a list of tools in the "
                    f"{plate_bias.SCRIPT_NAME} module. Note that 'run' "
                    f"executes the main {plate_bias.SCRIPT_NAME} pipeline.",
    )

    # Main plate-bias
    subparser = subparsers.add_parser(
        "run", description=plate_bias.DESCRIPTION, help=plate_bias.DESCRIPTION,
    )
    plate_bias.add_args(subparser)
    subparser.set_defaults(func=plate_bias.main)


def add_related_samples_args(main_subparser):
    """Parser for the related_samples utility."""
    parser = main_subparser.add_parser(
        related_samples.SCRIPT_NAME, description=related_samples.DESCRIPTION,
    )

    parser.add_argument(
        "-v", "--version", action="version",
        version=f"pyGenClean {related_samples.SCRIPT_NAME} {__version__}",
    )

    subparsers = parser.add_subparsers(
        dest="subcommand", required=True,
        description=f"Below is a list of tools in the "
                    f"{related_samples.SCRIPT_NAME} module. Note that 'run' "
                    f"executes the main {related_samples.SCRIPT_NAME} "
                    f"pipeline.",
    )

    # Main related-samples
    subparser = subparsers.add_parser(
        "run", description=related_samples.DESCRIPTION,
        help=related_samples.DESCRIPTION,
    )
    related_samples.add_args(subparser)
    subparser.set_defaults(func=related_samples.main)


def add_flag_maf(main_subparser):
    """Parser for the flag_maf utility."""
    parser = main_subparser.add_parser(
        flag_maf.SCRIPT_NAME, description=flag_maf.DESCRIPTION,
    )

    parser.add_argument(
        "-v", "--version", action="version",
        version=f"pyGenClean {flag_maf.SCRIPT_NAME} {__version__}",
    )

    subparsers = parser.add_subparsers(
        dest="subcommand", required=True,
        description=f"Below is a list of tools in the "
                    f"{flag_maf.SCRIPT_NAME} module. Note that 'run' "
                    f"executes the main {flag_maf.SCRIPT_NAME} pipeline.",
    )

    # Main flag-maf
    subparser = subparsers.add_parser(
        "run", description=flag_maf.DESCRIPTION, help=flag_maf.DESCRIPTION,
    )
    flag_maf.add_args(subparser)
    subparser.set_defaults(func=flag_maf.main)


def add_sample_call_rate(main_subparser):
    """Parser for the sample-call-rate utility."""
    parser = main_subparser.add_parser(
        sample_call_rate.SCRIPT_NAME, description=sample_call_rate.DESCRIPTION,
    )

    parser.add_argument(
        "-v", "--version", action="version",
        version=f"pyGenClean {sample_call_rate.SCRIPT_NAME} {__version__}",
    )

    subparsers = parser.add_subparsers(
        dest="subcommand", required=True,
        description=f"Below is a list of tools in the "
                    f"{sample_call_rate.SCRIPT_NAME} module. Note that 'run' "
                    f"executes the main {sample_call_rate.SCRIPT_NAME} "
                    f"pipeline.",
    )

    # Main sample-call-rate
    subparser = subparsers.add_parser(
        "run",
        description=sample_call_rate.DESCRIPTION,
        help=sample_call_rate.DESCRIPTION,
    )
    sample_call_rate.add_args(subparser)
    subparser.set_defaults(func=sample_call_rate.main)


def add_marker_call_rate(main_subparser):
    """Parser for the marker-call-rate utility."""
    parser = main_subparser.add_parser(
        marker_call_rate.SCRIPT_NAME, description=marker_call_rate.DESCRIPTION,
    )

    parser.add_argument(
        "-v", "--version", action="version",
        version=f"pyGenClean {marker_call_rate.SCRIPT_NAME} {__version__}",
    )

    subparsers = parser.add_subparsers(
        dest="subcommand", required=True,
        description=f"Below is a list of tools in the "
                    f"{marker_call_rate.SCRIPT_NAME} module. Note that 'run' "
                    f"executes the main {marker_call_rate.SCRIPT_NAME} "
                    f"pipeline.",
    )

    # Main marker-call-rate
    subparser = subparsers.add_parser(
        "run",
        description=marker_call_rate.DESCRIPTION,
        help=marker_call_rate.DESCRIPTION,
    )
    marker_call_rate.add_args(subparser)
    subparser.set_defaults(func=marker_call_rate.main)


def add_flag_hw(main_subparser):
    """Parser for the marker-call-rate utility."""
    parser = main_subparser.add_parser(
        flag_hw.SCRIPT_NAME, description=flag_hw.DESCRIPTION,
    )

    parser.add_argument(
        "-v", "--version", action="version",
        version=f"pyGenClean {flag_hw.SCRIPT_NAME} {__version__}",
    )

    subparsers = parser.add_subparsers(
        dest="subcommand", required=True,
        description=f"Below is a list of tools in the "
                    f"{flag_hw.SCRIPT_NAME} module. Note that 'run' "
                    f"executes the main {flag_hw.SCRIPT_NAME} "
                    f"pipeline.",
    )

    # Main flag-hw
    subparser = subparsers.add_parser(
        "run",
        description=flag_hw.DESCRIPTION,
        help=flag_hw.DESCRIPTION,
    )
    flag_hw.add_args(subparser)
    subparser.set_defaults(func=flag_hw.main)


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
