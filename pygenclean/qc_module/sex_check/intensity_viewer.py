"""Small Dash application to view intensities."""


import logging
import argparse
from os import path
from typing import Optional, List

import dash
import plotly.graph_objects as go

import pandas as pd

from .intensity_plot import PLOT_CONFIG

from ...error import ProgramError
from ...version import pygenclean_version as __version__


SCRIPT_NAME = "intensity-viewer"
DESCRIPTION = "Small Dash application to view intensities"


logger = logging.getLogger(__name__)


def main(args: Optional[argparse.Namespace] = None,
         argv: Optional[List[str]] = None) -> None:
    """Create an interactive intensity viewer."""
    if args is None:
        args = parse_args(argv)
    check_args(args)

    # Reading the data
    df = pd.read_csv(args.summarized_intensities, sep="\t")\
        .set_index("sample_id", verify_integrity=True)

    # Creating the dash application
    application = create_application()

    # Adding the layout to the application
    add_layout(application, df.index.to_list())

    # Adding the callback
    add_callbacks(application, df)

    # Serving the Dash application
    application.run_server(debug=args.debug, port=args.port)


def create_application() -> dash.Dash:
    """Creates a simple Dash application."""
    application = dash.Dash(__name__)
    return application


def add_layout(application: dash.Dash, samples: List[str]) -> None:
    """Adds the layout to a Dash application."""
    # Layout for the sample selector
    application.layout = dash.html.Div(
        className="page-content",
        children=[
            dash.html.H1(
                "pyGenClean - Intensity viewer", className="page-title",
            ),
            dash.html.H2("Highlight samples", className="page-section"),
            dash.html.Div(
                className="sample-selector",
                children=[
                    dash.dcc.Dropdown(
                        id="selected-samples",
                        options=[
                            {"label": sample, "value": sample}
                            for sample in samples
                        ],
                        value=[],
                        multi=True,
                        className="sample-selector",
                        placeholder="Select one or more samples to "
                                    "highlight...",
                    ),
                ],
            ),
            dash.dcc.Graph(id="sex-check-graph"),
        ],
    )


def add_callbacks(app: dash.Dash, df: pd.DataFrame) -> None:
    """Adds callbacks to the Dash application."""
    @app.callback(
        dash.Output("sex-check-graph", "figure"),
        dash.Input("selected-samples", "value"),
    )
    def update_graph(selected_samples: List[str]) -> go.Figure:
        return create_sexcheck_figure(df, selected_samples)


def create_sexcheck_figure(
    df: pd.DataFrame, selected_samples: List[str]
) -> go.Figure:
    """Creates a Plotly figure containing sex-check intensities."""
    figure = go.Figure()

    # Adding all the points (excluding selected samples)
    remove = set(selected_samples)
    for status in ("OK", "Mismatch"):
        for sex in ("Male", "Female", "Unknown"):
            subset = (
                (df["sex"] == sex)
                & (df["status"] == status)
                & (~df.index.isin(remove))
            )
            sub_df = df.loc[subset, :]

            # Getting the points' config
            config = PLOT_CONFIG[(status, sex)]

            figure.add_trace(
                go.Scatter(
                    mode="markers",
                    x=sub_df["chr23"],
                    y=sub_df["chr24"],
                    name=f"{status} {sex} (n={sub_df.shape[0]:,d})",
                    text=sub_df.index.to_numpy(),
                    marker={
                        "color": config["color"],
                        "size": config["plotly_size"],
                        "symbol": config["plotly_symbol"],
                        "line": config.get("plotly_line", None),
                    },
                ),
            )

    for sample_id in selected_samples:
        sample_data = df.loc[[sample_id], :]
        assert sample_data.shape[0] == 1

        status = sample_data.iloc[0, :]["status"]
        sex = sample_data.iloc[0, :]["sex"]

        config = PLOT_CONFIG[(status, sex)]

        figure.add_trace(
            go.Scatter(
                mode="markers",
                x=sample_data["chr23"],
                y=sample_data["chr24"],
                name=sample_id,
                text=sample_data["status"] + " " + sample_data["sex"],
                marker={
                    "color": config["color"],
                    "size": config["highlighted_plotly_size"],
                    "symbol": config["plotly_symbol"],
                    "line": config["highlighted_plotly_line"],
                },
            ),
        )

    # Adding the title
    nb_samples = df.shape[0]
    figure.update_layout(
        title=f"<b>Summarized Intensities</b><br />"
              f"<sup><b>{nb_samples:,d} samples</b></sup>",
        title_x=0.5,
    )

    # Adding the x and y labels
    figure.update_xaxes(
        title_text="chrX",
    )
    figure.update_yaxes(
        title_text="chrY",
    )

    return figure


def check_args(args: argparse.Namespace) -> None:
    """Checks the arguments and options.

    Args:
        args (argparse.Namespace): the arguments and options.

    """
    if not path.isfile(args.summarized_intensities):
        raise ProgramError(f"{args.summarized_intensities}: no such file")


def parse_args(argv: Optional[List[str]] = None) -> argparse.Namespace:
    """Parses the arguments and options."""
    parser = argparse.ArgumentParser(description=DESCRIPTION)

    parser.add_argument(
        "-v", "--version", action="version",
        version=f"pyGenClean {SCRIPT_NAME} {__version__}",
    )

    # Adding the arguments
    add_args(parser)

    return parser.parse_args(argv)


def add_args(parser: argparse.ArgumentParser) -> None:
    """Adds arguments to the parser."""
    group = parser.add_argument_group("Input files")
    group.add_argument(
        "--summarized-intensities", type=str, metavar="FILE",
        help="The name of the file containing the summarized intensities.",
    )

    group = parser.add_argument_group("Options")
    group.add_argument(
        "--debug", action="store_true",
        help="Launch application in DEBUG mode.",
    )
    group.add_argument(
        "-p", "--port", type=int, default=8050,
        help="The port used by the Dash application.",
    )
