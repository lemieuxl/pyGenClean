"""Scatter plot utilities."""


import jinja2


__all__ = ["generate_html_scatter"]


def generate_html_scatter(filename, **kwargs):
    """Generates an HTML scatter plot using Plotly (with a template).

    Args:
        out_filename (str): the name of the HTML file.

    The extra keyword will be passed to the template. The following information
    is required:

    - ``page_title``: the title of the HTML page.
    - ``scatter_title``: the title of the plot.
    - ``xlabel``: the name of the x axis.
    - ``ylabel``: the name of the y axis.
    - ``plot_width``: the width of the plot (optional).
    - ``plot_height``: the height of the plot (optional).
    - ``labels``: a list of labels (one for each data group).
    - ``data``: a dictionary where keys are the labels, and values are ``x``,
      ``y``, and ``ids``.
    - ``colors``: a dictionary containing the color of each label.
    - ``sizes``: a dictionary containing the size of each label.
    - ``symbols``: a dictionary with marker's symbol for each label.
    - ``lines``: a dictionary describing lines (optional).

    """
    # Getting the template
    jinja2_env = jinja2.Environment(
        loader=jinja2.PackageLoader(__name__, "templates"),
    )
    template = jinja2_env.get_template("plotly_scatter.html")

    # Rendering the template
    with open(filename, "w") as f:
        print(template.render(**kwargs), file=f)
