"""Illumina utility tools."""


from ..utils import get_open_func


__all__ = ["get_data_position"]


def get_data_position(filename):
    """Finds the number of lines to skip before the [Data] tag.

    Args:
        filename (str): the filename to check.

    Returns:
        int: the number of line(s) to skip.

    """
    # Checking the first hundred lines to check if [Data] tag is present
    start = 0
    with get_open_func(filename)() as f:
        for i, line in enumerate(f):
            if i == 100:
                break

            if line.startswith("[Data]"):
                start = i + 1
                break

    return start
