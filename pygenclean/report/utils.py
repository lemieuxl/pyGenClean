"""Utility fonctions for reporting."""


import re


__all__ = ["format_numbers"]


def format_numbers(number: float) -> str:
    """Formats number in the scientific notation for LaTeX.

    Args:
        number (str): the number to format.

    Returns:
        str: a string containing the scientific notation of the number.

    """
    # Changing to str
    number = str(number)

    # Matching
    match = re.match(r"^([-+]?\d*\.\d+|\d+)e([-+]?\d+)$", number)

    # Nothing matched
    if not match:
        return number

    # Getting the coefficient and the exponent
    coefficient = match.group(1)
    exponent = int(match.group(2))

    return coefficient + r"\times 10^{" + str(exponent) + "}"
