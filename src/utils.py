import numpy as np
from typing import Dict, List, Tuple


def format_sig_figs(value: float, sig_figs=3) -> str:
    """
    this function formats a floating-point number to a specified number of significant figures

    args:
        value (float): the floating-point number to be formatted
        sig_figs (int, optional): the number of significant figures to display. defaults to 3

    returns:
        str: the formatted string representation of the number with the specified significant figures
    """
    if not np.isfinite(value): # check if the value is nan or inf
        return f"{value:.{sig_figs}g}" # use general format for non-finite values
    if value == 0:
        return f"{value:.{sig_figs}f}"
    else:
        digits = sig_figs - int(np.floor(np.log10(abs(value)))) - 1
        if digits < 0: # ensure we don't pass a negative value to the format string
            digits = 0
        return f"{value:.{digits}f}"
