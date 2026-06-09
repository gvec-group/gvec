# Copyright (c) 2025 GVEC Contributors, Max Planck Institute for Plasma Physics
# License: MIT
import functools

import matplotlib.pyplot as plt
import numpy as np


def _symbol_check(evaluations, quantities):
    """
    Check the 'symbol' attribute in the evaluations.
    """
    for quantity in quantities:
        # If symbol is not present use the quantity
        if "symbol" in evaluations[quantity].attrs:
            symbol = evaluations[quantity].symbol
        else:
            # Since it will be tex-ified make sure it will show as a string correctly
            symbol = "\\mathrm{" + quantity + "}"

        if "\\text" in symbol:
            # there may be issues with \\text if amsmath is not loaded
            symbol = symbol.replace("\\text", "\\mathrm")
        evaluations[quantity].attrs["symbol"] = symbol

    return evaluations


def _design_subgrid(nplots):
    """
    Find some numbers which divide the total number of subplots to make a subplot grid with.
    """
    divisor_a = 0
    divisor_b = 0

    divisor = np.sqrt(nplots)
    if np.floor(divisor) ** 2 != nplots:
        # If not a perfect square we will try and find the two numbers which give us this one
        i = 1
        while i < np.floor(nplots / 2):
            i += 1
            check = nplots / i
            if check == np.floor(check):
                divisor_a = check
                divisor_b = i
    elif np.floor(divisor) ** 2 == nplots:
        # If a perfect square we can finish
        divisor_a = divisor_b = divisor

    if divisor_a == 0:
        # If we couldn't find a divisor we will return one long row.
        divisor_a = 1
        divisor_b = nplots

    return [int(divisor_a), int(divisor_b)]


def _subplots(subplot_grid, sharex, sharey, **kwargs):
    """
    we will ensure this is always returns a `numpy.ndarray` of axis objects
    """

    nrow = subplot_grid[0]
    ncol = subplot_grid[1]

    f, ax = plt.subplots(
        nrow, ncol, layout="compressed", sharex=sharex, sharey=sharey, **kwargs
    )
    return f, ax


def _convience_axis_additions(state, ax, quantity):
    ax.gvec_quantity = quantity
    ax.nfp = state.nfp
