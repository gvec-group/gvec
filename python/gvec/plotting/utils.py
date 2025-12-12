# Copyright (c) 2025 GVEC Contributors, Max Planck Institute for Plasma Physics
# License: MIT
import numpy as np
from numpy import ndarray  # type checking


def _extrapolate_axis(state, quantities, zeta):
    """
    Generate an evaluations object and replace the axis values with those recovered from quadratic extrapolation
    """
    ev_axis = state.evaluate(
        *quantities,
        rho=[1e-4, 1.1e-4, 2.2e-4, 3.3e-4, 4.4e-4],  # must be >=1e-4
        theta=[0.0],
        zeta=zeta,
    )
    ev_axis = ev_axis.sel(theta=ev_axis.theta[0])
    # quadratic extrapolation to rho=0: y(0)= 3*(y(1) - y(2)) + y(3)
    ev_axis_1d_extrapol = ev_axis.sel(rho=ev_axis.rho[0])
    for var in ev_axis.data_vars:
        r1 = ev_axis.sel(rho=ev_axis.rho[1])[var].data
        r2 = ev_axis.sel(rho=ev_axis.rho[2])[var].data
        r3 = ev_axis.sel(rho=ev_axis.rho[3])[var].data
        # r4=ev_axis.sel(rho=ev_axis.rho[4])[var].data
        # print(f" check extrapolation, variable:{var} max diff:{np.amax(np.abs(r1 -  (3*(r2 - r3) + r4)))}")
        # OVERWRITE DATA WITH EXTRAPOLATION works with a 1d array...
        ev_axis_1d_extrapol[var].data = 3 * (r1 - r2) + r3

    return ev_axis_1d_extrapol


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


def _subplots(subplot_grid, **kwargs):
    """
    we will ensure this is always returns a `numpy.ndarray` of axis objects
    """
    import matplotlib.pyplot as plt

    nrow = subplot_grid[0]
    ncol = subplot_grid[1]

    f, ax = plt.subplots(nrow, ncol, layout="constrained", **kwargs)
    ax = np.asarray(ax).flatten()
    return f, ax
