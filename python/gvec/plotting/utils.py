import matplotlib.pyplot as plt
from numpy import array, asarray, floor, linspace, max, ndarray, pi, prod, sqrt


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


def _get_coord_range(coordinate, nfp, points, bounds=None):
    """
    Boiler plate code for building the arrays required for plotting

    Return is always of type `numpy.ndarray`.
    """

    # The output will be adjusted depending on the input type
    if isinstance(points, int):
        # Limits chosen based on the coordinate
        if bounds is None:
            if coordinate == "rho":
                auto_limit = [1e-4, 1.0]
            elif coordinate == "theta":
                auto_limit = [0.0, 2 * pi]
            elif coordinate == "zeta":
                auto_limit = [0.0, 2 * pi / nfp]
        else:
            auto_limit = bounds
        nodes = linspace(*auto_limit, points)
    elif isinstance(points, ndarray):
        nodes = points
    elif isinstance(points, float):
        nodes = array([points])

    return nodes


def _get_scalars_for_plotting(evaluations, equilibrium_quantities, post_process, direction):
    """
    Move all eval quantities into a dict for plotting. Post process the quantities if required.

    Return is a `dict{"quantity": array([values])}

    TODO: Allow post processing even if the quantity is a scalar.
    TODO: Switch from using dicts to xarray object
    """

    plotting_quantities = {}
    for quantity in equilibrium_quantities:
        # Loop over the quantities and store what we're plotting in the dict
        if max(evaluations[quantity].shape) != prod(evaluations[quantity].shape):
            # If this quantity is not a scalar
            #   we check to see if it is going to be remapped, or we error out
            if post_process is None:
                raise TypeError("The plotted quantities must be scalars.")
            else:
                # Post process the quantity and give the new quantity the given name in plotting
                tmp = post_process[quantity][0](evaluations[quantity])
                if len(tmp.shape) == 1:
                    plotting_quantities[post_process[quantity][1]] = tmp
                else:
                    raise ValueError("Post-processing function does not return a 1D array.")
        else:
            plotting_quantities[quantity] = evaluations[quantity].data.flatten()

    return plotting_quantities


def _add_postprocess_to_xarray(evaluations, quantities, post_process: dict):
    """
    Post process the quantities if required.
    """
    for value, process in post_process.items():
        process_check = all([val in process.keys() for val in ["name", "function"]])
        if not process_check:
            raise KeyError(
                f"post_process entries must be a dict containing at least keys 'name' and 'function'. Check entries for {value}."
            )

        # Add the new value to the evaluations object
        evaluations[process["name"]] = process["function"](evaluations[value])

        # If `symbol` is present in the dict we'll use this for replacing the symbol
        # otherwise use the name
        if "symbol" in process.keys():
            evaluations[process["name"]].attrs["symbol"] = process["symbol"]
        else:
            evaluations[process["name"]].attrs["symbol"] = process["name"]

        if len(evaluations[process["name"]]) != prod(evaluations[process["name"]].size):
            # Make sure that we are plotting a scalar
            raise ValueError("Post-processing function does not return a 1D array.")
        for i in range(len(quantities)):
            # swap out the plotting values in the quantities list so that we plot the correct thing
            if quantities[i] == value:
                quantities[i] = process["name"]

    return evaluations, quantities


def _subplots(nrow, ncol, **kwargs):
    """
    we will ensure this is always returns a `numpy.ndarray` of axis objects
    """
    f, ax = plt.subplots(nrow, ncol, **kwargs)
    ax = asarray(ax).flatten()
    return f, ax


def _design_subgrid(nplots):
    """
    Find some numbers which divide the total number of subplots to make a subplot grid with.
    """
    divisor_a = 0
    divisor_b = 0

    divisor = sqrt(nplots)
    if floor(divisor) ** 2 != nplots:
        # If not a perfect square we will try and find the two numbers which give us this one
        i = 1
        while i < floor(nplots / 2):
            i += 1
            check = nplots / i
            if check == floor(check):
                divisor_a = check
                divisor_b = i
    elif floor(divisor) ** 2 == nplots:
        # If a perfect square we can finish
        divisor_a = divisor_b = divisor

    if divisor_a == 0:
        # If we couldn't find a divisor we will return one long row.
        divisor_a = 1
        divisor_b = nplots

    return [int(divisor_a), int(divisor_b)]
