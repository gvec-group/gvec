import matplotlib.pyplot as plt
from numpy import array, asarray, linspace, max, ndarray, pi, prod


def _get_coord_range(coordinate, nfp, points, min_val=0.0):
    """
    Boiler plate code for building the arrays required for plotting

    Return is always of type `numpy.ndarray`.
    """

    # Limits chosen based on the coordinate
    if coordinate == "rho":
        auto_limit = [min_val, 1]
    elif coordinate == "theta":
        auto_limit = [min_val, 2 * pi]
    elif coordinate == "zeta":
        auto_limit = [min_val, 2 * pi / nfp]

    # The output will be adjusted depending on the input type
    if isinstance(points, int):
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


def _subplots(nrow, ncol, **kwargs):
    """
    we will ensure this is always returns a `numpy.ndarray` of axis objects
    """
    f, ax = plt.subplots(nrow, ncol, **kwargs)
    ax = asarray(ax).flatten()
    return f, ax
