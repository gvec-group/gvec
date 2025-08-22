# Copyright (c) 2025 GVEC Contributors, Max Planck Institute for Plasma Physics
# License: MIT
from pathlib import Path
import logging

from pyevtk.hl import gridToVTK
import xarray as xr
import numpy as np


def ev2vtk(
    filename: Path | str,
    xrds: xr.Dataset,
    quiet: bool = True,
):
    """
    Write a GVEC evaluation dataset to a VTS file.

    Parameters
    ----------
    filename : str
        The name of the output file without the '.vts' extension.
    xrds : xr.Dataset
        The dataset containing the evaluation data.
    quiet : bool, optional
        If False, return the path to the output file, by default True.

    Returns
    -------
    Path
        The path to the output file.

    Notes
    -----
    The following dimension are expected to be in the dataset:
    - 'pos' : the cartesian components of grid points
    - 'xyz' : the dimension name for the cartesian components of grid points
    - 'rad' : the radial dimension name
    - 'pol' : the poloidal dimension name
    - 'tor' : the toroidal dimension name

    Scalar variables without the 'xyz' dimension are broadcasted to the 'rad', 'pol', 'tor' dimensions.

    If a variable does not have the expected dimensions, it is ignored.

    Examples
    --------
    >>> from gvec.vtk import ev2vtk
    >>> import xarray as xr
    >>> filename = "my_evaluation"
    >>> xrds = xr.Dataset({"pos": (["xyz", "rad", "pol", "tor"], np.random.rand(3, 10, 10, 10))})
    >>> ev2vtk(filename, xrds)
    """
    # pyevtk expects a string
    if isinstance(filename, Path):
        filename = str(filename)

    # name of the cartesian components of grid points
    position_vector = "xyz"
    cart_pos_vector = "pos"

    # make sure dimensions are in the expected order
    dimension_order = ["xyz", "rad", "pol", "tor"]

    assert (
        "pos" in xrds
    ), """Expected 'pos' in 'xrds', please make sure you are working with a pygvec evaluation dataset
    or rename your variable for the  cartesian components of grid points to 'pos'."""

    expected_dimension = {"rad": "radial", "pol": "poloidal", "tor": "toroidal"}
    for dim in expected_dimension:
        assert (
            dim in xrds.dims
        ), f"""Expected '{dim}' in 'xrds' dimensions, please make sure you are working with a pygvec evaluation dataset
        or rename your {expected_dimension[dim]} dimension to '{dim}'."""

    outvars = []
    ignored_variables = []
    for var in xrds.data_vars:
        if set(xrds[var].dims).issubset(dimension_order) and len(xrds[var].dims) >= 1:
            outvars.append(var)
        else:
            ignored_variables.append(var)

    # variables without the "xyz" dimension
    scalar_vars = [var for var in outvars if (position_vector not in xrds[var].dims)]

    broadcast_like_scalar_var = xr.DataArray(
        np.zeros((xrds.sizes["rad"], xrds.sizes["pol"], xrds.sizes["tor"])),
        dims=("rad", "pol", "tor"),
    )

    # variables with the "xyz" dimension
    vector_vars = [var for var in outvars if (position_vector in xrds[var].dims)]

    # vector of the cartesian components of grid points
    xcoord, ycoord, zcoord = xrds[cart_pos_vector].transpose(*dimension_order).values

    # point data handed to gridToVTK
    ptdata = {}

    # broadcasting of the coordinates to rad, pol, tor
    for coord in xrds.coords:
        if position_vector == coord:
            continue

        coord_reshaped = xrds[coord].broadcast_like(broadcast_like_scalar_var)
        coord_reshaped = coord_reshaped.transpose(*dimension_order[1:])
        ptdata[coord] = np.ascontiguousarray(coord_reshaped.values)

    # broadcasting and storing of the scalar variables to rad, pol, tor
    for var in scalar_vars:
        if var == cart_pos_vector:
            continue
        if len(xrds[var].dims) < 3:
            var_values = xrds[var]
            var_values = var_values.broadcast_like(broadcast_like_scalar_var)
        else:
            var_values = xrds[var]
        var_values = var_values.transpose(*dimension_order[1:]).values
        ptdata[var] = np.ascontiguousarray(var_values)

    # storing of the vector variables
    for var in vector_vars:
        if var == cart_pos_vector:
            continue
        vx, vy, vz = xrds[var].transpose(*dimension_order).values
        ptdata[var] = (
            np.ascontiguousarray(vx),
            np.ascontiguousarray(vy),
            np.ascontiguousarray(vz),
        )

    # NOTE: gridToVTK expects C_contiguous or F_contiguous arrays and does not support Path for filenames
    fn = gridToVTK(
        filename,
        np.ascontiguousarray(xcoord),
        np.ascontiguousarray(ycoord),
        np.ascontiguousarray(zcoord),
        pointData=ptdata,
    )

    if len(ignored_variables) != 0:
        logging.warning(
            f"The following varivables are ignored and not written to {filename}.vts: {ignored_variables}."
        )

    if not quiet:
        return Path(fn)
