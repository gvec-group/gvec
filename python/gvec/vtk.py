# Copyright (c) 2025 GVEC Contributors, Max Planck Institute for Plasma Physics
# License: MIT
from pathlib import Path
import logging

from pyevtk.hl import gridToVTK
import xarray as xr
import numpy as np
from gvec import fourier


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


def gframe_to_vtk(
    file: str | Path,
    prefix="visu",
    zeta_visu: np.ndarray = None,
    theta_visu: np.ndarray = None,
    phi_visu: np.ndarray = None,
    box_axis=None,
    filetype="vts",
):
    """
    Reads a netcdf file that defines the G-Frame and possibly the boundary X1,X2 in that frame.
    The file is for example produced by the GVEC quasr script, and then used in GVEC for initialization.
    Writes vtk-visualization files from the data.
    Input:
        * file : netcdf file that contains axis and boundary data. Format as produced by the GVEC quasr script.
        * prefix : prefix of the output files. Default is `visu_`.
        * zeta_visu : 1d zeta positions of the axis and boundary surface output. If not specified, the ones from the file are used.
        * theta_visu : 1d theta positions of the boundary surface output. If not specified, the ones from the file are used.
        * box_axis : if =[a,b], visualize G-Frame additionally as a box of with distances +a -a in N direction and +b -b in B direction.
        * filetype : can be "vts"  (VTK) or "nc" (netcdf)
    Output:
        * writes `prefix_axis.filetype` : if 'axis' group exists in `file`, provides the origin curve position in 3D and N,B vectors on that curve. On full torus or on given `zeta_visu` positions
        * writes `prefix_boundary.filetype` : if 'boundary' group exists `file`, provides the boundary surface position in 3D. On one field period, or on given `zeta_visu` positions
        * writes `prefix_axis_box.filetype` : if box_axis=[a,b], G-Frame is visualized as a box aroud the axis.
    """

    ds_main = xr.open_dataset(file, engine="netcdf4")
    nfp = ds_main.NFP.data
    try:
        ds_axis = xr.open_dataset(file, engine="netcdf4", group="axis")
    except Exception as e:
        raise RuntimeError(f"Could not open axis group in {file}") from e

    zeta_fp = ds_axis["zeta(:)"].data
    pos_axis = ds_axis["xyz(::)"].data
    N_axis = ds_axis["Nxyz(::)"].data
    B_axis = ds_axis["Bxyz(::)"].data
    ds_axis.close()

    nzeta_fp_axis = len(zeta_fp)
    nzeta_full = N_axis.shape[1]
    assert nzeta_full == nzeta_fp_axis * nfp, (
        f"axis data must be given on a full turn, with nfp being a factor in the number of points! nfp={nfp}, nzeta_full={nzeta_full}, nzeta of one fp={nzeta_fp_axis}"
    )
    zeta_axis = zeta_fp[0] + np.linspace(0, 2 * np.pi, nzeta_full, endpoint=False)

    assert np.allclose(zeta_fp, zeta_axis[0:nzeta_fp_axis]), "zeta on axis must be equidistant"
    if zeta_visu is not None:
        zdft = fourier.real_dft_mat(zeta_axis, zeta_visu, nfp=1)
        zeta_out = zeta_visu
        pos = pos_axis @ zdft["BF"].T
        N = N_axis @ zdft["BF"].T
        B = B_axis @ zdft["BF"].T
    else:
        zeta_out = zeta_axis
        pos = pos_axis
        N = N_axis
        B = B_axis

    # convert to xarray
    xr_axis = xr.Dataset(
        coords=dict(
            rho=("rad", [0.0]),
            theta=("pol", [0]),
            zeta=("tor", zeta_out),
            xyz=("xyz", [0, 1, 2]),
        ),
        data_vars=dict(
            pos=(["xyz", "rad", "pol", "tor"], pos.reshape((3, 1, 1, len(zeta_out)))),
            N=(["xyz", "rad", "pol", "tor"], N.reshape((3, 1, 1, len(zeta_out)))),
            B=(["xyz", "rad", "pol", "tor"], B.reshape((3, 1, 1, len(zeta_out)))),
        ),
    )
    if filetype == "vts":
        ev2vtk(f"{prefix}_axis", xr_axis)
    elif filetype == "nc":
        xr_axis.to_netcdf(f"{prefix}_axis.nc", mode="w")
    else:
        raise ValueError(f"unknown filetype {filetype}, only 'vts' and 'nc' supported.")

    # optional box visualization
    if box_axis is not None:
        assert len(box_axis) == 2, "box_axis input must be a list of two values"
        X = np.array([-1, 1]) * box_axis[0]
        Y = np.array([-1, 1]) * box_axis[1]

        pos_box = (
            pos[:, None, None, :]
            + N[:, None, None, :] * X[None, :, None, None]
            + B[:, None, None, :] * Y[None, None, :, None]
        )
        xr_box = xr.Dataset(
            coords=dict(
                rho=("rad", [-1, 1]),  # =X
                theta=("pol", [-1, 1]),  # =Y
                zeta=("tor", zeta_out),
                xyz=("xyz", [0, 1, 2]),
            ),
            data_vars=dict(
                pos=(["xyz", "rad", "pol", "tor"], pos_box),
            ),
        )
        if filetype == "vts":
            ev2vtk(f"{prefix}_box_axis", xr_box)
        elif filetype == "nc":
            xr_box.to_netcdf(f"{prefix}_box_axis.nc", mode="w")

    # read boundary group
    try:
        ds_boundary = xr.open_dataset(file, engine="netcdf4", group="boundary")
    except Exception as e:
        print(
            f" {e}\n Warning: boundary group not found in {file}. boundary visualization skipped."
        )
        return

    theta_bnd = ds_boundary["theta(:)"].data
    zeta_bnd = ds_boundary["zeta(:)"].data
    X = ds_boundary["X(::)"].data
    Y = ds_boundary["Y(::)"].data
    ds_boundary.close()

    # if necessary, apply fourier.real_dft_mat to get axis and boundary positions
    if zeta_visu is None:
        zeta_out = zeta_bnd
        if len(zeta_bnd) == nzeta_fp_axis:
            # same zeta points in axis and boundary:
            pos = pos_axis[:, None, 0:nzeta_fp_axis]
            N = N_axis[:, None, 0:nzeta_fp_axis]
            B = B_axis[:, None, 0:nzeta_fp_axis]
            XX = X
            YY = Y
        else:
            # sample axis on zeta_bnd points:
            zdft = fourier.real_dft_mat(zeta_axis, zeta_bnd)
            pos = pos_axis @ zdft["BF"].T
            N = N_axis @ zdft["BF"].T
            B = N_axis @ zdft["BF"].T
            XX = X
            YY = Y
    else:
        zeta_out = zeta_visu
        # sample axis
        zdft = fourier.real_dft_mat(zeta_axis, zeta_visu, nfp=1)
        pos = pos_axis @ zdft["BF"].T
        N = N_axis @ zdft["BF"].T
        B = B_axis @ zdft["BF"].T
        # sample boundary
        zdft = fourier.real_dft_mat(zeta_bnd, zeta_visu, nfp=nfp)
        XX = X @ zdft["BF"].T
        YY = Y @ zdft["BF"].T

    if theta_visu is None:
        theta_out = theta_bnd
    else:
        theta_out = theta_visu
        tdft = fourier.real_dft_mat(theta_bnd, theta_visu)
        XX = tdft["BF"] @ XX
        YY = tdft["BF"] @ YY
    pos_bnd = pos[:, None, :] + XX[None, :, :] * N[:, None, :] + YY[None, :, :] * B[:, None, :]

    # convert to xarray
    xr_bnd = xr.Dataset(
        coords=dict(
            rho=("rad", [1.0]),
            theta=("pol", theta_out),
            zeta=("tor", zeta_out),
            xyz=("xyz", [0, 1, 2]),
        ),
        data_vars=dict(
            pos=(
                ["xyz", "rad", "pol", "tor"],
                pos_bnd.reshape((3, 1, len(theta_out), len(zeta_out))),
            ),
            X1=(["pol", "tor"], XX),
            X2=(["pol", "tor"], YY),
        ),
    )
    # write visualization file
    if filetype == "vts":
        ev2vtk(f"{prefix}_boundary", xr_bnd)
    elif filetype == "nc":
        xr_bnd.to_netcdf(f"{prefix}_boundary.nc", mode="w")
