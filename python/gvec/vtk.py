# Copyright (c) 2025 GVEC Contributors, Max Planck Institute for Plasma Physics
# License: MIT
from pathlib import Path
import logging

from pyevtk.hl import gridToVTK
import xarray as xr
import numpy as np
from gvec import fourier
from gvec import gframe


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
    The following data / dimensions are expected to be in the dataset:
    - 'xyz' : the dimension name for the cartesian components of grid points
    - 'rad' /'rho'   : the radial dimension name. If not present in xrds, it is added as {"rad":[0]}
    - 'pol' /'theta' : the poloidal dimension name. If not present in xrds, it is added as {"pol":[0]}
    - 'tor' /'zeta'  : the toroidal dimension name. If not present in xrds, it is added as {"tor":[0]}
    - 'pos' : datarray with the cartesian components of grid points, with dimension 'xyz' and at least one of ['rad', 'pol', 'tor']

    Scalar variables are dataarrays without the 'xyz' dimension, and  are broadcasted to the 'rad', 'pol', 'tor' dimensions.
    Vector variables are dataarrays with the 'xyz' dimension, and are broadcasted to the 'xyz', 'rad', 'pol', 'tor' dimensions.

    If a dataarray (except 'pos') does not have the expected dimensions, it is ignored.

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
    ds_out = xrds.copy()
    assert (
        "pos" in ds_out
    ), """Expected 'pos' in 'xrds', please make sure you are working with a pygvec evaluation dataset
    or rename your variable for the  cartesian components of grid points to 'pos'."""

    # rename dimensions
    for dim, dimnew in (("rho", "rad"), ("theta", "pol"), ("zeta", "tor")):
        if dim in ds_out.pos.dims:
            ds_out = ds_out.rename_dims({dim: dimnew})

    assert np.any([dim in ds_out.pos.dims for dim in ("rad", "pol", "tor")]), (
        'pos data array must contain at least one of "rad", "pol", "tor" dimensions.'
    )
    # add missing dimensions to pos
    expanded_dim = []
    for dim in ("rad", "pol", "tor"):
        if dim not in ds_out.pos.dims:
            expanded_dim.append(dim)
            ds_out["pos"] = ds_out.pos.expand_dims({dim: [0.0]})

    expected_dimension = {"rad": "radial", "pol": "poloidal", "tor": "toroidal"}

    for dim in expected_dimension:
        assert (
            dim in ds_out.dims
        ), f"""Expected '{dim}' in 'xrds' dimensions, please make sure you are working with a pygvec evaluation dataset
        or rename your {expected_dimension[dim]} dimension to '{dim}'."""

    outvars = []
    ignored_variables = []
    for var in ds_out.data_vars:
        if set(ds_out[var].dims).issubset(dimension_order) and len(ds_out[var].dims) >= 1:
            outvars.append(var)
        else:
            ignored_variables.append(var)

    # variables without the "xyz" dimension
    scalar_vars = [var for var in outvars if (position_vector not in ds_out[var].dims)]

    broadcast_like_scalar_var = xr.DataArray(
        np.zeros((ds_out.sizes["rad"], ds_out.sizes["pol"], ds_out.sizes["tor"])),
        dims=("rad", "pol", "tor"),
    )
    broadcast_like_vector_var = xr.DataArray(
        np.zeros((ds_out.sizes["rad"], ds_out.sizes["pol"], ds_out.sizes["tor"], 3)),
        dims=("rad", "pol", "tor", "xyz"),
    )

    # variables with the "xyz" dimension
    vector_vars = [var for var in outvars if (position_vector in ds_out[var].dims)]

    # vector of the cartesian components of grid points
    xcoord, ycoord, zcoord = ds_out[cart_pos_vector].transpose(*dimension_order).values

    # point data handed to gridToVTK
    ptdata = {}

    # broadcasting of the coordinates to rad, pol, tor
    for coord in ds_out.coords:
        if position_vector == coord:
            continue
        if coord in expanded_dim:
            continue
        coord_reshaped = ds_out[coord].broadcast_like(broadcast_like_scalar_var)
        coord_reshaped = coord_reshaped.transpose(*dimension_order[1:])
        ptdata[coord] = np.ascontiguousarray(coord_reshaped.values)

    # broadcasting and storing of the scalar variables to rad, pol, tor
    for var in scalar_vars:
        if var == cart_pos_vector:
            continue
        var_values = ds_out[var]
        if len(ds_out[var].dims) < 3:
            var_values = var_values.broadcast_like(broadcast_like_scalar_var)
        var_values = var_values.transpose(*dimension_order[1:]).values
        ptdata[var] = np.ascontiguousarray(var_values)

    # storing of the vector variables
    for var in vector_vars:
        if var == cart_pos_vector:
            continue
        var_values = ds_out[var]
        if len(var_values.dims) < 4:
            var_values = var_values.broadcast_like(broadcast_like_vector_var)
        vx, vy, vz = var_values.transpose(*dimension_order).values
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
    visu_boundary=True,
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
        * box_axis : if =[a>0,b>0], visualize G-Frame additionally as a box of with distances +a -a in N direction and +b -b in B direction.
        * visu_boundary : if True, visualize the boundary surface
        * filetype : can be "vts"  (VTK) or "nc" (netcdf)
    Output:
        * writes `prefix_axis.filetype` : if 'axis' group exists in `file`, provides the origin curve position in 3D and N,B vectors on that curve. On full torus or on given `zeta_visu` positions
        * writes `prefix_boundary.filetype` : if 'boundary' group exists `file`, provides the boundary surface position in 3D. On one field period, or on given `zeta_visu` positions
        * writes `prefix_axis_box.filetype` : if box_axis=[a,b], G-Frame is visualized as a box aroud the axis.
    """
    dnc = gframe.read_Gframe_ncfile(file)
    nfp = dnc["nfp"]
    if "axis" in dnc:
        d_axis = dnc["axis"]
    else:
        raise RuntimeError(f"Could not open axis group in {file}")

    zeta_fp = d_axis["zeta"]
    zeta_axis = d_axis["zetafull"]
    pos_axis = d_axis["xyz"]
    N_axis = d_axis["Nxyz"]
    B_axis = d_axis["Bxyz"]
    nzeta_fp_axis = len(zeta_fp)

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
            zeta=("tor", zeta_out),
            xyz=("xyz", [0, 1, 2]),
        ),
        data_vars=dict(
            pos=(["xyz", "tor"], pos),
            N=(["xyz", "tor"], N),
            B=(["xyz", "tor"], B),
        ),
    )
    if filetype == "vts":
        ev2vtk(f"{prefix}_axis", xr_axis)
    elif filetype == "nc":
        xr_axis.to_netcdf(f"{prefix}_axis.nc", mode="w", engine="netcdf4")
    else:
        raise ValueError(f"unknown filetype {filetype}, only 'vts' and 'nc' supported.")

    # optional box visualization
    if box_axis is not None:
        assert len(box_axis) == 2, "box_axis input must be a list of two values"
        X = np.array([-1, 1]) * box_axis[0]
        Y = np.array([-1, 1]) * box_axis[1]

        pos_out = (
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
                pos=(["xyz", "rad", "pol", "tor"], pos_out),
            ),
        )
        if filetype == "vts":
            ev2vtk(f"{prefix}_box_axis", xr_box)
        elif filetype == "nc":
            xr_box.to_netcdf(f"{prefix}_box_axis.nc", mode="w", engine="netcdf4")

    # read boundary group
    if not visu_boundary:
        return
    if "boundary" in dnc:
        d_bnd = dnc["boundary"]
    else:
        raise RuntimeError(f"Could not open boundary group in {file}")

    theta_bnd = d_bnd["theta"]
    zeta_bnd = d_bnd["zeta"]
    X = d_bnd["X1"]
    Y = d_bnd["X2"]

    # if necessary, apply fourier.real_dft_mat to get axis and boundary positions
    same_zeta = False
    if zeta_visu is None:
        zeta_out = zeta_bnd
        if zeta_bnd.shape[0] == nzeta_fp_axis:
            same_zeta = np.allclose(zeta_axis[0:nzeta_fp_axis], zeta_bnd)
    else:
        zeta_out = zeta_visu

    if same_zeta:
        # same zeta points in axis and boundary:
        pos = pos_axis[:, 0:nzeta_fp_axis]
        N = N_axis[:, 0:nzeta_fp_axis]
        B = B_axis[:, 0:nzeta_fp_axis]
    else:
        # sample axis on zeta_bnd points:
        zdft = fourier.real_dft_mat(zeta_axis, zeta_out, nfp=1)
        pos = pos_axis @ zdft["BF"].T
        N = N_axis @ zdft["BF"].T
        B = B_axis @ zdft["BF"].T

    if same_zeta:
        XX = X
        YY = Y
    else:
        # sample boundary
        zdft = fourier.real_dft_mat(zeta_bnd, zeta_out, nfp=nfp)
        XX = X @ zdft["BF"].T
        YY = Y @ zdft["BF"].T

    if theta_visu is None:
        theta_out = theta_bnd
    else:
        theta_out = theta_visu
        tdft = fourier.real_dft_mat(theta_bnd, theta_out)
        XX = tdft["BF"] @ XX
        YY = tdft["BF"] @ YY

    pos_out = pos[:, None, :] + XX[None, :, :] * N[:, None, :] + YY[None, :, :] * B[:, None, :]
    assert pos_out.ndim == 3, (
        f"problem with dimensions, pos.shape={pos.shape}, XX.shape={XX.shape}, YY.shape={YY.shape}, N.shape={N.shape}, B.shape={B.shape}"
    )
    # convert to xarray
    xr_bnd = xr.Dataset(
        coords=dict(
            theta=("pol", theta_out),
            zeta=("tor", zeta_out),
            xyz=("xyz", [0, 1, 2]),
        ),
        data_vars=dict(
            pos=(["xyz", "pol", "tor"], pos_out),
            X1=(["pol", "tor"], XX),
            X2=(["pol", "tor"], YY),
        ),
    )
    # write visualization file
    if filetype == "vts":
        ev2vtk(f"{prefix}_boundary", xr_bnd)
    elif filetype == "nc":
        xr_bnd.to_netcdf(f"{prefix}_boundary.nc", mode="w", engine="netcdf4")
