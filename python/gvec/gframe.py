# Copyright (c) 2025 GVEC Contributors, Max Planck Institute for Plasma Physics
# License: MIT
r"""pyGVEC G-frame

This module provides functions for computing the G-Frame from a surface, cutting a surface at a given plane, change from xyz to xyz_hat coordinates and read and write the G-Frame from/to a file.
"""

import logging
from pathlib import Path
from typing import Literal

import numpy as np
import xarray as xr
from scipy.optimize import root_scalar

from gvec.util import logging_setup, write_parameters, linking_number, compute_FD
from gvec import fourier


def check_field_periodicity(xyz: np.ndarray, nfp: int, atol=1e-12):
    """
    checks if all  xyz positions of the surface on a full turn, ``xyz[0:nz*nfp,0:nt,0:2]``, have the field periodicity with ``nfp``.
    returns the sign of the rotation ``+2pi/nfp`` or ``-2pi/nfp``
    """
    assert xyz.shape[-1] == 3, (
        "last dimension must be the cartesian components surface positions!"
    )
    sign_rot = 1
    if nfp == 1:
        return sign_rot  # nothing to check
    assert np.mod(xyz.shape[0], nfp) == 0, (
        "number of points in zeta direction  must be divisible by nfp!"
    )
    nzeta_fp = xyz.shape[0] // nfp

    # find rotation direction:
    xyz_rot_pos = rodrigues(xyz[0, 0, :], 2 * np.pi / nfp)
    dist_pos = np.sqrt(np.sum((xyz_rot_pos - xyz[nzeta_fp, 0, :]) ** 2))
    xyz_rot_neg = rodrigues(xyz[0, 0, :], -2 * np.pi / nfp)
    dist_neg = np.sqrt(np.sum((xyz_rot_neg - xyz[nzeta_fp, 0, :]) ** 2))
    if dist_pos < atol:
        sign_rot = 1
    elif dist_neg < atol:
        sign_rot = -1
    else:
        raise ValueError(
            f"the first point of the surface [0,0] is not rotationally symmetric with the next field period. nfp={nfp}, absolute distance={np.amin([dist_pos, dist_neg])} not within tolerance {atol}!"
        )

    # rotate first fp and compare with next
    for ifp in range(1, nfp):
        xyz_rot = rodrigues(xyz[0:nzeta_fp, :, :], ifp * sign_rot * 2 * np.pi / nfp)
        maxdist = np.amax(
            np.sqrt(
                np.sum(
                    (xyz_rot - xyz[ifp * nzeta_fp : (ifp + 1) * nzeta_fp, :, :]) ** 2, axis=-1
                )
            )
        )
        if maxdist > atol:
            raise ValueError(
                f"the surface points of the first field period are not rotationally symmetric with the points in the other field periods. nfp={nfp}, maxdist={maxdist} not within tolerance {atol}!"
            )
    return sign_rot


def rodrigues(
    xyz: np.ndarray, angle, origin=np.array([0.0, 0.0, 0.0]), rot_axis=np.array([0.0, 0.0, 1.0])
):
    """
    Rodrigues rotation function.

    Parameters
    ----------
    xyz : ndarray
        cartesian point positions. if multidimensional, the LAST dimension must contain the cartesian components.
    angle : float
        The rotation angle.

    Returns
    -------
    pos_rot : ndarray
        The rotated position vector (cartesian components).
    """
    assert np.array(xyz.shape[-1]) == 3, (
        "last dimension must be the cartesian components surface positions!"
    )
    vec = xyz - origin  # origin of rotation

    cos_angle = np.cos(angle)
    sin_angle = np.sin(angle)
    cross_product = np.cross(rot_axis, vec, axis=-1)
    dot_product = np.expand_dims(np.sum(rot_axis * vec, axis=-1), axis=-1)
    vec_rot = (
        cos_angle * vec + cross_product * sin_angle + dot_product * rot_axis * (1.0 - cos_angle)
    )
    pos_rot = vec_rot + origin  # origin of rotation
    return pos_rot


def xyz_to_xyz_hat(xyz_in: np.ndarray, zeta: np.ndarray, sign_rot: float):
    """
    Change from cartesian xyz positions of the surface (can be a full torus or a single field period) on to 'hat' coordinates,
    which are periodic on field period.

    Parameters
    ----------
    xyz_in : ndarray
        xyz positions of the surface, xyz[0:nz,0:nt,0:2], sampled at zeta positions, must exclude the endpoint
    zeta : ndarray
        1d array of zeta positions belonging to the surface (without endpoint), size [0:nz]
    sign_rot : float
        direction of zeta for the rotation into hat coordinates, +1 or -1

    Returns
    -------
    xhat, yhat, zhat : ndarray
        hat coordinates, periodic on one field period, size [0:nz,0:nt,0:2]
    """

    assert sign_rot == -1 or sign_rot == 1, f"sign_rot must be -1 or 1, but is {sign_rot}"
    assert xyz_in.shape[0] == zeta.shape[0], (
        f"xyz and zeta must have the same length, but are {xyz_in.shape[0]} and {zeta.shape[0]} respectively"
    )
    sinzeta = sign_rot * np.sin(zeta)
    coszeta = np.cos(zeta)
    xhat = xyz_in[:, :, 0] * coszeta[:, None] + xyz_in[:, :, 1] * sinzeta[:, None]
    yhat = xyz_in[:, :, 1] * coszeta[:, None] - xyz_in[:, :, 0] * sinzeta[:, None]
    zhat = xyz_in[:, :, 2]
    return xhat, yhat, zhat


def xyz_hat_to_xyz(
    xhat: np.ndarray, yhat: np.ndarray, zhat: np.ndarray, zeta: np.ndarray, sign_rot: float
):
    """
    Change from xyz 'hat' coordinates to cartesian xyz positions.

    Parameters
    ----------
    xhat : ndarray
        Hat coordinates, size [0:nz,0:nt,0:2].
    yhat : ndarray
        Hat coordinates, size [0:nz,0:nt,0:2].
    zhat : ndarray
        Hat coordinates, size [0:nz,0:nt,0:2].
    zeta : ndarray
        1d array of zeta positions corresponding to xyz positions (without endpoint), size [0:nz].
    sign_rot : float
        Direction of zeta for the rotation into hat coordinates, +1 or -1.

    Returns
    -------
    xyz : ndarray
        Cartesian xyz positions of the surface, size [0:nz,0:nt,0:2].
    """
    assert sign_rot == -1 or sign_rot == 1, f"sign_rot must be -1 or 1, but is {sign_rot}"
    assert xhat.shape[0] == zeta.shape[0], (
        f"xyz and zeta must have the same length, but are {xhat.shape[0]} and {zeta.shape[0]} respectively"
    )
    assert xhat.shape == yhat.shape == zhat.shape, (
        f"xhat, yhat and zhat must have the same shape, but are {xhat.shape}, {yhat.shape} and {zhat.shape} respectively"
    )
    sinzeta = sign_rot * np.sin(zeta)
    coszeta = np.cos(zeta)
    xyz = np.zeros((xhat.shape[0], xhat.shape[1], 3))
    xyz[:, :, 0] = xhat * coszeta[:, None] - yhat * sinzeta[:, None]
    xyz[:, :, 1] = yhat * coszeta[:, None] + xhat * sinzeta[:, None]
    xyz[:, :, 2] = zhat
    return xyz


def get_X0_N_B(xyz):
    """Get guiding curve and two guiding vectors from the cartesian coordinates of a surface."""
    nt = xyz.shape[1]
    t1d = np.linspace(0, 2 * np.pi, nt, endpoint=False)
    ## STEP 2: Project to surface with elliptical cross-sections
    m0 = t1d * 0 + 1
    m1c = np.cos(t1d)
    m1s = np.sin(t1d)

    # m=0 and m=1 fourier modes. give an ellipse. theta is not necessarily the geometric angle!
    xyz0 = np.sum(xyz * m0[None, :, None], axis=1) / nt
    xyz1c = np.sum(xyz * m1c[None, :, None], axis=1) * 2 / nt
    xyz1s = np.sum(xyz * m1s[None, :, None], axis=1) * 2 / nt

    ## STEP 3: Compute the plane of the ellipse cross-sections

    N = xyz1c / np.linalg.norm(xyz1c, axis=1, keepdims=True)
    B = xyz1s / np.linalg.norm(xyz1s, axis=1, keepdims=True)

    K = np.cross(N, B, axis=1)  # not tangent of curve here, only to  ortogonalize N,B.
    K = K / np.linalg.norm(K, axis=1, keepdims=True)
    B = np.cross(K, N, axis=1)

    xyz_ell = (
        xyz0[:, None, :]
        + xyz1c[:, None, :] * m1c[None, :, None]
        + xyz1s[:, None, :] * m1s[None, :, None]
    )

    # 2D plane positions:

    x1_ell = np.sum((xyz_ell - xyz0[:, None, :]) * N[:, None, :], axis=-1)
    x2_ell = np.sum((xyz_ell - xyz0[:, None, :]) * B[:, None, :], axis=-1)

    ## STEP 4: compute fourier coefficients of the ellipse

    x1_c = np.sum(x1_ell * m1c[None, :], axis=-1) * 2 / nt
    x1_s = np.sum(x1_ell * m1s[None, :], axis=-1) * 2 / nt
    x2_c = np.sum(x2_ell * m1c[None, :], axis=-1) * 2 / nt
    x2_s = np.sum(x2_ell * m1s[None, :], axis=-1) * 2 / nt

    # deduce (rotation - theta_0):

    gam_m_thet0 = np.arctan2(x1_s - x2_c, x2_s + x1_c)

    ## STEP 5: Final frame

    Nnew = N * np.cos(-gam_m_thet0[:, None]) + B * np.sin(-gam_m_thet0[:, None])
    Bnew = -N * np.sin(-gam_m_thet0[:, None]) + B * np.cos(-gam_m_thet0[:, None])
    return xyz0, Nnew, Bnew


def eval_curve(zeta_in, xyz, dft_dict):
    """evaluates the curve at a single point zeta_in
    given by cartesian positions of a periodic curve xyz[0:len(zeta1d)+1,0:2], evaluated zeta1d[0:2pi[,
    """
    B = fourier.get_B_dft(
        np.asarray([zeta_in]).flatten(),
        deriv=dft_dict["deriv"],
        nfp=dft_dict["nfp"],
        modes=dft_dict["modes"],
    )
    return (B @ dft_dict["F"]).real @ xyz


def eval_distance_to_curve(zeta_in, origin, normal, xyz, dft_dict):
    return eval_distance_to_plane(eval_curve(zeta_in, xyz, dft_dict), origin, normal)


def eval_distance_to_plane(xyz_in, origin, normal):
    return np.sum((xyz_in - origin) * normal, axis=-1)


def find_zeta_cuts(zeta_in, origin, normal, xyz, dft_dict, zeta_bracket: float):
    def eval_dist(zeta_in):
        return eval_distance_to_curve(zeta_in, origin, normal, xyz, dft_dict).item()

    for factor in [0.1, 0.5, 0.99]:
        try:
            return root_scalar(
                eval_dist,
                bracket=[zeta_in - zeta_bracket * factor, zeta_in + zeta_bracket * factor],
                xtol=1e-10,
                rtol=1e-14,
            ).root
        except ValueError:
            pass
    raise RuntimeError(
        f"Could not find zeta cuts with the given bracket (or 1/2,1/10), for initial guess zeta_in={zeta_in}"
    )


def get_xyz_cut(zeta_start, origins, normals, xyz_in, dft_dict, nfp):
    """
    Find the xyz positions of the cut of a surface with a given number of planes.

    Parameters
    ----------
    zeta_start : float, shape (nz_out,)
        Initial guess for zeta position for each cut.
    origins : float, shape (nz_out, 3)
        Origin of the cutting plane.
    normals : float, shape (nz_out, 3)
        Normal of the cutting plane.
    dft_dict : dict
        Dictionary containing the dft from the zeta points of xyz_in, to be used for
        evaluation of theta=const curves at arbitrary zeta.

    Returns
    -------
    xyz_cut : float, shape (nz_out, nt, 3)
        XYZ positions from the cuts.
    """
    nz_out = origins.shape[0]
    nt = xyz_in.shape[1]
    zeta_out = np.zeros(nz_out)
    xyz_cut = np.zeros((nz_out, nt, 3))
    for it in range(0, nt):
        for iz in range(0, nz_out):
            zeta_out[iz] = find_zeta_cuts(
                zeta_start[iz],
                origins[iz, :],
                normals[iz, :],
                xyz_in[:, it, :],
                dft_dict,
                zeta_bracket=np.pi / nfp,
            )

        xyz_cut[:, it, :] = eval_curve(zeta_out, xyz_in[:, it, :], dft_dict)
        # check result
        assert np.allclose(eval_distance_to_plane(xyz_cut[:, it, :], origins, normals), 0)
    return xyz_cut


def cut_surf(xyz, nfp, xyz0_in, N_in, B_in):
    r"""
    Given a surface :math:`\mathbf{x}(\zeta,\theta)` on the full torus, find intersection point of lines of :math:`\theta=` const with all cutting planes defined by an origin curve ``xyz0`` and two base vectors ``N_in`` and ``B_in`` (normal of the cutting plane is then :math:`N \times B`). The intersection points are projected on the base vectors, yielding the :math:`X^1,X^2` coordinates in all cutting planes. Cutting planes must be given on the full torus and must be a multiple of the number of field periods.

    The input surface is assumed to have the same number of field periods, and must be given on the full torus. It is described by an array of cartesian point positions :math:`\mathbf{x}(\zeta_i,\vartheta_j)`, assumed to be evaluated at a meshgrid of a toroidal and poloidal parameterization :math:`(\zeta,\theta)`, excluding the periodic endpoints:

    :math:`\zeta_i=2\pi \frac{i}{n_\zeta n_{FP}},i=0\dots,(n_\zeta n_{FP})-1,\quad \theta_j=2\pi\frac{j}{n_\theta},j=0,\dots,n_\theta-1`

    The surface can thus be evaluated using a Fourier transform.

    Note: If the number of cutting planes does not match the number of toroidal points of the surface, the cutting planes are resampled to match the number of toroidal points, before cutting.

    Parameters
    ----------
        xyz :  cartesian coordinates of the surface on the full torus, excluding the periodic endpoints: shape is ``(nzeta,ntheta,3)``
        nfp : number of field periods
        xyz0_in : origin of the cutting plane, shape is ``(nzeta_cut*nfp,3)``.
        N_in : first base vector to span the cutting plane, shape is ``(nzeta_cut*nfp,3)``
        B_in : second base vector to span the cutting plane, shape is ``(nzeta_cut*nfp,3)``.

    Returns
    -------
        x1_cut,x2_cut:
            coordinates of the intersection points for cutting planes in one field period.
            shape is ``(nzeta,ntheta)``.
    """
    assert xyz.shape[2] == 3, "xyz must have shape (nzeta,ntheta,3)"
    nz = xyz.shape[0]
    assert np.mod(nz, nfp) == 0, (
        "number of surface points in zeta direction must be divisible by nfp!"
    )

    nz_gframe = xyz0_in.shape[0]
    if not xyz0_in.shape[0] == N_in.shape[0] == B_in.shape[0]:
        raise ValueError(
            "xyz0,N,B must have the same number of points, but they have different lengths!"
        )
    assert xyz0_in.shape[1] == 3, "second dimension of xyz0 must be 3"
    assert N_in.shape[1] == 3, "second dimension of N must be 3"
    assert B_in.shape[1] == 3, "second dimension of B must be 3"
    # cut geometry with new frame (xyz0,N,B)
    zeta1d = np.linspace(0.0, 2 * np.pi, nz, endpoint=False)
    if nz != nz_gframe:
        zeta_gframe = np.linspace(0.0, 2 * np.pi, nz_gframe, endpoint=False)
        zdft_gframe = fourier.real_dft_mat(zeta_gframe, zeta1d, nfp=1)
        xyz0 = zdft_gframe["BF"] @ xyz0_in
        N = zdft_gframe["BF"] @ N_in
        B = zdft_gframe["BF"] @ B_in
    else:
        xyz0 = xyz0_in
        N = N_in
        B = B_in

    zdft = fourier.real_dft_mat(zeta1d, zeta1d, nfp=1)  # must be on the full torus

    # only over one field period:
    xyz_cut = get_xyz_cut(
        zeta1d[0 : nz // nfp],
        xyz0[0 : nz // nfp, :],
        np.cross(N[0 : nz // nfp, :], B[0 : nz // nfp, :], axis=-1),
        xyz,
        zdft,
        nfp,
    )
    x1_cut = np.sum(
        (xyz_cut - xyz0[0 : nz // nfp, None, :]) * N[0 : nz // nfp, None, :], axis=-1
    )
    x2_cut = np.sum(
        (xyz_cut - xyz0[0 : nz // nfp, None, :]) * B[0 : nz // nfp, None, :], axis=-1
    )
    return x1_cut, x2_cut


def write_Gframe_ncfile(filename: str | Path, dict_in):
    """
    Write the G-Frame [and boundary] to a GVEC-compatible netCDF file.

    Parameters
    ----------
    filename : str | Path
        name of the file to be written, should include ``.nc`` ending, overwrites existing file
    dict_in : dict
        dictionary containing

        - ``nfp``: number of field periods
        - ``axis``: dictionary for the G-Frame data:

            - ``nzeta``: number of zeta positions for one field-period
            - ``zetafull``: equidistant positions along the curve parameter zeta, without endpoint, length is ``nfp*nzeta``
            - ``xyz``: cartesian positions along the curve, shape is ``(3,nfp*nzeta)``
            - ``Nxyz``: cartesian components of the 'normal' vector, shape is ``(3,nfp*nzeta)``
            - ``Bxyz``: cartesian components of the 'bi-normal' vector, shape is ``(3,nfp*nzeta)``

        - ``boundary``: dictionary for the boundary data (written to file if exists):

            - ``nzeta``: number of (toroidal ) zeta positions for one field-period
            - ``zeta``: zeta positions where boundary was sampled, length is ``nzeta``
            - ``n_max: maximum toroidal mode number ``>=(nzeta-1)//2``
            - ``ntheta``: number of (poloidal) theta positions for one field-period
            - ``theta``: theta positions where boundary was sampled, length is ``ntheta``
            - ``m_max: maximum poloidal mode number ``>=(ntheta-1)//2``
            - ``X1``: position of the boundary in the G-Frame, first component, in direction of N, shape is ``(ntheta,nzeta)``
            - ``X2``: position of the boundary in the G-Frame, second component, in direction of B, shape is ``(ntheta,nzeta)``
            - ``lasym``: logical for asymmetry, ``=False`` if boundary is stellarator-symmetric

    Returns
    -------
    None

    """
    import netCDF4 as nc

    if Path(filename).exists():
        Path(filename).unlink()

    ncfile = nc.Dataset(str(filename), "w")
    ncvars = {}
    ncfile.createDimension("vec", 3)
    ncfile.createDimension("nzeta_axis", dict_in["axis"]["nzeta"])
    nzetaFull = dict_in["axis"]["nzeta"] * dict_in["nfp"]
    assert len(dict_in["axis"]["zetafull"]) == nzetaFull, (
        f"zeta of axis must be of length nfp*nzeta={nzetaFull}, but it is {len(dict_in['axis']['zetafull'])}"
    )
    ncfile.createDimension("nzetaFull_axis", nzetaFull)
    version = 301
    axis_n_max = (dict_in["axis"]["nzeta"] - 1) // 2
    for ivar, ival in zip(
        ["VERSION", "NFP", "axis/n_max", "axis/nzeta"],
        [version, dict_in["nfp"], axis_n_max, dict_in["axis"]["nzeta"]],
    ):
        ncvars[ivar + "_var"] = ncfile.createVariable(ivar, "i8")
        ncvars[ivar + "_var"].assignValue(ival)

    if "generatedFrom" in dict_in:
        ncfile.generatedFrom = dict_in["generatedFrom"]

    ncvars["zeta_var"] = ncfile.createVariable("axis/zeta(:)", "double", ("nzeta_axis"))
    ncvars["zeta_var"][:] = dict_in["axis"]["zetafull"][0 : dict_in["axis"]["nzeta"]]

    for vecvar, vecval in zip(["axis/xyz", "axis/Nxyz", "axis/Bxyz"], ["xyz", "Nxyz", "Bxyz"]):
        assert np.all(dict_in["axis"][vecval].shape == (3, nzetaFull)), (
            f"shape of axis/{vecval} must be (3,{nzetaFull}), but it is {dict_in['axis'][vecval].shape}"
        )
        ncvars[vecvar + "_var"] = ncfile.createVariable(
            vecvar + "(::)", "f8", ("vec", "nzetaFull_axis")
        )
        ncvars[vecvar + "_var"][:, :] = dict_in["axis"][vecval]

    if "boundary" in dict_in:
        for ivar in ["ntheta", "nzeta", "m_max", "n_max", "lasym"]:
            ncvars["boundary/" + ivar + "_var"] = ncfile.createVariable(
                "boundary/" + ivar, "i8"
            )
            ncvars["boundary/" + ivar + "_var"].assignValue(1 * dict_in["boundary"][ivar])

        ncfile.createDimension("ntheta_boundary", dict_in["boundary"]["ntheta"])

        ncvars["theta_var"] = ncfile.createVariable(
            "boundary/theta(:)", "double", ("ntheta_boundary")
        )
        assert len(dict_in["boundary"]["theta"]) == dict_in["boundary"]["ntheta"]
        ncvars["theta_var"][:] = dict_in["boundary"]["theta"]

        ncfile.createDimension("nzeta_boundary", dict_in["boundary"]["nzeta"])

        ncvars["zeta_var"] = ncfile.createVariable(
            "boundary/zeta(:)", "double", ("nzeta_boundary")
        )
        assert len(dict_in["boundary"]["zeta"]) == dict_in["boundary"]["nzeta"]
        ncvars["zeta_var"][:] = dict_in["boundary"]["zeta"]

        for vecvar, vecval in zip(["boundary/X", "boundary/Y"], ["X1", "X2"]):
            assert np.all(
                dict_in["boundary"][vecval].shape
                == (dict_in["boundary"]["ntheta"], dict_in["boundary"]["nzeta"])
            ), f"shape of boundary/{vecval} must be (ntheta_boundary,nzeta_boundary)"
            ncvars[vecvar + "_var"] = ncfile.createVariable(
                vecvar + "(::)", "f8", ("ntheta_boundary", "nzeta_boundary")
            )
            ncvars[vecvar + "_var"][:, :] = dict_in["boundary"][vecval]

    ncfile.title = "== File that containts axis and boundary information, used in GVEC with the hmap_axisNB module"
    hdr = "======= HEADER OF THE NETCDF FILE VERSION 3.0.1 ==================================="
    hdr += "\n    Note: This file was generated from QUASR data using pyGVEC."
    hdr += "\n=== FILE DESCRIPTION:"
    hdr += "\n  * axis, normal and binormal of the frame are given in cartesian coordinates along the curve parameter zeta [0,2pi]."
    hdr += "\n  * The curve is allowed to have a field periodicity NFP, but the curve must be provided on a full turn."
    hdr += (
        "\n  * The data is given in real space, sampled along equidistant zeta point positions:"
    )
    hdr += "\n      zeta(i)=(i+0.5)/nzeta * (2pi/NFP), i=0,...,nzeta-1"
    hdr += "\n    always shifted by (2pi/NFP) for the next field period."
    hdr += "\n    Thus the number of points along the axis for a full turn is NFP*nzeta"
    hdr += "\n  * definition of the axis-following frame in cartesian coordinates ( boundary surface at rho=1):"
    hdr += "\n"
    hdr += "\n     {x,y,z}(rho,theta,zeta)=axis_{x,y,z}(zeta) + X(rho,theta,zeta)*N_{x,y,z}(zeta)+Y(rho,theta,zeta)*B_{x,y,z}(zeta)  "
    hdr += "\n"
    hdr += "\n=== DATA DESCRIPTION"
    hdr += "\n- general data"
    hdr += "\n  * NFP: number of field periods"
    hdr += "\n  * VERSION: version number as integer: V3.0.1 => 301"
    hdr += "\n- 'axis' data group:"
    hdr += "\n  * 'axis/n_max'   : maximum mode number in zeta (in one field period)"
    hdr += "\n  * 'axis/nzeta'   : number of points along the axis, in one field period (>=2*n_max+1)"
    hdr += "\n  * 'axis/zeta(:)' : zeta positions, 1D array of size 'axis/nzeta', for one field period. zeta[i]=zeta[1] + (i-1)/nzeta*(2pi/nfp). starting value arbitrary"
    hdr += "\n  * 'axis/xyz(::)' : cartesian positions along the axis for ONE FULL TURN, 2D array of size (3,NFP* nzeta ), sampled at zeta positions,"
    hdr += "\n                     xyz[:,j+fp*nzeta]=axis(zeta[j]+fp*2pi/NFP), for j=0,..nzeta-1 and  fp=0,...,NFP-1"
    hdr += "\n  * 'axis/Nxyz(::)': cartesian components of the normal vector of the axis frame, 2D array of size (3, NFP* nzeta), evaluated analogously to the axis"
    hdr += "\n  * 'axis/Bxyz(::)': cartesian components of the bi-normal vector of the axis frame, 2D array of size (3, NFP*nzeta), evaluated analogously to the axis"
    if "boundary" in dict_in:
        hdr += "\n- 'boundary' data group:"
        hdr += "\n  * 'boundary/m_max'    : maximum mode number in theta "
        hdr += "\n  * 'boundary/n_max'    : maximum mode number in zeta (in one field period)"
        hdr += "\n  * 'boundary/lasym'    : asymmetry, logical. "
        hdr += "\n                           if lasym=0, boundary surface position X,Y in the N-B plane of the axis frame can be represented only with"
        hdr += "\n                             X(theta,zeta)=sum X_mn*cos(m*theta-n*NFP*zeta), with {m=0,n=0...n_max},{m=1...m_max,n=-n_max...n_max}"
        hdr += "\n                             Y(theta,zeta)=sum Y_mn*sin(m*theta-n*NFP*zeta), with {m=0,n=1...n_max},{m=1...m_max,n=-n_max...n_max}"
        hdr += "\n                           if lasym=1, full fourier series is taken for X,Y"
        hdr += "\n  * 'boundary/ntheta'    : number of points in theta (>=2*m_max+1)"
        hdr += "\n  * 'boundary/nzeta'     : number of points in zeta  (>=2*n_max+1), can be different to 'axis/nzeta' !"
        hdr += "\n  * 'boundary/theta(:)'  : theta positions, 1D array of size 'boundary/ntheta',  theta[i]=theta[1] + (i-1)/ntheta*(2pi), starting value arbitrary"
        hdr += "\n  * 'boundary/zeta(:)'   : zeta positions, 1D array of size 'boundary/nzeta', for one field period! zeta[i]=zeta[1] + (i-1)/nzeta*(2pi/nfp). starting value arbitrary"
        hdr += "\n  * 'boundary/X(::)',"
        hdr += "\n    'boundary/Y(::)'     : boundary position X,Y in the N-B plane of the axis frame, in one field period, 2D array of size(ntheta, nzeta),  with"
        hdr += "\n                              X[i, j]=X(theta[i],zeta[j])"
        hdr += "\n                              Y[i, j]=Y(theta[i],zeta[j]), i=0...ntheta-1,j=0...nzeta-1"

    ncfile.header = hdr
    ncfile.close()


def read_Gframe_ncfile(ncfile: str | Path):
    """
    Read G-frame netcdf file and store data in a dictionary.

    Parameters
    ----------
    ncfile : str or Path
        Name/path to netcdf file

    Returns
    -------
    dict_out : dict
        Dictionary with the data (with 'axis' and 'boundary' groups, if they exist in the netcdf file)
    """
    with xr.open_datatree(ncfile, engine="netcdf4") as ds:
        ds.load()
    nfp = ds.NFP.data
    dict_out = {"nfp": nfp}
    if "generatedFrom" in ds:
        dict_out["generatedFrom"] = ds.generatedFrom

    if "axis" in ds:
        dict_out["axis"] = {}
        for dvar, ncvar in [
            ("zeta", "zeta(:)"),
            ("n_max", "n_max"),
            ("nzeta", "nzeta"),
            ("xyz", "xyz(::)"),
            ("Nxyz", "Nxyz(::)"),
            ("Bxyz", "Bxyz(::)"),
        ]:
            dict_out["axis"][dvar] = ds["axis"][ncvar].data
        # sizecheck
        assert dict_out["axis"]["zeta"].shape[0] == dict_out["axis"]["nzeta"], (
            "nzeta and len(zeta) must be equal!"
        )

        zeta_fp = dict_out["axis"]["zeta"]
        nzeta_fp = zeta_fp.shape[0]
        nzetaFull = dict_out["axis"]["xyz"].shape[1]
        assert nzetaFull == nfp * nzeta_fp, (
            f"axis data must be given on a full turn, with nfp being a factor in the number of points! nfp={nfp}, nzetaFull={nzetaFull}, nzeta of one fp={nzeta_fp}"
        )
        assert dict_out["axis"]["xyz"].shape == (3, nzetaFull), (
            f"shape of xyz must be (3, nzeta*nfp), but is {dict_out['axis']['xyz'].shape}"
        )
        assert dict_out["axis"]["xyz"].shape == dict_out["axis"]["Nxyz"].shape, (
            "xyz and Nxyz must have same shape"
        )
        assert dict_out["axis"]["xyz"].shape == dict_out["axis"]["Bxyz"].shape, (
            "xyz and Bxyz must have same shape"
        )

        dict_out["axis"]["nzetaFull"] = nzetaFull
        zetafull = zeta_fp[0] + np.linspace(0, 2 * np.pi, nzetaFull, endpoint=False)
        assert np.allclose(zeta_fp, zetafull[0:nzeta_fp]), "zeta on axis must be equidistant"
        dict_out["axis"]["zetafull"] = zetafull

    if "boundary" in ds:
        dict_out["boundary"] = {}
        for dvar, ncvar in [
            ("theta", "theta(:)"),
            ("zeta", "zeta(:)"),
            ("nzeta", "nzeta"),
            ("ntheta", "ntheta"),
            ("lasym", "lasym"),
            ("m_max", "m_max"),
            ("n_max", "n_max"),
            ("X1", "X(::)"),
            ("X2", "Y(::)"),
        ]:
            dict_out["boundary"][dvar] = ds["boundary"][ncvar].data

        assert dict_out["boundary"]["nzeta"] == dict_out["boundary"]["zeta"].shape[0], (
            "nzeta and len(zeta) must be equal!"
        )
        assert dict_out["boundary"]["ntheta"] == dict_out["boundary"]["theta"].shape[0], (
            "ntheta and len(theta) must be equal!"
        )
        assert dict_out["boundary"]["X1"].shape == (
            dict_out["boundary"]["ntheta"],
            dict_out["boundary"]["nzeta"],
        ), "shape of X and Y must be (ntheta, nzeta)"
        assert dict_out["boundary"]["X1"].shape == dict_out["boundary"]["X2"].shape, (
            "X and Y must have same shape"
        )

    return dict_out


def construct_gframe_from_surface(
    xyz_in: np.ndarray,
    nfp: int,
    name: str,
    tolerance_output: float = 1e-8,
    format: Literal["yaml", "toml"] = "yaml",
    tolerance_clean_surface: float = 0.0,
    impose_stell_symmetry: bool = False,
    theta0: float = 0.0,
    zeta0: float = 0.0,
    cutoff_gframe: int = -1,
    atol_field_periodicity: float = 1e-12,
    boundary_coefficients: bool = False,
    writeFiles: bool = True,
    logger: logging.Logger | None = None,
):
    """
    Construct a G-Frame from a surface given by its cartesian points.

    Parameters
    ----------
    xyz_in : ndarray, shape (nz*nfp,nt,3)
        Cartesian points xyz[0:nz*nfp,0:nt,0:2] of the surface.
    nfp : int
        Number of field periods.
    name : str
        Name of the output file.
    tolerance_output : float, optional
        Tolerance for the output surface, computes the necessary modes in X1,X2
        without changing the output data. Defaults to 1e-8.
    format : str, optional
        Output parameter file format. Defaults to "yaml".
    tolerance_clean_surface : float, optional
        Tolerance to reduce input surface resolution, computes the necessary modes
        for one field period, and recomputes the input surface with these modes.
        Defaults to 0.0.
    impose_stell_symmetry : bool, optional
        If set, imposes stellarator symmetry for the input surface. Use this
        with great care! Defaults to False.
    theta0 : float, optional
        First point in logical theta direction where xyz was sampled. Defaults to 0.
    zeta0 : float, optional
        First point in logical zeta direction where xyz was sampled. Defaults to 0.
    cutoff_gframe : int, optional
        Maximum mode number (`>=0`) to be used along the toroidal direction to
        construct the G-frame. Default `-1` means no cutoff.
    atol_field_periodicity: float, optional
        Absolute tolerance for field-periodicity check. Default is 1e-12
    boundary_coefficients : bool, optional
        Write the boundary data as fourier coefficients in the parameter file instead of pointing to the netCDF file.
    writeFiles : bool, optional
        If True, write the GVEC parameters to file and the G-frame data to netcdf file.
    logger : logging.Logger, optional
        Logger for logging messages. If None, a new logger is created.

    Returns
    -------
    parameters : dict
        Dictionary containing the GVEC parameters.
    dict_out : dict
        Dictionary containing the G-frame data.
    """
    if logger is None:
        logging_setup()
        logger = logging.getLogger(__name__)

    assert (xyz_in.shape[0] // nfp) * nfp == xyz_in.shape[0], (
        "xyz_in must be sampled on the full torus and nfp must be a factor in the number of points!"
    )
    nz_in = xyz_in.shape[0] // nfp
    nt_in = xyz_in.shape[1]
    # make the number of points odd
    nt = (nt_in // 2) * 2 + 1
    nz = (nz_in // 2) * 2 + 1
    # correct for odd numbers of points and the shift in theta
    xyz_surf = xyz_in.copy()
    if nt_in % 2 == 0 or theta0 != 0.0:
        xyz_surf = fourier.shift_1d(xyz_surf, theta0, 1, newshape=nt)
    if nz_in % 2 == 0 or zeta0 != 0.0:
        xyz_surf = fourier.shift_1d(xyz_surf, zeta0, 0, newshape=nz * nfp)

    zetafull = np.linspace(0, 2 * np.pi, nz * nfp, endpoint=False)
    theta = np.linspace(0, 2 * np.pi, nt, endpoint=False)
    # check field periodicity
    logger.info(". check field periodicity")
    sign_rot = check_field_periodicity(xyz_surf, nfp, atol_field_periodicity)

    # analyze input surface
    logger.info(". analyze input surface")
    xhat, yhat, zhat = xyz_to_xyz_hat(xyz_surf[0:nz, :, :], zetafull[0:nz], sign_rot)
    xhat_c, xhat_s = fourier.fft2d(xhat.T)
    yhat_c, yhat_s = fourier.fft2d(yhat.T)
    zhat_c, zhat_s = fourier.fft2d(zhat.T)
    Min, Nin = xhat_c.shape[0] - 1, xhat_c.shape[1] // 2
    Mmax, Nmax = Min, Nin
    recompute_xyz = False
    if tolerance_clean_surface > 0.0:
        logger.info(
            f"  - Finding minimal mode numbers for input surface with (M={Min}, N={Nin}), with tolerance {tolerance_clean_surface:.1e}"
        )
        Mmax, Nmax = minimal_modes(xhat.T, yhat.T, Z=zhat.T, tolerance=tolerance_clean_surface)
        logger.info(
            f"     Found minimal (M={Mmax}, N={Nmax}) for (xhat,yhat,zhat) in one field period."
        )

        xhat_c = fourier.scale_modes2d(xhat_c, Mmax, Nmax)
        xhat_s = fourier.scale_modes2d(xhat_s, Mmax, Nmax)
        yhat_c = fourier.scale_modes2d(yhat_c, Mmax, Nmax)
        yhat_s = fourier.scale_modes2d(yhat_s, Mmax, Nmax)
        zhat_c = fourier.scale_modes2d(zhat_c, Mmax, Nmax)
        zhat_s = fourier.scale_modes2d(zhat_s, Mmax, Nmax)
        recompute_xyz = True

    max_xhat_c = np.amax(np.abs(xhat_c))
    max_xhat_s = np.amax(np.abs(xhat_s))
    max_yhat_c = np.amax(np.abs(yhat_c))
    max_yhat_s = np.amax(np.abs(yhat_s))
    max_zhat_c = np.amax(np.abs(zhat_c))
    max_zhat_s = np.amax(np.abs(zhat_s))
    # check stellarator-symmetry of the surface:
    # is symmetric if xhat even (only cosine), yhat and zhat  odd (only sine) -> lasym =False
    norm = np.amax([np.amax(np.abs(xhat_c)), np.amax(np.abs(yhat_s)), np.amax(np.abs(zhat_s))])
    lasym = not (
        np.amax(np.abs(xhat_s)) < 1e-12 * norm
        and np.amax(np.abs(yhat_c)) < 1e-12 * norm
        and np.amax(np.abs(zhat_c)) < 1e-12 * norm
    )
    if not lasym:
        logger.info(
            "  - input surface is stellarator-symmetric (symmetric if xhat~cos,yhat~sin,zhat~sin)"
        )
    else:
        logger.info(
            "  - input surface is NOT stellarator-symmetric (symmetric if xhat~cos,yhat~sin,zhat~sin):"
        )
        logger.info(f"    max|xhat_c|={max_xhat_c}, max|xhat_s|={max_xhat_s}, ")
        logger.info(f"    max|yhat_c|={max_yhat_c}, max|yhat_s|={max_yhat_s}, ")
        logger.info(f"    max|zhat_c|={max_zhat_c}, max|zhat_s|={max_zhat_s}.")
    if impose_stell_symmetry:
        logger.info("  => Imposing stellarator-symmetry to the input surface!")
        xhat_s *= 0
        yhat_c *= 0
        zhat_c *= 0
        recompute_xyz = True

    if recompute_xyz:
        t_in = np.linspace(0, 2 * np.pi, nt_in, endpoint=False)
        z_in = np.linspace(0, 2 * np.pi, nz_in * nfp, endpoint=False)
        t, z = np.meshgrid(t_in, z_in, indexing="ij")
        xhatfull = fourier.eval2d(xhat_c, xhat_s, t, z, nfp=nfp).T
        yhatfull = fourier.eval2d(yhat_c, yhat_s, t, z, nfp=nfp).T
        zhatfull = fourier.eval2d(zhat_c, zhat_s, t, z, nfp=nfp).T
        xyz_tmp = xyz_hat_to_xyz(xhatfull, yhatfull, zhatfull, z_in, sign_rot)
        logger.info("  - maximum distance of cleaned and input surface:")
        logger.info(
            f"    max(sqrt(|xyz_old-xyz|^2))={np.amax(np.sum((xyz_in - xyz_tmp) ** 2, axis=-1) ** 0.5)}"
        )
        t, z = np.meshgrid(theta, zetafull, indexing="ij")
        xhatfull = fourier.eval2d(xhat_c, xhat_s, t, z, nfp=nfp).T
        yhatfull = fourier.eval2d(yhat_c, yhat_s, t, z, nfp=nfp).T
        zhatfull = fourier.eval2d(zhat_c, zhat_s, t, z, nfp=nfp).T
        xyz_surf = xyz_hat_to_xyz(xhatfull, yhatfull, zhatfull, zetafull, sign_rot)

    logger.info(". Constructing the G-Frame")
    # check linking number of the surface (curves at theta=0 and theta=pi):
    Lk = linking_number(xyz_surf[:, 0, :], xyz_surf[:, nt // 2, :])
    assert np.abs(Lk - np.rint(Lk)) < 1e-8, "Linking number of the surface is not integer!"
    logger.info(
        f"  - linking number of the surface = {Lk:.0f} (number of poloidal windings around its center)"
    )

    if cutoff_gframe < 0:
        nz_gframe = nz
        zetafull_gframe = np.linspace(0, 2 * np.pi, nz_gframe * nfp, endpoint=False)
        xyz_gframe = xyz_surf
    else:
        logger.info(f"  - filter surface for G-frame construction with cutoff {cutoff_gframe}")
        xhat, yhat, zhat = xyz_to_xyz_hat(xyz_surf[0:nz, :, :], zetafull[0:nz], sign_rot)
        xhat_c, xhat_s = fourier.fft2d(xhat.T)
        yhat_c, yhat_s = fourier.fft2d(yhat.T)
        zhat_c, zhat_s = fourier.fft2d(zhat.T)
        nz_gframe = max(3, 2 * cutoff_gframe + 1)
        zetafull_gframe = np.linspace(0, 2 * np.pi, nz_gframe * nfp, endpoint=False)
        xhat_c = fourier.scale_modes2d(xhat_c, Min, cutoff_gframe)
        xhat_s = fourier.scale_modes2d(xhat_s, Min, cutoff_gframe)
        yhat_c = fourier.scale_modes2d(yhat_c, Min, cutoff_gframe)
        yhat_s = fourier.scale_modes2d(yhat_s, Min, cutoff_gframe)
        zhat_c = fourier.scale_modes2d(zhat_c, Min, cutoff_gframe)
        zhat_s = fourier.scale_modes2d(zhat_s, Min, cutoff_gframe)
        t, z = np.meshgrid(theta, zetafull_gframe, indexing="ij")
        xhatfull = fourier.eval2d(xhat_c, xhat_s, t, z, nfp=nfp).T
        yhatfull = fourier.eval2d(yhat_c, yhat_s, t, z, nfp=nfp).T
        zhatfull = fourier.eval2d(zhat_c, zhat_s, t, z, nfp=nfp).T
        xyz_gframe = xyz_hat_to_xyz(xhatfull, yhatfull, zhatfull, zetafull_gframe, sign_rot)

    xyz0, N, B = get_X0_N_B(xyz_gframe)

    Wr, Lk, Tw = writhe(xyz0, N)
    logger.info(
        f"  - G-frame linking number = {Lk:.0f}, twist = {Tw:.3f}, writhe = (Lk-Tw)= {Wr:.3f}"
    )

    logger.info(". Cutting the surface")
    x1_cut, x2_cut = cut_surf(xyz_surf, nfp, xyz0, N, B)

    logger.info(
        f". Finding minimal modes for X^1,X^2, (M={Mmax}, N={Nmax}), with tolerance {tolerance_output:.1e}"
    )
    Mmax, Nmax = minimal_modes(x1_cut.T, x2_cut.T, tolerance=tolerance_output)
    logger.info(f" Found minimal (M={Mmax}, N={Nmax}) for (X^1,X^2)")

    X1c, X1s = fourier.fft2d(x1_cut.T)
    X2c, X2s = fourier.fft2d(x2_cut.T)

    X1c_minimal = fourier.scale_modes2d(X1c, Mmax, Nmax)
    X1s_minimal = fourier.scale_modes2d(X1s, Mmax, Nmax)
    X2c_minimal = fourier.scale_modes2d(X2c, Mmax, Nmax)
    X2s_minimal = fourier.scale_modes2d(X2s, Mmax, Nmax)
    t, z = np.meshgrid(theta, zetafull[0:nz], indexing="ij")
    X1_minimal = fourier.eval2d(X1c_minimal, X1s_minimal, t, z, nfp=nfp).T
    X2_minimal = fourier.eval2d(X2c_minimal, X2s_minimal, t, z, nfp=nfp).T
    logger.info(
        f"- with maximum distance in the poloidal plane: {np.amax(np.sqrt((x1_cut - X1_minimal) ** 2 + (x2_cut - X2_minimal) ** 2)):.1e}"
    )

    lasym = not (
        np.amax(np.abs(X1s)) < 1e-12 * np.amax(np.abs(X1c))
        and np.amax(np.abs(X2c)) < 1e-12 * np.amax(np.abs(X2s))
    )
    if not lasym:
        logger.info(". output X^1,X^2 coordinates are stellarator-symmetric")
    else:
        logger.info(". output X^1,X^2 coordinates are NOT stellarator-symmetric:")
        logger.info(f"  max|X1_c|={np.amax(np.abs(X1c))}, max|X1_s|={np.amax(np.abs(X1s))},")
        logger.info(f"  max|X2_c|={np.amax(np.abs(X2c))}, max|X2_s|={np.amax(np.abs(X2s))}.")
    logger.info(". Exporting h-map & boundary")
    dict_out = {"nfp": nfp}
    dict_out["generatedFrom"] = "Generated from xyz surface using following parameters:"
    dict_out["generatedFrom"] += f"tolerance_output = {tolerance_output}"
    dict_out["generatedFrom"] += f", tolerance_clean_surface = {tolerance_clean_surface}"
    dict_out["generatedFrom"] += f", impose_stell_symmetry = {impose_stell_symmetry}"
    dict_out["generatedFrom"] += f", theta0 = {theta0}"
    dict_out["generatedFrom"] += f", zeta0 = {zeta0}"
    dict_out["generatedFrom"] += f", cutoff_gframe = {cutoff_gframe}"
    dict_out["generatedFrom"] += f", atol_field_periodicity = {atol_field_periodicity}"

    dict_out["axis"] = {
        "nzeta": nz_gframe,
        "nzetaFull": nz_gframe * nfp,
        "zetafull": zetafull_gframe,
        "xyz": xyz0.T,
        "Nxyz": N.T,
        "Bxyz": B.T,
    }

    # ToDo: if not boundary_coefficients:
    dict_out["boundary"] = {
        "ntheta": nt,
        "nzeta": nz,
        "theta": theta,
        "zeta": zetafull[0:nz],
        "lasym": lasym,
        "m_max": Mmax,
        "n_max": Nmax,
        "X1": x1_cut.T,
        "X2": x2_cut.T,
    }

    parameters = dict(
        ProjectName=name,
        which_hmap=21,
        hmap_ncfile=f"{name}-Gframe.nc",
        X1X2_deg=5,
        LA_deg=5,
        sgrid=dict(
            grid_type=0,
            nElems=5,
        ),
        X1_mn_max=(Mmax, Nmax),
        X2_mn_max=(Mmax, Nmax),
        LA_mn_max=(Mmax, Nmax),
        X1_sin_cos="_cos_" if not lasym else "_sincos_",
        X2_sin_cos="_sin_" if not lasym else "_sincos_",
        LA_sin_cos="_sin_" if not lasym else "_sincos_",
        minimize_tol=1e-7,
        totalIter=100000,
        logIter=100,
        pres=dict(
            type="polynomial",
            coefs=[0.0],
        ),
        I_tor=dict(
            type="polynomial",
            coefs=[0.0],
        ),
        picard_current="auto",
    )
    if boundary_coefficients:
        if not lasym:
            parameters["X1_b_cos"] = fourier.fft2d_to_parameters(X1c_minimal)
            parameters["X2_b_sin"] = fourier.fft2d_to_parameters(X2s_minimal)
        else:
            parameters["X1_b_cos"] = fourier.fft2d_to_parameters(X1c_minimal)
            parameters["X1_b_sin"] = fourier.fft2d_to_parameters(X1s_minimal)
            parameters["X2_b_cos"] = fourier.fft2d_to_parameters(X2c_minimal)
            parameters["X2_b_sin"] = fourier.fft2d_to_parameters(X2s_minimal)
    else:
        parameters["getBoundaryFromFile"] = 1
        parameters["boundary_filename"] = f"{name}-Gframe.nc"

    if writeFiles:
        logger.info(f". Writing files: {name}-Gframe.nc , {name}-parameters.{format}")
        write_Gframe_ncfile(f"{name}-Gframe.nc", dict_out)
        write_parameters(parameters, f"{name}-parameters.{format}")
        logger.info("Done")

    return parameters, dict_out


def minimal_modes(X, Y, Z=None, tolerance=1e-8):
    """
    Find the minimal maximum mode numbers (M, N) such that the error is below the tolerance.
    First dimension of X and Y is assumed to be theta (starting at 0., without endpoint),
    second dimension is assumed to be zeta (starting at 0., without endpoint).
    """
    Xcos, Xsin = fourier.fft2d(X)
    Ycos, Ysin = fourier.fft2d(Y)
    if Z is None:
        Zcos, Zsin = np.zeros_like(Xcos), np.zeros_like(Xsin)
    else:
        Zcos, Zsin = fourier.fft2d(Z)
    M, N = Xcos.shape[0] - 1, Xcos.shape[1] // 2

    m, n = fourier.fft2d_modes(M, N, grid=True)
    Mrange, Nrange = np.arange(M + 1), np.arange(0, N + 1)
    error = np.full((M + 1, N + 1), np.inf)
    norm = np.sqrt(np.sum(Xcos**2 + Xsin**2 + Ycos**2 + Ysin**2 + Zcos**2 + Zsin**2))
    for Mnew in Mrange:
        for Nnew in Nrange:
            # sum magnitudes of all modes above the cutoff
            mask = (m > Mnew) | (n > Nnew) | (n < -Nnew)
            err = (
                Xcos[mask] ** 2
                + Xsin[mask] ** 2
                + Ycos[mask] ** 2
                + Ysin[mask] ** 2
                + Zcos[mask] ** 2
                + Zsin[mask] ** 2
            )
            error[Mnew, Nnew] = np.sqrt(np.sum(err)) / norm

    # select candidates with error below the tolerance
    mcan, ncan = np.meshgrid(Mrange, Nrange, indexing="ij")
    mask = error < tolerance
    # restrict candidates to those with minimum DoFs
    dofs = ncan + 1 + mcan * (2 * ncan + 1)
    mask &= dofs == dofs[mask].min()
    # select candidate with minimum error
    mask &= error == error[mask].min()

    return mcan[mask].item(), ncan[mask].item()


def to_axis(dict_in: dict, nzeta: int = 81):
    """
    Convert the "axis" of a gframe file to a different resolution.

    Parameters
    ----------
    dict_in : dict
        dictionary of the Gframe netcdf file, from `gvec.gframe.read_Gframe_ncfile(filename)`
    nzeta : int, optional
        number of zeta positions for the output, to sample on one field period (default: 81)

    Returns
    -------
    dict_axis : dict
        dictionary containing

        - ``nfp``: number of field periods
        - ``axis``: dictionary for the G-Frame data:

            - ``nzeta``: number of zeta positions for one field-period
            - ``zetafull``: equidistant positions along the curve parameter zeta, without endpoint, length is ``nfp*nzeta``
            - ``xyz``: cartesian positions along the curve, shape is ``(3,nfp*nzeta)``
            - ``Nxyz``: cartesian components of the 'normal' vector, shape is ``(3,nfp*nzeta)``
            - ``Bxyz``: cartesian components of the 'bi-normal' vector, shape is ``(3,nfp*nzeta)``
    """
    nfp = dict_in["nfp"]
    zetafull_out = np.linspace(0, 2 * np.pi, nzeta * nfp, endpoint=False)
    zdft = fourier.real_dft_mat(dict_in["axis"]["zetafull"], zetafull_out, nfp=1)
    xyz = zdft["BF"] @ dict_in["axis"]["xyz"].T  # [0:nz*nfp,0:2]
    Nxyz = zdft["BF"] @ dict_in["axis"]["Nxyz"].T
    Bxyz = zdft["BF"] @ dict_in["axis"]["Bxyz"].T
    return {
        "nfp": nfp,
        "axis": {
            "nzeta": nzeta,
            "zetafull": zetafull_out,
            "xyz": xyz.T,
            "Nxyz": Nxyz.T,
            "Bxyz": Bxyz.T,
        },
    }


def to_surface(dict_in: dict, nzeta: int = 81, ntheta: int = 81, tolerance: float = 1e-08):
    """
    Convert a gframe file with axis+boundary to boundary surface in cartesian coordinates.

    Parameters
    ----------
    dict_in : dict
        dictionary of the Gframe netcdf file, from `gvec.gframe.read_Gframe_ncfile(filename)`
    nzeta : int, optional
        number of zeta positions for the output, to sample on one field period (default: 81)
    ntheta : int, optional
        number of theta positions for the output (default: 81)

    Returns
    -------
    dict
        dictionary with:

        - ``xyz`` : boundary surface in cartesian coordinates, with shape ``(nzeta*nfp,ntheta,3)``
        - ``X1``, ``X2`` : boundary in G-Frame, shape ``(ntheta,nzeta)``
        - ``zetafull`` : zeta values of the boundary surface
        - ``theta`` : theta values of the boundary surface
        - ``lasym`` : logical for asymmetry, ``=False`` if stellarator symmetry is found
        - ``nfp`` : number of field periods
        - ``Mmax``, ``Nmax`` : maximum mode numbers needed for the given tolerance
        - ``X1c``, ``X1s``, ``X2c``, ``X2s`` : boundary modes in G-Frame, up to ``Mmax``, ``Nmax``
        - ``m_modes``: poloidal mode numbers (m) in first dimension of ``Rc, Rs,Zc, Zs``
        - ``n_modes``: toroidal mode numbers (n) in second dimension of ``Rc, Rs,Zc, Zs``
    """
    nfp = dict_in["nfp"]
    theta_out = np.linspace(0, 2 * np.pi, ntheta, endpoint=False)
    zetafull_out = np.linspace(0, 2 * np.pi, nzeta * nfp, endpoint=False)

    zdft = fourier.real_dft_mat(dict_in["axis"]["zetafull"], zetafull_out, nfp=1)
    xyz0 = zdft["BF"] @ dict_in["axis"]["xyz"].T  # [0:nz*nfp,0:2]
    N = zdft["BF"] @ dict_in["axis"]["Nxyz"].T
    B = zdft["BF"] @ dict_in["axis"]["Bxyz"].T

    zdft = fourier.real_dft_mat(dict_in["boundary"]["zeta"], zetafull_out, nfp=nfp)
    tdft = fourier.real_dft_mat(dict_in["boundary"]["theta"], theta_out)
    X1 = zdft["BF"] @ dict_in["boundary"]["X1"].T @ tdft["BF"].T  # [0:nzeta*nfp,0:ntheta]
    X2 = zdft["BF"] @ dict_in["boundary"]["X2"].T @ tdft["BF"].T
    xyz = xyz0[:, None, :] + X1[:, :, None] * N[:, None, :] + X2[:, :, None] * B[:, None, :]
    # transpose for output, and restrict to one field period
    X1 = X1[0:nzeta, :].T
    X2 = X2[0:nzeta, :].T
    Mmax, Nmax = minimal_modes(X1, X2, tolerance=tolerance)
    X1c, X1s = fourier.fft2d(X1)
    X2c, X2s = fourier.fft2d(X2)
    X1c = fourier.scale_modes2d(X1c, Mmax, Nmax)
    X1s = fourier.scale_modes2d(X1s, Mmax, Nmax)
    X2c = fourier.scale_modes2d(X2c, Mmax, Nmax)
    X2s = fourier.scale_modes2d(X2s, Mmax, Nmax)
    mn_modes = fourier.fft2d_modes(Mmax, Nmax)
    lasym = not (
        np.amax(np.abs(X1s)) < 1e-12 * np.amax(np.abs(X1c))
        and np.amax(np.abs(X2c)) < 1e-12 * np.amax(np.abs(X2s))
    )
    return {
        "xyz": xyz,
        "X1": X1,
        "X2": X2,
        "zetafull": zetafull_out,
        "theta": theta_out,
        "nfp": nfp,
        "lasym": lasym,
        "Mmax": Mmax,
        "Nmax": Nmax,
        "X1c": X1c,
        "X1s": X1s,
        "X2c": X2c,
        "X2s": X2s,
        "m_modes": mn_modes[0],
        "n_modes": mn_modes[1],
        "tolerance": tolerance,
    }


def to_RZ(
    xyz: np.ndarray,
    nfp: int,
    nzeta=81,
    ntheta=81,
    zeta0: float = 0.0,
    theta0: float = 0.0,
    tolerance: float = 1e-8,
):
    """
    Cut a xyz surface to yield a R,Z positions on one field period.

    Parameters
    ----------
    xyz : ndarray
        Boundary surface in cartesian coordinates, with shape ``(nzeta*nfp,ntheta,3)``
    nfp : int
        Number of field periods
    zeta0 : float, optional
        First point in logical zeta direction where xyz was sampled. Defaults to 0.
    theta0 : float, optional
        First point in logical theta direction where xyz was sampled. Defaults to 0.
    nzeta : int, optional
        Number of zeta positions (=geometric angle -phi) for the output, to sample on one field period. Defaults to 81.
    ntheta : int, optional
        Number of theta positions for the output. Defaults to 81.
    tolerance : float, optional
        Tolerance for finding minimal mode numbers. Defaults to 1e-8.

    Returns
    -------
    dict
        Dictionary with:
        - ``zeta`` : zeta positions on one field period, length ``nzeta``
        - ``theta`` : theta positions, length ``ntheta``
        - ``R`` : R positions on one field period, with shape ``(ntheta,nzeta)``
        - ``Z`` : Z positions on one field period, with shape ``(ntheta,nzeta)``
        - ``nfp`` : number of field periods
        - ``lasym`` : logical for asymmetry, ``=False`` if stellarator symmetry is found
        - ``Mmax``,``Nmax`` : maximum mode numbers needed for the given tolerance
        - ``Rc``,``Rs``,``Zc``,``Zs`` : R and Z cosine and sine Fourier mode coefficients,  shape is ``(Mmax+1,2*Nmax+1)``
        - ``m_modes``: poloidal mode numbers (m) in first dimension of ``Rc, Rs,Zc, Zs``
        - ``n_modes``: toroidal mode numbers (n) in second dimension of ``Rc, Rs,Zc, Zs``
        - ``tolerance``: input ``tolerance``
    """
    assert xyz.shape[2] == 3, "xyz must have shape [nzeta*nfp, ntheta, 3]"
    nzetafull_in, ntheta_in = xyz.shape[0], xyz.shape[1]
    nzeta_in = nzetafull_in // nfp
    assert nzeta_in * nfp == nzetafull_in, "nfp must be a factor in the number of zeta points"
    zetafull = zeta0 + np.linspace(0, 2 * np.pi, nzetafull_in, endpoint=False)
    zeta_out = zeta0 + np.linspace(0, 2 * np.pi / nfp, nzeta, endpoint=False)
    theta = theta0 + np.linspace(0, 2 * np.pi, ntheta_in, endpoint=False)
    theta_out = np.linspace(0, 2 * np.pi, ntheta, endpoint=False)
    tdft = fourier.real_dft_mat(theta, theta_out)
    zdft = fourier.real_dft_mat(zetafull, zeta_out, nfp=1)

    if theta_out.shape == theta.shape:
        xyz_t = xyz
    else:
        xyz_t = tdft["BF"][np.newaxis, :, :] @ xyz

    origins = np.zeros((nzeta, 3))
    Ncirc = np.zeros((nzeta, 3))
    Bcirc = np.zeros((nzeta, 3))
    Ncirc[:, 0] = np.cos(zeta_out)
    Ncirc[:, 1] = np.sin(-zeta_out)
    Bcirc[:, 2] = 1.0
    # try to find the rotation direction of zeta from zeta[0] and zeta[1]:
    # midpoint of zeta[0]
    xyz_mid = np.mean(xyz[0, :, :], axis=0)
    e_zeta = np.mean(xyz[1, :, :], axis=0) - xyz_mid
    sign_zeta = np.sign(
        np.cross(xyz_mid, e_zeta)[2]
    )  # +1 if zeta is counterclockwise from above

    xyz_RZcut = get_xyz_cut(
        (-sign_zeta) * zeta_out, origins, np.cross(Ncirc, Bcirc, axis=-1), xyz_t, zdft, nfp
    )
    R = np.sum((xyz_RZcut - origins[:, None, :]) * Ncirc[:, None, :], axis=-1).T
    Z = np.sum((xyz_RZcut - origins[:, None, :]) * Bcirc[:, None, :], axis=-1).T
    Mmax, Nmax = minimal_modes(R, Z, tolerance=tolerance)
    Rc, Rs = fourier.fft2d(R)
    Zc, Zs = fourier.fft2d(Z)
    Rc = fourier.scale_modes2d(Rc, Mmax, Nmax)
    Rs = fourier.scale_modes2d(Rs, Mmax, Nmax)
    Zc = fourier.scale_modes2d(Zc, Mmax, Nmax)
    Zs = fourier.scale_modes2d(Zs, Mmax, Nmax)
    mn_modes = fourier.fft2d_modes(Mmax, Nmax)
    lasym = not (
        np.amax(np.abs(Rs)) < 1e-12 * np.amax(np.abs(Rc))
        and np.amax(np.abs(Zc)) < 1e-12 * np.amax(np.abs(Zs))
    )

    return {
        "R": R,
        "Z": Z,
        "zeta": zeta_out,
        "theta": theta_out,
        "Mmax": Mmax,
        "Nmax": Nmax,
        "lasym": lasym,
        "nfp": nfp,
        "Rc": Rc,
        "Rs": Rs,
        "Zc": Zc,
        "Zs": Zs,
        "m_modes": mn_modes[0],
        "n_modes": mn_modes[1],
        "tolerance": tolerance,
    }


def plot_cross_section_comparison(dict_surf, dict_RZ, step=1, halfperiod=True):
    """
    Plot cross-sections from two dictionaries

    Parameters
    ----------
    dict_surf : dict
        Dictionary from `to_surface` function
    dict_RZ : dict
        Dictionary from `to_RZ` function
    step : int, optional
        Step in zeta array. Defaults to 1.
    halfperiod : bool, optional
        If True, only half of the zeta array is plotted. Defaults to True.

    Returns
    -------
    fig : matplotlib.figure.Figure
        Figure object
    """
    import matplotlib as mpl
    import matplotlib.pyplot as plt

    cmap = mpl.colormaps["viridis"]
    fig, axs = plt.subplots(1, 2, figsize=(14, 6))
    ax = axs[1]

    hp = 2 if halfperiod else 1
    p1 = 1 if halfperiod else 0
    nz_in = dict_surf["X1"].shape[1]
    nz = nz_in // hp + p1
    iz_s = np.arange(0, nz, step)
    c_s = hp * np.arange(nz_in) / nz_in
    for iz in iz_s:
        ax.plot(-dict_surf["X1"][:, iz], dict_surf["X2"][:, iz], color=cmap(c_s[iz]))
    ax.set_xlabel(r"$-X^1$")
    ax.set_ylabel(r"$X^2$")
    ax.set(
        title=f"N-B cross-sections, (M={dict_surf['Mmax']},N={dict_surf['Nmax']}, for tol={dict_surf['tolerance']:.0e}) "
    )
    # ax.set_aspect('equal', adjustable='box')
    ax.axis("equal")

    ax = axs[0]
    nz_in = dict_RZ["R"].shape[1]
    nz = nz_in // hp + p1
    iz_s = np.arange(0, nz, step)
    c_s = hp * np.arange(nz_in) / nz_in
    for iz in iz_s:
        ax.plot(dict_RZ["R"][:, iz], dict_RZ["Z"][:, iz], color=cmap(c_s[iz]))
    ax.set_xlabel(r"$R$")
    ax.set_ylabel(r"$Z$")
    ax.axis("equal")
    ax.set(
        title=f"R-Z cross-sections, (M={dict_RZ['Mmax']},N={dict_RZ['Nmax']}, for tol={dict_RZ['tolerance']:.0e}) "
    )
    axs[1].figure.colorbar(
        plt.cm.ScalarMappable(norm=plt.Normalize(vmin=0, vmax=1 / hp), cmap=cmap),
        ax=axs[1],
        label=r"$\zeta/(2\pi/N_{FP})$",
    )
    return fig


def twist_of_ribbon(X0: np.ndarray, N: np.ndarray, nint: int = None) -> float:
    r"""
    Compute the twist of a closed (!) ribbon defined by a centerline curve $X_0(\zeta)$ and non-vanishing "normal" vector $N(\zeta)$. The vector ``N`` only needs to be linearly independent of the tangent vector $|X_0^\prime(\zeta) \times N| >0$

    The twist is computed as

    $$
    \text{Tw} = \frac{1}{2\pi}\int_0^{2\pi}\frac{\left ({N\times \left [{N^\prime\left|{\Xp}\right|^2 - \left({N \cdot X_0^\prime}\right) X_0^{\prime\prime}}\right]}\right) \cdot X_0^\prime}{\left|{N\left|{X_0^\prime}\right|^2-\left({N \cdot X_0^\prime}\right) X_0^\prime}\right|^2}\left|{X_0^\prime}\right| d\zeta
    $$

    Derivatives of $X_0$ and $N$ are computed via fft, so the curve is assumed to be given on an equispaced grid in ``zeta=np.linspace(0,2*np.pi,npoints,endpoint=False)``, excluding the endpoint.

    Parameters
    ----------
    X0 : np.ndarray
        positions along closed curve, in cartesian coordinates, shape ``(npoints,3)``,  excluding periodic endpoint!
    N  : np.ndarray
        linearly independent, non-normalized "normal" vector at curve positions, in cartesian coordinates, same shape as ``X0``
    nint : int, optional
        number of points for integration, None sets default ``=npoints``

    Returns
    -------
    Tw : float
        twist of the ribbon defined by ``X0`` and ``N``
    """
    nzeta = X0.shape[0]
    assert nzeta == N.shape[0], (
        f"X0 and N must have the same length, but are {X0.shape[0]} and {N.shape[0]} respectively"
    )
    assert X0.shape[1] == 3, f"X0 must have shape (npoints,3), but has shape {X0.shape}"
    assert N.shape[1] == 3, f"X0 must have shape (npoints,3), but has shape {N.shape}"
    assert np.sqrt(np.sum((X0[0, :] - X0[-1, :]) ** 2)) > 1e-8 * np.sum(X0[0, :] ** 2), (
        "X0 must exclude endpoint, but first and last point coincide"
    )

    nint_zeta = nzeta
    if nint is not None and nint > nzeta:
        nint_zeta = nint
    X0_c, X0_s = fourier.fft1d(X0, axis=0)
    N_c, N_s = fourier.fft1d(N, axis=0)

    Xp = fourier.ifft1d(X0_c, X0_s, nint_zeta, deriv=1, axis=0)
    Xpp = fourier.ifft1d(X0_c, X0_s, nint_zeta, deriv=2, axis=0)
    Np = fourier.ifft1d(N_c, N_s, nint_zeta, deriv=1, axis=0)

    if nint_zeta == nzeta:
        N_f = N
    else:
        N_f = fourier.ifft1d(N_c, N_s, nint_zeta, axis=0)

    N_x_Xp = np.sqrt(np.sum((np.cross(N_f, Xp)) ** 2, axis=-1))
    N_dot_Xp = np.vecdot(N_f, Xp)
    N_Xp_angle = np.atan2(N_x_Xp, N_dot_Xp)
    assert np.all(N_Xp_angle > 1e-8), (
        f"N must be linearly independent of Xp, min angle={np.amin(N_Xp_angle)}"
    )

    # cross and vecdot apply only on last axis!
    Xp_dot_Xp = np.vecdot(Xp, Xp)

    num = np.vecdot(np.cross(N_f, (Np * Xp_dot_Xp[:, None] - N_dot_Xp[:, None] * Xpp)), Xp)
    denom = np.sum((N_f * Xp_dot_Xp[:, None] - N_dot_Xp[:, None] * Xp) ** 2, axis=-1)
    # intergate over zeta and divide by 2pi, using trapezoidal rule
    return np.average(num / denom * np.sqrt(Xp_dot_Xp))


def writhe(X0: np.ndarray, N: np.ndarray = None, width: float = 1e-6, nint: int = None):
    r"""
    Compute the writhe of a closed (!) curve defined by $X_0(\zeta)$, by attaching a ribbon defined along a non-vanishing "normal" vector $N(\zeta)$, and computing writhe as the difference of the linking number of the ribbon and the twist of the ribbon, ``Wr = Lk - Tw`.
    The vector ``N`` only needs to be linearly independent of the tangent vector $|X_0^\prime(\zeta) \times N| >0$.

    If ``N`` is not provided, the centroid frame of the curve is used.

    Note that this method is more accurate for a smooth curve than using a polygonal approximation ``writhe_from_polygon``

    Parameters
    ----------
    X0 : np.ndarray
        positions along closed curve, in cartesian coordinates, shape ``(npoints,3)``,  excluding periodic endpoint!
    N  : np.ndarray, optional
        linearly independent, non-normalized "normal" vector at curve positions, in cartesian coordinates, same shape as ``X0``. If   not provided, `N` is computed from the centroid frame of ``X0``
    width : float, optional
        In order to compute the linking number, the width of the ribbon must be specified, generating a second curve `X0+N*width`. Default width is `1e-6`
    nint : int, optional
        number of points for integration for twist, None sets default =npoints

    Returns
    -------
    Wr : float
        Write of the closed curve ``X0``
    Lk : float
        Linking number of the ribbon defined by ``X0`` and ``N``
    Tw : float
        twist of the ribbon defined by ``X0`` and ``N``
    """
    if N is None:
        N = X0 - np.mean(X0, axis=0)  # centroid frame
    Tw = twist_of_ribbon(X0, N, nint=nint)
    Lk = linking_number(X0, X0 + N * width)
    Wr = Lk - Tw
    return Wr, Lk, Tw


def frenet_frame(X0: np.ndarray) -> dict:
    r"""
    Compute the Frenet frame of a closed (!) curve defined by $X_0(\zeta)$, by computing the tangent vector $T$, normal vector $N$ and binormal vector $B$.

    Derivatives of $X_0$ are computed via fft, so the curve is assumed to be given on an equispaced grid in ``zeta=zeta0+np.linspace(0,2*np.pi,npoints,endpoint=False)``, excluding the endpoint.

    Warnings
    --------

    The curvature ``kappa`` scales with the second derivative. If it is zero at any point, normal and bi-normal vectors are not defined. ``N`` and ``B`` and torsion ``tau`` will return as ``None``.


    Parameters
    ----------
    X0 : np.ndarray
        closed curve positions, in cartesian coordinates, shape ``(npoints,3)``,  excluding periodic endpoint!

    Returns
    -------
    dict_frame : dict
        dictionary with:

        - ``X0``: input curve positions, shape ``(npoints,3)``
        - ``T``: tangent vector, shape ``(npoints,3)``
        - ``N``: normal vector, shape ``(npoints,3)``. If curvature is zero at any point, ``None`` is returned
        - ``B``: binormal vector, shape ``(npoints,3)``.  If curvature is zero at any point, ``None`` is returned
        - ``lp``: arc-length element of the curve $|X_0^{\prime}(\zeta)|$, shape ``(npoints,)``
        - ``kappa``: curvature, shape ``(npoints,)``
        - ``tau``: torsion, shape ``(npoints,)``. If curvature is zero at any point, ``None`` is returned

    """
    X0_c, X0_s = fourier.fft1d(X0, axis=0)
    Xp = fourier.ifft1d(X0_c, X0_s, axis=0, deriv=1)
    Xpp = fourier.ifft1d(X0_c, X0_s, axis=0, deriv=2)
    Xppp = fourier.ifft1d(X0_c, X0_s, axis=0, deriv=3)
    return frenet_frame_evaluate(X0, Xp, Xpp, Xppp)


def frenet_frame_evaluate(X0, X0p, X0pp, X0ppp) -> dict:
    r"""
    Compute the Frenet frame from the first three derivatives of the curve $X_0(\zeta)$ in cartesian coordinates.

    The tangent vector, binormal vector and normal vector are computed as

    .. math::
        T = \frac{X_0^\prime}{|X_0^\prime|}, \quad
        B = \frac{X_0^\prime \times X_0^{\prime\prime}}{|X_0^\prime \times X_0^{\prime\prime}|}, \quad
        N = B \times T

    The arc-length element is :math:`\ell^\prime = |X_0^\prime|`, and curvature and torsion are computed as

    .. math::

        \kappa = \frac{|X_0^\prime \times X_0^{\prime\prime}|}{(\ell^\prime)^3},\quad
        \tau = \frac{((X_0^\prime \times X_0^{\prime\prime}) \cdot X_0^{\prime\prime\prime})}{|X_0^\prime \times X_0^{\prime\prime}|^2}\,.



    Warnings
    --------

    The curvature ``kappa`` scales with the second derivative. If it is zero at any point, normal and bi-normal vectors are not defined. ``N`` and ``B`` and torsion ``tau`` will return as ``None``.

    Parameters
    ----------
    X0 : np.ndarray
        closed curve positions, in cartesian coordinates, shape ``(npoints,3)``,  excluding periodic endpoint!
    X0p : np.ndarray
        first derivative of the curve, shape ``(npoints, 3)``
    X0pp : np.ndarray
        second derivative of the curve, shape ``(npoints, 3)``
    X0ppp : np.ndarray
        third derivative of the curve, shape ``(npoints, 3)``

    Returns
    -------
    dict_frame : dict
        dictionary with:

        - ``T``: tangent vector, shape ``(npoints,3)``
        - ``N``: normal vector, shape ``(npoints,3)``. If any ``kappa*lp <1e-8`, ``None`` is returned
        - ``B``: binormal vector, shape ``(npoints,3)``.  If any ``kappa*lp <1e-8`, ``None`` is returned
        - ``lp``: arc-length element of the curve $|X_0^\prime(\zeta)|$, shape ``(npoints,)``
        - ``kappa``: curvature, shape ``(npoints,)``
        - ``tau``: torsion, shape ``(npoints,)``. If any ``kappa*lp <1e-8`, ``None`` is returned
    """
    lp = np.sqrt(np.sum(X0p**2, axis=-1))
    T = X0p / lp[:, None]
    B_full = np.cross(X0p, X0pp, axis=-1)
    Bnorm = np.sqrt(np.sum(B_full**2, axis=-1))
    kappa = Bnorm / (lp**3)
    if np.any(kappa * lp < 1e-8):
        return {"lp": lp, "kappa": kappa, "T": T, "N": None, "B": None, "tau": None}
    tau = np.sum(B_full * X0ppp, axis=-1) / Bnorm**2
    B = B_full / Bnorm[:, None]
    N = np.cross(B, T, axis=-1)
    return {"X0": X0, "T": T, "N": N, "B": B, "lp": lp, "kappa": kappa, "tau": tau}
