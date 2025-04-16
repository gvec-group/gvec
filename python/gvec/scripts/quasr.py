# Copyright (c) 2025 GVEC Contributors, Max Planck Institute for Plasma Physics
# License: MIT
"""The pyGVEC script for converting a QUASR configuration to a G-Frame for use with GVEC.

QUASR is the A QUAsi-symmetric Stellarator Repository: https://quasr.flatironinstitute.org/


### STEP 1: Evaluate surface in cartesian space:

1. using a json file from QUASR (from `https://quasr.flatironinstitute.org/`) with the simsopt interface
2. Evaluate the surface cartesian position $(x,y,z)(\vartheta_i,\zeta_j)$ at a meshgrid  on the full torus:

   $\vartheta_i=2\pi \frac{i}{n_t},i=0\dots,n_t-1,\quad \zeta_j=2\pi\frac{j}{n_z},j=0,\dots,n_z-1$

   where the angles $\vartheta,\zeta$ are just the parametrization of the given surface. (For the quasr surfaces, its a boozer angle parameterization!).

   The number $n_z$ is chosen as a multiple of the number of field periods $n_{FP}$, to be able to reduce the discrete dataset exactly to one field period.

### STEP 2: Project to surface with elliptical cross-sections

Project ${\bm x}_m(\zeta)=\frac{1}{2\pi}\int_{\vartheta=0}^{2\pi}{\bm x}(\vartheta,\zeta)\sigma_m(\vartheta)d\vartheta$ with $\sigma_0(\vartheta)=1,\sigma_{s}(\vartheta)=2\sin(\vartheta),\sigma_{c}(\vartheta)=2\cos(\vartheta)$, leading to a surface

${\bm x}_m(\vartheta,\zeta)={\bm x}_0(\zeta) + {\bm x}_s(\zeta)\sin(\vartheta)+ {\bm x}_c(\zeta)\cos(\vartheta)$

where cross-sections of $\zeta=\text{const.}$ are planar elliptical curves.


### STEP 3: Compute the plane of the ellipse cross-sections

First choice is to set the first basis unit vector $N$ from the center point of the ellipse to a point on the boundary at $\vartheta=0$ position. Then use a first guess for the second basis unit vector $B$ from the center point to $\theta=\frac{\pi}{2}$ position to span the unit normal of the plane $K=(N \times B)$, and then set the second unit vector $B=K\times N$, such that $N$ and $B$ are orthonormal and describe the plane of the ellipse.
###  STEP 4: compute fourier coefficients of the ellipse

The ellipse in a single $N,B$ plane is defined as $  $X^k(\vartheta)= x^k_c\cos(\vartheta)+x^k_s\sin(\vartheta),k=1,2$

We can deduce from the four coefficients the shift  $\vartheta_0$ and the rotation angle $\Gamma$.

### STEP 5: Final frame

The final frame is obtained by rotating the $N$ and $B$ vectors by $-(\Gamma-\vartheta_0)$, which yields constant rotation speed along $\zeta$.
Thus, in the final frame, the rotating ellipse is represented by a single poloidal and a single toroidal Fourier mode.

### STEP 6: cut the original surface with the planes of the frame

For each discrete $N,B$ plane, we compute the intersection of the all curves $\bm x(\vartheta_i,\zeta)$ and compute its position $X^1,X^2$ in the $N,B$ plane. This gives the final surface.
"""

import numpy as np
from sympy import Q
from scipy.optimize import root_scalar

# === Functions === #


def real_dft_mat(x_in, x_out, nfp=1, modes=None, deriv=0):
    """
    Flexible Direct Fourier Transform for real data
    takes an input array of equidistant points in [0,2pi/nfp[ (exclude endpoint!),
    evaluate the discrete fourier transform with the given 1d mode vector (all >=0) using the input points x_in, then evaluate the inverse transform (or its derivative deriv>0) on the output points x_out anywhere...
    len(x_in) must be > 2*max(modes)
    output is the matrix that transforms real function to real function [derivative]:
     f^deriv(x_out) = Mat f(x_in) (can then be used to do 2d transforms with matmul!)

    nfp is the number of field periods, default 1 (int), all modes are multiples of nfp

    """
    if modes is None:
        modes = np.arange((len(x_in) - 1) // 2 + 1)  # all modes up to Nyquist
    assert np.allclose(x_in[-1] + (x_in[1] - x_in[0]) - x_in[0], 2 * np.pi / nfp)
    assert np.all(modes >= 0), "modes must be positive"
    zeromode = np.where(modes == 0)
    assert len(zeromode) <= 1, "only one zero mode allowed"
    maxmode = np.amax(modes)
    assert len(x_in) > 2 * maxmode, (
        f"number of sampling points ({len(x_in)}) > 2*maxmodenumber ({maxmode})"
    )
    # matrix for forward transform
    Fmat = np.exp(1j * nfp * (modes[:, None] * x_in[None, :]))
    mass_re = Fmat.real @ Fmat.real.T
    mass_im = Fmat.imag @ Fmat.imag.T
    diag_re = np.copy(np.diag(mass_re))
    diag_im = np.copy(np.diag(mass_im))

    assert np.all(np.abs(mass_re - np.diag(diag_re)) < 1.0e-8), (
        "massre must be diagonal"
    )
    assert np.all(np.abs(mass_im - np.diag(diag_im)) < 1.0e-8), (
        "massim must be diagonal"
    )
    diag_im[zeromode] = 1  # imag (=sin) is zero at zero mode
    assert np.all(diag_re > 0.0)
    assert np.all(diag_im > 0.0)

    # inverse mass matrix applied (for real and imag)
    Fmat_mod = np.diag(1 / diag_re) @ Fmat.real + np.diag(1j / diag_im) @ Fmat.imag
    Bmat = get_B(x_out=x_out, nfp=nfp, modes=modes, deriv=deriv)
    Mat = (Bmat @ Fmat_mod).real
    return {
        "F": Fmat_mod,
        "B": Bmat,
        "BF": Mat,
        "modes": modes,
        "x_in": x_in,
        "x_out": x_out,
        "deriv": deriv,
        "nfp": nfp,
    }


def get_B(x_out, deriv, nfp, modes):
    modes_back = np.exp(-1j * nfp * (modes[None, :] * x_out[:, None]))
    if deriv > 0:
        modes_back *= (-1j * nfp * modes[None, :]) ** deriv
    return modes_back


def get_xyz_from_quasr(nt_in, nz_in, quasr_json):
    """
    load surface object using simsopt, sample surface at nt_in,nz_in*nfp point positions on the full torus.
    gives cartesian positions xyz[0:nz_in*nfp,0:nt_in,0:2]
    """
    from simsopt._core import load

    [surfaces, coils] = load(quasr_json + ".json")
    surf = surfaces[-1]
    nfp = surf.nfp
    nt = nt_in
    nz = nz_in * nfp
    t1d = np.linspace(0, 1, nt, endpoint=False)
    z1d = np.linspace(0, 1.0, nz, endpoint=False)
    t, z = np.meshgrid(t1d, z1d)

    xyz = np.zeros((nz, nt, 3))
    surf.gamma_lin(xyz, z.flatten(), t.flatten())
    return {"t1d": t1d, "z1d": z1d, "nfp": nfp, "nt": nt, "nz": nz, "xyz": xyz}


def get_X0_N_B(xyz):
    nt = xyz.shape[1]
    t1d = np.linspace(0, 1, nt, endpoint=False)
    ## STEP 2: Project to surface with elliptical cross-sections
    m0 = t1d * 0 + 1
    m1c = np.cos(2 * np.pi * t1d)
    m1s = np.sin(2 * np.pi * t1d)

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
    B = get_B(np.asarray([zeta_in]).flatten(), **dft_dict)
    return (B @ dft_dict["F"]).real @ xyz


def eval_distance_to_curve(zeta_in, origin, normal, xyz, dft_dict):
    return eval_distance_to_plane(eval_curve(zeta_in, xyz, dft_dict), origin, normal)


def eval_distance_to_plane(xyz_in, origin, normal):
    return np.sum((xyz_in - origin) * normal, axis=-1)


def find_zeta_cuts(zeta_in, origin, normal, xyz, dft_dict):
    def eval_dist(zeta_in):
        return eval_distance_to_curve(zeta_in, origin, normal, xyz, dft_dict)

    return root_scalar(
        eval_dist, bracket=[zeta_in - 0.2 * np.pi, zeta_in + 0.2 * np.pi], xtol=1e-10
    ).root


def get_xyz_cut(zeta_start, origins, normals, xyz_in, dft_dict):
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
            )

        xyz_cut[:, it, :] = eval_curve(zeta_out, xyz_in[:, it, :], dft_dict)
        # check result
        assert np.allclose(
            eval_distance_to_plane(xyz_cut[:, it, :], origins, normals), 0
        )
    return xyz_cut


def cut_surf(xyz, nfp, xyz0, N, B):
    """
    given xyz(zeta,theta) on the full torus, find intersection point of lines of theta=const with all N-B planes with origin xyz0. then project these points to find x1,x2 coordinates in each N-B cross-section
    """
    nz = xyz.shape[0]
    if not nz == xyz0.shape[0] == N.shape[0] == B.shape[0]:
        raise ValueError(
            "xyz0,N,B must have the same number of points, but they have different lengths!"
        )
    # cut geometry with new frame (xyz0,N,B)
    zeta1d = np.linspace(0.0, 2 * np.pi, nz, endpoint=False)
    zdft = real_dft_mat(zeta1d, zeta1d, nfp=1)  # must be on the full torus

    # only over one field period:
    xyz_cut = get_xyz_cut(
        zeta1d[0 : nz // nfp],
        xyz0[0 : nz // nfp, :],
        np.cross(N[0 : nz // nfp, :], B[0 : nz // nfp, :], axis=-1),
        xyz,
        zdft,
    )
    x1_cut = np.sum(
        (xyz_cut - xyz0[0 : nz // nfp, None, :]) * N[0 : nz // nfp, None, :], axis=-1
    )
    x2_cut = np.sum(
        (xyz_cut - xyz0[0 : nz // nfp, None, :]) * B[0 : nz // nfp, None, :], axis=-1
    )
    return x1_cut, x2_cut


def write_Gframe_ncfile(ncout_file, dict_in):
    """
    write the data from dictionary to netcdf
    """
    import netCDF4 as nc

    ncfile = nc.Dataset(f"{ncout_file}.nc", "w")
    ncvars = {}
    ncfile.createDimension("vec", 3)
    ncfile.createDimension("nzeta_axis", dict_in["axis"]["nzeta"])
    assert (
        len(dict_in["axis"]["zetafull"]) == dict_in["axis"]["nzeta"] * dict_in["nfp"]
    ), "zeta of axis must be of length nfp*nzeta!"
    ncfile.createDimension("nzetaFull_axis", dict_in["axis"]["nzeta"] * dict_in["nfp"])
    version = 300
    axis_n_max = (
        dict_in["axis"]["nzeta"] * dict_in["axis"]["nzeta"] * dict_in["nfp"] - 1
    ) // 2
    for ivar, ival in zip(
        ["VERSION", "NFP", "axis/n_max", "axis/nzeta"],
        [version, dict_in["nfp"], axis_n_max, dict_in["axis"]["nzeta"]],
    ):
        ncvars[ivar + "_var"] = ncfile.createVariable(ivar, "i8")
        ncvars[ivar + "_var"].assignValue(ival)

    ncvars["zeta_var"] = ncfile.createVariable("axis/zeta(:)", "double", ("nzeta_axis"))
    ncvars["zeta_var"][:] = dict_in["axis"]["zetafull"][0 : dict_in["axis"]["nzeta"]]

    for vecvar, vecval in zip(
        ["axis/xyz", "axis/Nxyz", "axis/Bxyz"], ["xyz", "Nxyz", "Bxyz"]
    ):
        assert np.all(
            dict_in["axis"][vecval].shape == (3, dict_in["axis"]["nzetaFull"])
        ), (
            f"shape of axis/{vecval} must be (3,nzetaFull_axis), but it is {dict_in['axis'][vecval].shape}"
        )
        ncvars[vecvar + "_var"] = ncfile.createVariable(
            vecvar + "(::)", "f8", ("vec", "nzetaFull_axis")
        )
        ncvars[vecvar + "_var"][:, :] = dict_in["axis"][vecval]

    boundary_ntheta = dict_in["boundary"]["ntheta"]
    boundary_m_max = (boundary_ntheta - 1) // 2
    boundary_nzeta = dict_in["boundary"]["nzeta"]  # same for now
    boundary_n_max = (boundary_nzeta - 1) // 2

    boundary_lasym = 1 * dict_in["boundary"]["lasym"]

    for ivar, ival in zip(
        [
            "boundary/ntheta",
            "boundary/nzeta",
            "boundary/m_max",
            "boundary/n_max",
            "boundary/lasym",
        ],
        [
            boundary_ntheta,
            boundary_nzeta,
            boundary_m_max,
            boundary_n_max,
            boundary_lasym,
        ],
    ):
        ncvars[ivar + "_var"] = ncfile.createVariable(ivar, "i8")
        ncvars[ivar + "_var"].assignValue(ival)

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
            dict_in["boundary"][vecval].shape == (boundary_ntheta, boundary_nzeta)
        ), f"shape of boundary/{vecval} must be (ntheta_boundary,nzeta_boundary)"
        ncvars[vecvar + "_var"] = ncfile.createVariable(
            vecvar + "(::)", "f8", ("ntheta_boundary", "nzeta_boundary")
        )
        ncvars[vecvar + "_var"][:, :] = dict_in["boundary"][vecval]

    ncfile.title = "== File that containts axis and boundary information, used in GVEC with the hmap_axisNB module"
    hdr = "======= HEADER OF THE NETCDF FILE VERSION 3.0 ==================================="
    hdr += "\n    Note: This file was generated from QUASR data using pyGVEC."
    hdr += "\n=== FILE DESCRIPTION:"
    hdr += "\n  * axis, normal and binormal of the frame are given in cartesian coordinates along the curve parameter zeta [0,2pi]."
    hdr += "\n  * The curve is allowed to have a field periodicity NFP, but the curve must be provided on a full turn."
    hdr += "\n  * The data is given in real space, sampled along equidistant zeta point positions:"
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
    hdr += "\n  * VERSION: version number as integer: V3.0 => 300"
    hdr += "\n- 'axis' data group:"
    hdr += "\n  * 'axis/n_max'   : maximum mode number in zeta (in one field period)"
    hdr += "\n  * 'axis/nzeta'   : number of points along the axis, in one field period (>=2*n_max+1)"
    hdr += "\n  * 'axis/zeta(:)' : zeta positions, 1D array of size 'axis/nzeta', for one field period. zeta[i]=zeta[1] + (i-1)/nzeta*(2pi/nfp). starting value arbitrary"
    hdr += "\n  * 'axis/xyz(::)' : cartesian positions along the axis for ONE FULL TURN, 2D array of size (3,NFP* nzeta ), sampled at zeta positions,"
    hdr += "\n                     xyz[:,j+fp*nzeta]=axis(zeta[j]+fp*2pi/NFP), for j=0,..nzeta-1 and  fp=0,...,NFP-1"
    hdr += "\n  * 'axis/Nxyz(::)': cartesian components of the normal vector of the axis frame, 2D array of size (3, NFP* nzeta), evaluated analogously to the axis"
    hdr += "\n  * 'axis/Bxyz(::)': cartesian components of the bi-normal vector of the axis frame, 2D array of size (3, NFP*nzeta), evaluated analogously to the axis"
    hdr += "\n- 'boundary' data group:"
    hdr += "\n  * 'boundary/m_max'    : maximum mode number in theta "
    hdr += (
        "\n  * 'boundary/n_max'    : maximum mode number in zeta (in one field period)"
    )
    hdr += "\n  * 'boundary/lasym'    : asymmetry, logical. "
    hdr += "\n                           if lasym=0, boundary surface position X,Y in the N-B plane of the axis frame can be represented only with"
    hdr += "\n                             X(theta,zeta)=sum X_mn*cos(m*theta-n*NFP*zeta), with {m=0,n=0...n_max},{m=1...m_max,n=-n_max...n_max}"
    hdr += "\n                             Y(theta,zeta)=sum Y_mn*sin(m*theta-n*NFP*zeta), with {m=0,n=1...n_max},{m=1...m_max,n=-n_max...n_max}"
    hdr += (
        "\n                           if lasym=1, full fourier series is taken for X,Y"
    )
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

    print(('NETCDF FILE "%s" WRITTEN! ' % (ncout_file + ".nc")))


def main():
    nt_in = 81
    nz_in = 81
    # quasr_json="quasr-QH-N3-serial1942772"
    # quasr_json="serial1601280"
    # quasr_json="serial1601292"
    # quasr_json="serial1603880"
    quasr_json = "serial2457904"

    surf = get_xyz_from_quasr(nt_in, nz_in, quasr_json)
    # nz = surf["nz"]
    nfp = surf["nfp"]
    xyz = surf["xyz"]
    print("get frame")
    xyz0, N, B = get_X0_N_B(xyz)
    print("cut surface")
    x1_cut, x2_cut = cut_surf(xyz, nfp, xyz0, N, B)

    zetafull = np.linspace(0, 2 * np.pi, nz_in * nfp, endpoint=False)
    zeta = np.linspace(0, 2 * np.pi / nfp, nz_in, endpoint=False)
    theta = np.linspace(0, 2 * np.pi, nt_in, endpoint=False)
    dict_out = {"nfp": nfp, "axis": {}, "boundary": {}}
    dict_out["axis"] = {
        "nzeta": nz_in,
        "nzetaFull": nz_in * nfp,
        "zetafull": zetafull,
        "xyz": xyz0.T,
        "Nxyz": N.T,
        "Bxyz": B.T,
    }
    dict_out["boundary"] = {
        "ntheta": nt_in,
        "nzeta": nz_in,
        "theta": theta,
        "zeta": zeta,
        "lasym": False,
        "X1": x1_cut.T,
        "X2": x2_cut.T,
    }
    write_Gframe_ncfile(quasr_json + "_gvec", dict_out)
