# Copyright (c) 2025 GVEC Contributors, Max Planck Institute for Plasma Physics
# License: MIT
"""gist.py - generate GIST files to be used with GENE."""

from typing import Literal
from collections.abc import Sequence, Mapping
from pathlib import Path
import logging
import argparse
import datetime

import numpy as np
import xarray as xr
import f90nml

import gvec

parser = argparse.ArgumentParser(
    prog="pygvec-to-gist",
    description="Produce a GENE-GIST input file from a GVEC state.",
)
parser.add_argument(
    "--rundir",
    type=Path,
    help="GVEC run directory",
    default=Path("."),
)
srho = parser.add_mutually_exclusive_group(required=True)
srho.add_argument(
    "-s",
    type=float,
    help="position of the target flux surface (in normalized toroidal flux, 0 < s <= 1)",
)
srho.add_argument(
    "-r",
    "--rho",
    type=float,
    help="position of the target flux surface (in square root of the normalized toroidal flux, 0 < rho <= 1)",
)
parser.add_argument(
    "--npol",
    type=float,
    default=1,
    help="number of poloidal turns (default 1)",
)
parser.add_argument(
    "--gridpoints",
    type=int,
    default=128,
    help="number of grid points along the fieldline (default 128)",
)
parser.add_argument(
    "--MNfactor",
    type=int,
    default=3,
    help="multiplication factor for the maximum fourier modes for the boozer transform (default 3)",
)
parser.add_argument(
    "-x",
    "--flip",
    choices=["auto", "none", "pol", "tor", "both"],
    default="auto",
    help="flip the poloidal or toroidal direction with respect to GVEC's Boozer coordinates; 'auto' determines the necessary flips to get positive toroidal and poloidal flux (default: 'auto')",
)
parser.add_argument(
    "-o",
    "--outputfile",
    type=Path,
    help="output file name (default: '{projectname}_s{s}.gist.txt')",
    default=None,
)
parser.add_argument(
    "-p",
    "--plot",
    action="store_true",
    help="plot the output quantities ('{projectname}_s{s}.gist.png')",
)
parser.add_argument(
    "--projectname",
    type=str,
    help="override the project name for the output files (default: use GVEC state name)",
    default=None,
)
verbosity = parser.add_mutually_exclusive_group()
verbosity.add_argument(
    "-v",
    "--verbose",
    action="count",
    default=0,
    help="verbosity level: -v for info, -vv for debug",
)
verbosity.add_argument("-q", "--quiet", action="store_true", help="suppress output")

# === Main function === #


def warning_assert_allclose(logger, name, a, b, atol=1e-16):
    try:
        np.testing.assert_allclose(a, b, atol=atol)
    except AssertionError as e:
        logger.warning(f"Assertion for {name} failed{e}")


def determine_flip(state):
    """Determine the necessary flips to get coordiantes with positive toroidal and poloidal flux."""
    ev = state.evaluate("Phi", "chi", rho=1.0).squeeze()
    if ev.Phi.item() < 0 and ev.chi.item() < 0:
        return "both"
    elif ev.Phi.item() < 0:
        return "tor"
    elif ev.chi.item() < 0:
        return "pol"
    else:
        return "none"


def generate_fieldline_coordinates(
    state, s0, gridpoints, n_pol, flip, alpha=0.0, boozer_kwargs={}
):
    rho = np.array([np.sqrt(s0)])  # radial position
    alpha = np.array([alpha])  # fieldline label
    # (poloidal) angle along the fieldline
    theta_B = np.linspace(-np.pi * n_pol, np.pi * n_pol, gridpoints, endpoint=False)

    # flip theta_B: [-pi, pi) -> [pi, -pi)
    # this does not change the coordinate, only the order of evaluation
    if flip in ("pol", "both"):
        theta_B = -theta_B

    # evaluate the rotational transform (fieldline angle) on the desired surfaces
    iota = state.evaluate("iota", rho=rho, theta=None, zeta=None).iota

    # 3D toroidal and poloidal arrays that correspond to fieldline coordinates for each surface
    zeta_B = theta_B[None, :, None] / iota.data[:, None, None] - alpha[None, None, :]

    # create the grid
    ev = gvec.EvaluationsBoozer(
        rho=rho,
        theta_B=theta_B,
        zeta_B=zeta_B,
        state=state,
        radial_derivative=True,
        **boozer_kwargs,
    )

    # set the fiedline label as toroidal coordinate & index (not necessary, but good practice)
    ev["alpha"] = ("tor", alpha)
    ev["alpha"].attrs = dict(symbol=r"\alpha", long_name="fieldline label")
    ev = ev.set_coords("alpha").set_xindex("alpha")
    return ev


def compute_gist_quantities(ev, state, flip):
    # === extract necessary quantities and convert to left-handed coordinates === #

    logger = logging.getLogger(__name__)

    state.compute(
        ev,
        "grad_rho",
        "grad_theta_B",
        "grad_zeta_B",
        "iota",
        "diota_dr",
        "p",
        "dp_dr",
        "r_minor",
        "Jac_B",
        "mod_B",
        "dmod_B_dr_B",
        "dmod_B_dt_B",
        "dmod_B_dz_B",
        "mu0",
    )

    ev = ev.squeeze()

    # === GVEC coordinates === #

    rho = ev.rho.item()
    theta = ev.theta_B
    a = ev.r_minor.item()
    iota = ev.iota.item()
    diota_dr = ev.diota_dr.item()
    p = ev.p.item()
    dp_dr = ev.dp_dr.item()
    grad_rho = ev.grad_rho
    grad_theta = ev.grad_theta_B
    grad_zeta = ev.grad_zeta_B
    dmodB_drho = ev.dmod_B_dr_B
    dmodB_dtheta = ev.dmod_B_dt_B
    dmodB_dzeta = ev.dmod_B_dz_B
    Phi_edge = state.evaluate("Phi", rho=1.0, theta=None, zeta=None).Phi.item()
    mu0 = ev.mu0.item()
    modB = ev.mod_B
    Jac_B = ev.Jac_B  # only used for assertion

    # === flipping coordinates === #

    if flip in ("tor", "both"):
        iota = -iota
        diota_dr = -diota_dr
        grad_zeta = -grad_zeta
        dmodB_dzeta = -dmodB_dzeta
        Phi_edge = -Phi_edge
        Jac_B = -Jac_B
    if flip in ("pol", "both"):
        theta = -theta
        iota = -iota
        diota_dr = -diota_dr
        grad_theta = -grad_theta
        dmodB_dtheta = -dmodB_dtheta
        Jac_B = -Jac_B

    # === fieldline coordinates === #

    grad_alpha = -theta / iota**2 * diota_dr * grad_rho + 1 / iota * grad_theta - grad_zeta
    grad_phi = grad_theta

    dmodB_drhoa = dmodB_drho - theta / iota**2 * diota_dr * dmodB_dzeta
    dmodB_dalpha = -dmodB_dzeta
    dmodB_dphi = dmodB_dtheta + 1 / iota * dmodB_dzeta

    grr = xr.dot(grad_rho, grad_rho, dims="xyz")  # g^{\rho\rho}
    gra = xr.dot(grad_rho, grad_alpha, dims="xyz")  # g^{\rho\alpha}
    grp = xr.dot(grad_rho, grad_phi, dims="xyz")  # g^{\rho\phi_\alpha}
    gap = xr.dot(grad_alpha, grad_phi, dims="xyz")  # g^{\alpha\phi_\alpha}
    gaa = xr.dot(grad_alpha, grad_alpha, dims="xyz")  # g^{\alpha\alpha}

    # === GIST quantities === #

    ngrad_x1 = a * grad_rho
    ngrad_x2 = a * rho * iota * grad_alpha
    ngrad_x3 = a * grad_phi

    g11 = xr.dot(ngrad_x1, ngrad_x1, dims="xyz")
    warning_assert_allclose(logger, "g11", g11, a**2 * grr)
    g12 = xr.dot(ngrad_x1, ngrad_x2, dims="xyz")
    warning_assert_allclose(logger, "g12", g12, a**2 * rho * iota * gra)
    g22 = xr.dot(ngrad_x2, ngrad_x2, dims="xyz")
    warning_assert_allclose(logger, "g22", g22, a**2 * rho**2 * iota**2 * gaa)
    Jac = 1 / xr.dot(ngrad_x1, xr.cross(ngrad_x2, ngrad_x3, dim="xyz"), dim="xyz")
    warning_assert_allclose(logger, "Jac", Jac, Jac_B / (rho * iota * a**3))

    Bref = 2 * Phi_edge / a**2
    L1 = -1 / (Bref * rho * iota)
    L1 *= dmodB_dalpha + dmodB_dphi * ((grr * gap - gra * grp) / (grr * gaa - gra * gra))
    L2 = 1 / Bref
    L2 *= dmodB_drhoa + dmodB_dphi * ((gaa * grp - gra * gap) / (grr * gaa - gra * gra))
    bhat = modB / Bref
    bhat_phi = dmodB_dphi / Bref
    q0 = 1 / iota
    shat = -rho / iota * diota_dr
    betahat = p * mu0 / Bref**2
    phat = -(a**4) * rho * mu0 / (2 * Phi_edge * np.abs(Phi_edge)) * dp_dr

    params = dict(
        s0=rho**2,  # radial position in normalized toroidal flux
        my_dpdx=phat,
        q0=q0,  # safety factor
        shat=shat,
        Lref=a,  # effective minor radius
        Bref=Bref,  # reference magnetic field
        n0_global=state.nfp,  # number of field periods
        beta=betahat,  # plasma beta / normalized pressure
        # gridpoints, n_pol - added later
    )
    data = np.array([g11, g12, g22, bhat, Jac, L2, -L1, bhat_phi])
    return params, data


def plot_gist(projectname: str, params: Mapping, data: np.ndarray):
    import matplotlib.pyplot as plt

    n_pol = params["n_pol"]
    gridpoints = params["gridpoints"]
    s = params["s0"]

    names = [
        r"$g^{11}$",
        r"$g^{12}$",
        r"$g^{22}$",
        r"$B/B_{ref}$",
        r"$\mathcal{J}$",
        r"$L_2$",
        r"$-L_1$",
        r"$\frac{\partial B}{\partial \theta}/B_{ref}$",
    ]
    theta = np.linspace(-np.pi * n_pol, np.pi * n_pol, gridpoints)
    fig, axs = plt.subplots(2, 4, figsize=(12, 4), layout="constrained", sharex=True)
    for ax, name, d in zip(axs.ravel(), names, data):
        ax.plot(theta, d, "k-")
        ax.set(
            title=name,
            xlabel=r"$\theta$",
            xticks=np.linspace(-np.pi * n_pol, np.pi * n_pol, 5),
            xticklabels=[f"${x:.1f}\pi$" for x in np.linspace(-n_pol, n_pol, 5)],
        )
    fig.suptitle(f"GENE-GIST output for {projectname} at $s={s:.2f}$")
    return fig


def gvec_to_gist(
    state: gvec.State,
    filename: str | Path,
    s0: float,
    gridpoints: int = 128,
    n_pol: int = 1,
    flip: Literal["auto", "none", "pol", "tor", "both"] = "auto",
    plotfile: str | Path | None = None,
    boozer_kwargs={},
):
    logger = logging.getLogger(__name__)

    # === parse arguments === #
    if flip not in ("auto", "none", "pol", "tor", "both"):
        raise ValueError("flip must be 'auto', 'none', 'pol', 'tor' or 'both'")
    if not (0 < s0 <= 1):
        raise ValueError("s0 must be in (0, 1]")
    if not (isinstance(filename, str) or isinstance(filename, Path)):
        raise ValueError("name must be a string or Path")

    if flip == "auto":
        flip = determine_flip(state)
        logger.info(f"determined flip='{flip}' for positive fluxes")
    ev = generate_fieldline_coordinates(
        state, s0, gridpoints, n_pol, flip, boozer_kwargs=boozer_kwargs
    )
    logger.info("generated fieldline coordinates")

    params, data = compute_gist_quantities(ev, state, flip)
    params["gridpoints"] = gridpoints  # number of points along fieldline / parallel resolution
    params["n_pol"] = n_pol  # number of poloidal turns
    params["gvec_version"] = gvec.__version__
    params["gvec_projectname"] = state.name  # project name of the input GVEC state
    params["gvec_datetime"] = datetime.datetime.now().isoformat(sep=" ", timespec="seconds")
    params["gvec_flip"] = flip
    logger.info("computed GIST quantities")

    if plotfile is not None:
        try:
            import matplotlib.pyplot
        except ImportError:
            logger.error("matplotlib not available, unable to generate plot")
        else:
            plot_gist(state.name, params, data).savefig(plotfile)
            logger.info(f"saved plot to '{plotfile}'")

    nml = f90nml.Namelist(dict(parameters=params))
    with open(filename, "w") as file:
        nml.write(file)
        np.savetxt(file, np.asarray(data).T, fmt="%20.10E", delimiter="")
    logger.info(f"wrote GIST file to '{filename}'")


def main(args: Sequence[str] | argparse.Namespace | None = None):
    if isinstance(args, argparse.Namespace):
        pass
    else:
        args = parser.parse_args(args)

    gvec.util.logging_setup()
    logger = logging.getLogger("gvec")
    if args.quiet:
        logging.disable()
    elif args.verbose >= 2:
        logger.setLevel(logging.DEBUG)
    elif args.verbose == 1:
        logger.setLevel(logging.INFO)
    logger.debug(f"parsed args: {args}")

    state = gvec.find_state(args.rundir)

    if args.s is not None:
        s = args.s
    else:
        s = args.rho**2

    if args.projectname is None:
        args.projectname = state.name
    if args.outputfile is None:
        args.outputfile = f"{args.projectname}_s{int(s * 100):03d}.gist.txt"

    plotfile = f"{args.projectname}_s{int(s * 100):03d}.gist.png" if args.plot else None

    gvec_to_gist(
        state,
        args.outputfile,
        s,
        gridpoints=args.gridpoints,
        n_pol=args.npol,
        flip=args.flip,
        plotfile=plotfile,
        boozer_kwargs=dict(MNfactor=args.MNfactor),
    )


if __name__ == "__main__":
    main()
