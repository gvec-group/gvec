# Copyright (c) 2025 GVEC Contributors, Max Planck Institute for Plasma Physics
# License: MIT
"""gene_gist.py - generate GIST files to be used with GENE."""

from typing import Literal
from collections.abc import Sequence, Mapping
from pathlib import Path
import logging
import argparse

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
parser.add_argument(
    "-o",
    "--outputfile",
    type=Path,
    help="output file name (default: '{projectname}_s{s}.gist.txt')",
    default=None,
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
    choices=[None, "pol", "tor", "both"],
    default=None,
    help="flip the poloidal or toroidal direction with respect (left-handed Boozer coordinates).",
)
parser.add_argument(
    "-p",
    "--plot",
    action="store_true",
    help="plot the output quantities ('{projectname}_s{s}.gist.png')",
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


def warning_assert_allclose(logger, name, a, b):
    try:
        np.testing.assert_allclose(a, b)
    except AssertionError as e:
        logger.warning(f"Assertion for {name} failed{e}")


def generate_fieldline_coordinates(state, s0, gridpoints, n_pol, flip, boozer_kwargs={}):
    rho = np.array([np.sqrt(s0)])  # radial position
    alpha = np.array([0.0])  # fieldline label
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
        "minor_radius",
        "Jac_B",
        "mod_B",
        "dmod_B_dr_B",
        "dmod_B_dt_B",
        "dmod_B_dz_B",
        "dp_dr",
        "mu0",
    )

    ev = ev.squeeze()

    rho = ev.rho
    theta = ev.theta_B
    a = ev.minor_radius.item()
    s0 = ev.rho**2
    sqrt_s0 = ev.rho
    iota = ev.iota
    diota_dr = ev.diota_dr
    diota_ds = diota_dr / (2 * ev.rho)
    p0 = ev.p.item()
    grad_rho = ev.grad_rho
    grad_vartheta = ev.grad_theta_B
    grad_zeta = ev.grad_zeta_B
    dB_drho = ev.dmod_B_dr_B
    dB_dvartheta = ev.dmod_B_dt_B
    dB_dzeta = ev.dmod_B_dz_B
    Phi_lcfs = state.evaluate("Phi", rho=1.0, theta=None, zeta=None).Phi.item()
    mu0 = ev.mu0.item()
    B = ev.mod_B
    Jac_B = ev.Jac_B  # only used for assertion

    if flip in ("tor", "both"):
        iota = -iota
        diota_dr = -diota_dr
        diota_ds = -diota_ds
        grad_zeta = -grad_zeta
        dB_dzeta = -dB_dzeta
        Phi_lcfs = -Phi_lcfs
        Jac_B = -Jac_B
    if flip in ("pol", "both"):
        theta = -theta
        iota = -iota
        diota_dr = -diota_dr
        diota_ds = -diota_ds
        grad_vartheta = -grad_vartheta
        dB_dvartheta = -dB_dvartheta
        Jac_B = -Jac_B

    # === compute GIST quantities === #

    grad_s = 2 * rho * grad_rho
    grad_alpha = -diota_dr * theta / iota**2 * grad_rho + 1 / iota * grad_vartheta - grad_zeta
    grad_theta = grad_vartheta

    q0 = 1 / iota
    dq0_ds = -diota_dr / iota**2 / (2 * rho)
    Bref = 2 * np.abs(Phi_lcfs) / a**2
    beta = p0 * mu0 / Bref**2
    sqrt_g = 1 / xr.dot(grad_s, xr.cross(grad_theta, grad_zeta, dim="xyz"), dim="xyz")
    warning_assert_allclose(logger, "sqrt_g", sqrt_g, Jac_B / (2 * ev.rho))

    dB_ds = 1 / (2 * sqrt_s0) * dB_drho - diota_ds * theta / iota**2 * dB_dzeta
    dB_dalpha = -dB_dzeta
    dB_dtheta = dB_dvartheta + 1 / iota * dB_dzeta
    dp_ds = 1 / (2 * sqrt_s0) * ev.dp_dr

    ngrad_x1 = a / (2 * sqrt_s0) * grad_s
    ngrad_x2 = a * sqrt_s0 / q0 * grad_alpha
    ngrad_x3 = a * grad_theta

    gss = xr.dot(grad_s, grad_s, dims="xyz")
    gsa = xr.dot(grad_s, grad_alpha, dims="xyz")
    gst = xr.dot(grad_s, grad_theta, dims="xyz")
    gat = xr.dot(grad_alpha, grad_theta, dims="xyz")
    gaa = xr.dot(grad_alpha, grad_alpha, dims="xyz")

    g11 = xr.dot(ngrad_x1, ngrad_x1, dims="xyz")
    warning_assert_allclose(logger, "g11", g11, a**2 / (4 * s0) * gss)
    g12 = xr.dot(ngrad_x1, ngrad_x2, dims="xyz")
    warning_assert_allclose(logger, "g12", g12, a**2 / (2 * q0) * gsa)
    g22 = xr.dot(ngrad_x2, ngrad_x2, dims="xyz")
    warning_assert_allclose(logger, "g22", g22, a**2 * s0 / q0**2 * gaa)
    Jac = 1 / xr.dot(ngrad_x1, xr.cross(ngrad_x2, ngrad_x3, dim="xyz"), dim="xyz")
    warning_assert_allclose(logger, "Jac", Jac, 2 * q0 / a**3 * sqrt_g)

    L1 = -(
        q0
        / (Bref * sqrt_s0)
        * (dB_dalpha + dB_dtheta * ((gss * gat - gsa * gst) / (gss * gaa - gsa * gsa)))
    )
    L2 = (
        2
        * sqrt_s0
        / Bref
        * (dB_ds + dB_dtheta * ((gaa * gst - gsa * gat) / (gss * gaa - gsa * gsa)))
    )
    my_dpdx = -(a**4) * rho * mu0 * dp_ds / (Phi_lcfs * np.abs(Phi_lcfs))
    dBdt = dB_dtheta / Bref
    shat = 2 * s0 / q0 * dq0_ds
    B_Bref = B / Bref

    params = dict(
        Lref=a,  # minor radius
        Bref=Bref,  # reference magnetic field
        s0=s0.item(),  # radial position in normalized toroidal flux
        q0=q0.item(),  # safety factor
        shat=shat.item(),
        my_dpdx=my_dpdx.item(),
        beta=beta,  # plasma beta / normalized pressure
        n_pol=state.nfp,  # number of field periods
        # sign_Ip_CW=sgn,  # sign of the poloidal current / poloidal magnetic field, in clockwise direction when viewed from the front
        # sign_Bt_CW=np.sign(q0.item()) * sgn,  # sign of the toroidal magnetic field, in clockwise direction when viewed from above
    )
    data = np.array([g11, g12, g22, B_Bref, Jac, L2, -L1, dBdt])
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
    flip: Literal[None, "pol", "tor", "both"] = "tor",
    plotfile: str | Path | None = None,
    boozer_kwargs={},
):
    logger = logging.getLogger(__name__)

    # === parse arguments === #
    if flip not in ("pol", "tor", "both", None):
        raise ValueError("flip must be 'pol', 'tor', 'both' or None")
    if not (0 < s0 <= 1):
        raise ValueError("s0 must be in (0, 1]")
    if not (isinstance(filename, str) or isinstance(filename, Path)):
        raise ValueError("name must be a string or Path")

    ev = generate_fieldline_coordinates(state, s0, gridpoints, n_pol, flip, boozer_kwargs)
    logger.info("generated fieldline coordinates")

    params, data = compute_gist_quantities(ev, state, flip)
    params["gridpoints"] = gridpoints  # number of points along fieldline / parallel resolution
    params["n_pol"] = n_pol  # number of poloidal turns
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
    logger = logging.getLogger(__name__)
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

    if args.outputfile is None:
        args.outputfile = args.rundir / f"{state.name}_s{int(s * 100):03d}.gist.txt"

    plotfile = f"{state.name}_s{int(s * 100):03d}.gist.png" if args.plot else None

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
