# Copyright (c) 2025 GVEC Contributors, Max Planck Institute for Plasma Physics
# License: MIT
"""cas3d.py - convert a GVEC equilibrium to be used in CAS3D"""

# === Imports === #

from pathlib import Path
import datetime
import argparse
from collections.abc import Sequence
from typing import Literal

import numpy as np
import xarray as xr
import tqdm

from gvec.core import compute
from gvec import State, EvaluationsBoozer, find_state, surface, __version__

# === Argument parser === #

parser = argparse.ArgumentParser(
    prog="pygvec-to-cas3d",
    description="Produce a CAS3D compatible input file from a GVEC state.",
)
parser.add_argument(
    "--rundir",
    type=Path,
    help="GVEC run directory",
    default=Path("."),
)
parser.add_argument(
    "--ns",
    type=int,
    help="number of flux surfaces (equally spaced in s=rho^2) (required)",
    required=True,
)
parser.add_argument(
    "--MN_out",
    type=int,
    nargs=2,
    help="maximum fourier modes in the output (M, N) (required)",
    required=True,
)
parser.add_argument(
    "--MNfactor",
    type=int,
    default=5,
    help="multiplication factor for the maximum fourier modes for the boozer transform (default 5)",
)
parser.add_argument(
    "--sampling",
    type=int,
    default=4,
    help="sampling factor for the fourier transform and surface reparametrization -> (S*M+1, S*N+1) points. (default 4)",
)
parser.add_argument(
    "--stellsym", action="store_true", help="filter the output for stellarator symmetry"
)
parser.add_argument(
    "--pointwise",
    action="store_true",
    help="output pointwise data to an additional file",
)
parser.add_argument(
    "-x",
    "--flip",
    choices=["pol", "tor"],
    default="pol",
    help="flip the poloidal or toroidal direction to match the CAS3D convention (left-handed Boozer coordinates).",
)
parser.add_argument(
    "-o",
    "--outputfile",
    type=Path,
    help="output file name (default: '{projectname}_GVEC2CAS3D')",
    default=None,
)
parser.add_argument(
    "--winding",
    type=int,
    default=1,
    help="winding number for the surface transform to (hat-coordinates) (default 1)",
)
parser.add_argument(
    "--reparam",
    action="store_true",
    help="reparametrize the surfaces in boozer angles to compute the geometric quantities",
)

# === Main function === #


def gvec_to_cas3d(
    state: State,
    outputfile: Path,
    ns: int,
    MN_out: tuple[int, int],
    MNfactor: int = 5,
    sampling: int = 2,
    stellsym: bool = False,
    pointwise: Path | None = None,
    flip: Literal["pol", "tor"] = "tor",
    winding: int = 1,
    grid: Literal["full", "half"] = "half",
    reparam: bool = False,
):
    if flip not in ["pol", "tor"]:
        raise ValueError(f"Invalid flip option: {flip}. Expected 'pol' or 'tor'.")
    with tqdm.tqdm(
        total=5,
        bar_format="{n_fmt}/{total_fmt} |{bar:25}| {desc}",
        desc="performing boozer transform...",
        ascii=True,
    ) as progress:
        # Boozer transform
        if grid == "full":
            rho = np.sqrt(np.linspace(0, 1.0, ns))
            rho[0] = 1e-4
        elif grid == "half":
            rho = np.sqrt(np.linspace(0, 1.0, 2 * ns + 1)[1::2])
        else:
            raise ValueError(f"Invalid grid option: {grid}. Expected 'full' or 'half'.")
        ev = EvaluationsBoozer(
            rho,
            sampling * MN_out[0] + 1,
            sampling * MN_out[1] + 1,
            state=state,
            MNfactor=MNfactor,
        )

        q_geo = [
            "g_tt_B",
            "g_tz_B",
            "g_zz_B",
            "II_tt_B",
            "II_tz_B",
            "II_zz_B",
        ]
        q_surf = [
            "xhat",
            "yhat",
            "zhat",
            "winding",
        ]
        # Surface reparametrization
        if reparam:
            progress.update(1)
            progress.set_description("reparametrizing surfaces...")
            state.compute(ev, "N_FP", "pos")
            surf = surface.init_surface(ev.pos, ev.N_FP, ift="fft", winding=winding)
            surface.compute(surf, *q_surf, *q_geo)
            surf = surf[q_surf + q_geo]
        else:
            progress.update(1)
            progress.set_description("computing geometric quantities...")
            state.compute(ev, "N_FP", "pos", *q_geo)
            surface.get_xyz_hat(ev, winding)

        # Quantities of interest (computed from equilibrium)
        progress.update(1)
        progress.set_description("computing equilibrium quantities...")
        q_vol = [
            "N_FP",
            "mod_B",
            "B_contra_t_B",
            "B_contra_z_B",
            "B_theta_avg",
            "B_zeta_avg",
            "iota",
            "p",
            "Phi",
            "chi",
            "Jac_B",
            "beta_avg",
        ]
        state.compute(ev, *q_vol)

        if reparam:
            ev = ev[q_vol]
            ds = xr.merge([ev, surf])
        else:
            ev = ev[q_vol + q_geo + q_surf]
            ds = ev

        # change coordinate convention: (s,θ,ζ), left-handed, with s,θ,ζ ∈ [0,1) and ζ normalized to one field period
        drho = 2 * ds.rho
        dtheta = 1 / (2 * np.pi)
        dzeta = 1 / (2 * np.pi) * ev.N_FP.item()
        if flip == "pol":
            dtheta *= -1
        elif flip == "tor":
            dzeta *= -1
            winding *= -1

        out = xr.Dataset()
        out.attrs = ds.attrs
        for var in ["N_FP", "mod_B", "p", "beta_avg", "winding", "xhat", "yhat", "zhat"]:
            out[var] = ds[var]

        # geometry
        out["Jac"] = ds.Jac_B * (drho * dtheta * dzeta) ** (-1)
        out["Jac"].attrs = dict(
            long_name="Jacobian determinant",
            symbol=r"\mathcal{J}",
            description=r"Jacobian determinant of the Boozer straight fieldline coordinates $s,\theta,\zeta$ with $s\propto\Phi$ and $s,\theta,\zeta \in [0,1)$",
        )
        out["g_tt"] = dtheta ** (-2) * ds.g_tt_B
        out["g_tt"].attrs = dict(
            long_name="poloidal component of the metric tensor",
            symbol=r"g_{\theta\theta}",
        )
        out["g_tz"] = dtheta ** (-1) * dzeta ** (-1) * ds.g_tz_B
        out["g_tz"].attrs = dict(
            long_name="poloidal-toroidal component of the metric tensor",
            symbol=r"g_{\theta\zeta}",
        )
        out["g_zz"] = dzeta ** (-2) * ds.g_zz_B
        out["g_zz"].attrs = dict(
            long_name="toroidal component of the metric tensor",
            symbol=r"g_{\zeta\zeta}",
        )
        out["II_tt"] = dtheta ** (-2) * ds.II_tt_B
        out["II_tt"].attrs = dict(
            long_name="poloidal component of the second fundamental form",
            symbol=r"\mathrm{II}_{\theta\theta}",
        )
        out["II_tz"] = dtheta ** (-1) * dzeta ** (-1) * ds.II_tz_B
        out["II_tz"].attrs = dict(
            long_name="poloidal-toroidal component of the second fundamental form",
            symbol=r"\mathrm{II}_{\theta\zeta}",
        )
        out["II_zz"] = dzeta ** (-2) * ds.II_zz_B
        out["II_zz"].attrs = dict(
            long_name="toroidal component of the second fundamental form",
            symbol=r"\mathrm{II}_{\zeta\zeta}",
        )
        # fields
        out["B_theta_avg"] = dtheta ** (-1) * ds.B_theta_avg
        out["B_theta_avg"].attrs = dict(
            long_name="covariant poloidal magnetic field",
            symbol=r"B_\theta",
        )
        out["B_zeta_avg"] = dzeta ** (-1) * ds.B_zeta_avg
        out["B_zeta_avg"].attrs = dict(
            long_name="covariant toroidal magnetic field",
            symbol=r"B_\zeta",
        )
        out["B_contra_t"] = dtheta * ds.B_contra_t_B
        out["B_contra_t"].attrs = dict(
            long_name="contravariant poloidal magnetic field",
            symbol=r"B^\theta",
        )
        out["B_contra_z"] = dzeta * ds.B_contra_z_B
        out["B_contra_z"].attrs = dict(
            long_name="contravariant toroidal magnetic field",
            symbol=r"B^\zeta",
        )
        # fluxes
        out["Phi"] = 2 * np.pi * ds.Phi
        if flip == "tor":
            out["Phi"] *= -1
        out["Phi"].attrs = dict(
            long_name="toroidal magnetic flux",
            symbol=r"\Phi",
        )
        out["chi"] = 2 * np.pi * ds.chi
        if flip == "pol":
            out["chi"] *= -1
        out["chi"].attrs = dict(
            long_name="poloidal magnetic flux",
            symbol=r"\chi",
        )
        out["iota"] = -ds.iota
        out["iota"].attrs = dict(
            long_name="rotational transform",
            symbol=r"\iota",
        )

        # flip theta
        if flip == "pol":
            out = out.assign_coords(theta_B=(-out.theta_B) % (2 * np.pi))
            out = out.sortby("theta_B")
        elif flip == "tor":
            out = out.assign_coords(zeta_B=(-out.zeta_B) % (2 * np.pi / state.nfp))
            out = out.sortby("zeta_B")

        # Fourier transform
        progress.update(1)
        progress.set_description("transforming to Fourier...")
        ft = compute.ev2ft(out)
        # Fourier truncation (remove extra modes from 'sampling' > 2)
        # Note: assumes that m=[0 ... M] and n=[0 ... N, -N ... -1]
        if sampling > 2:
            ft = ft.sel(
                m=slice(0, MN_out[0]),
                n=[*range(MN_out[1] + 1), *range(-MN_out[1], 0)],
                drop=False,
            )

        if stellsym:
            radial = [var for var in ft.data_vars if "m" not in ft[var].dims]
            odd = ["yhat", "zhat"]
            even = [var for var in out.data_vars if var not in odd and var not in radial]
            odd = [f"{var}_mns" for var in odd]
            even = [f"{var}_mnc" for var in even]
            ft = ft[radial + even + odd]
            # * xhat is even, yhat and zhat are odd-stellarator-symmetric
            # * all metric coefficients are even-stellarator-symmetric
            #   * as they are a scalar product of two basis vectors, each with a derivative of zhat in the z direction.
            #   * zhat is odd, the derivative is even, therefore the scalar product and the metric coefficients are even.
            # * modB and Jacobian are even-stellarator-symmetric
            #   * if they were odd, they would have to be 0 at the theta=zeta=0 point and flip sign there.
            #   * they both need to be > 0 everywhere though!

        out["s"] = out.rho**2
        out.s.attrs = dict(long_name="radial coordinate, normalized toroidal flux", symbol="s")
        ft = ft.swap_dims({"rad": "s"}).reset_coords("rho")

        # Set metadata
        ft.attrs["gvec_version"] = __version__
        ft.attrs["creator"] = "pygvec to-cas3d"
        ft.attrs["arguments"] = repr(
            dict(ns=ns, MN_out=MN_out, MNfactor=MNfactor, sampling=sampling, flip=flip)
        )
        ft.attrs["statefile"] = state.statefile.name
        ft.attrs["projectname"] = state.name
        ft.attrs["conversion_time"] = (
            datetime.datetime.now().astimezone().isoformat(timespec="seconds")
        )
        ft.attrs["fourier series"] = (
            "Assumes a fourier series of the form 'v(s, θ, ζ) = Σ v_mnc(s) cos(2π m θ - 2π n N_FP ζ) + v_mns(s) sin(2π m θ - 2π n N_FP ζ)'"
        )
        ft.attrs["coordinate_convention"] = (
            "Left-handed Boozer straight fieldline coordinates (s, θ, ζ), with s,θ,ζ ∈ [0,1]. "
            "s is the radial coordinate, proportional to the toroidal flux. "
            "s is spaced equidistantly on a half-mesh between 0 and 1. "
            "θ is the poloidal angle and ζ is the toroidal angle, normalized to one field period. "
        )
        ft.attrs["stellarator_symmetry"] = str(stellsym)

        # Save to netCDF
        progress.update(1)
        progress.set_description("Saving to netCDF...")
        if outputfile.exists():
            outputfile.unlink()
        ft.to_netcdf(outputfile)

        if pointwise is not None:
            out["s"] = out.rho**2
            out.s.attrs = dict(
                long_name="radial coordinate, normalized toroidal flux", symbol="s"
            )
            out["theta"] = out.theta_B * abs(dtheta)  # sign flip already done above
            out.theta.attrs = dict(
                long_name="poloidal coordinate, normalized to [0,1]", symbol=r"\theta"
            )
            out["zeta"] = out.zeta_B * abs(dzeta)
            out.zeta.attrs = dict(
                long_name="toroidal coordinate, normalized to [0,1] for one field period",
                symbol=r"\zeta",
            )
            out = (
                out.swap_dims({"rad": "s", "pol": "theta", "tor": "zeta"})
                .reset_coords("rho")
                .drop_vars(["theta_B", "zeta_B"])
            )

            # Set metadata
            for key in [
                "gvec_version",
                "creator",
                "arguments",
                "statefile",
                "projectname",
                "conversion_time",
                "coordinate_convention",
            ]:
                out.attrs[key] = ft.attrs[key]

            if pointwise.exists():
                pointwise.unlink()
            out.to_netcdf(pointwise)

        progress.update(1)
        progress.set_description("done")


def main(args: Sequence[str] | argparse.Namespace | None = None):
    if isinstance(args, argparse.Namespace):
        pass
    else:
        args = parser.parse_args(args)

    state = find_state(args.rundir)

    if args.outputfile is None:
        args.outputfile = args.rundir / f"{state.name}_GVEC2CAS3D"

    gvec_to_cas3d(
        state,
        args.outputfile.with_suffix(".nc"),
        args.ns,
        args.MN_out,
        args.MNfactor,
        args.sampling,
        args.stellsym,
        args.outputfile.parent / f"{args.outputfile.stem}_pw.nc" if args.pointwise else None,
        args.flip,
        args.winding,
        reparam=args.reparam,
    )


if __name__ == "__main__":
    main()
