# Copyright (c) 2025 GVEC Contributors, Max Planck Institute for Plasma Physics
# License: MIT
"""cas3d.py - convert a GVEC equilibrium to be used in CAS3D"""

# === Imports === #

from pathlib import Path
import datetime
import argparse
from collections.abc import Sequence

from gvec.core import compute
import numpy as np
import xarray as xr
import tqdm

from gvec import State, EvaluationsBoozer, util, surface, __version__

# === Argument parser === #

parser = argparse.ArgumentParser(
    prog="pygvec-to-cas3d",
    description="Convert a GVEC statefile to a CAS3D compatible input file.",
)
parser.add_argument("parameterfile", type=Path, help="input GVEC parameter-file")
parser.add_argument("statefile", type=Path, help="input GVEC state-file")
parser.add_argument("outputfile", type=Path, help="output netCDF file")
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
    "--MN_booz",
    type=int,
    nargs=2,
    help="maximum fourier modes for the boozer transform (M, N)",
)
parser.add_argument(
    "--sampling",
    type=int,
    default=4,
    help="sampling factor for the fourier transform and surface reparametrization -> (S*M+1, S*N+1) points.",
)
parser.add_argument(
    "--stellsym", action="store_true", help="filter the output for stellarator symmetry"
)
parser.add_argument(
    "--pointwise",
    type=Path,
    help="output pointwise data to a separate file",
)

# === Main function === #


def gvec_to_cas3d(
    parameterfile: Path,
    statefile: Path,
    outputfile: Path,
    ns: int,
    MN_out: tuple[int, int],
    MN_booz: tuple[int, int] | None = None,
    sampling: int = 4,
    stellsym: bool = False,
    pointwise: Path | None = None,
):
    if MN_booz is None:
        MN_booz = (sampling * MN_out[0], sampling * MN_out[1])
    if sampling < 2:
        raise ValueError("sampling factor must be at least 2")

    with tqdm.tqdm(
        total=5,
        bar_format="{n_fmt}/{total_fmt} |{bar:25}| {desc}",
        desc="performing boozer transform...",
        ascii=True,
    ) as progress:
        params = util.read_parameter_file_ini(parameterfile)
        name = params["ProjectName"]

        state = State(parameterfile, statefile)
        # Boozer transform
        rho = np.sqrt(np.linspace(0, 1.0, ns))
        rho[0] = 1e-4
        ev = EvaluationsBoozer(
            rho,
            sampling * MN_out[0] + 1,
            sampling * MN_out[1] + 1,
            M=MN_booz[0],
            N=MN_booz[1],
            state=state,
        )

        # Surface reparametrization
        progress.update(1)
        progress.set_description("reparametrizing surfaces...")
        state.compute(ev, "N_FP", "pos")
        surf = surface.init_surface(ev.pos, ev.N_FP, ift="fft")
        q_surf = [
            "xhat",
            "yhat",
            "zhat",
            "g_tt_B",
            "g_tz_B",
            "g_zz_B",
            "II_tt_B",
            "II_tz_B",
            "II_zz_B",
        ]
        surface.compute(surf, *q_surf)
        surf = surf[q_surf]

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
        ]
        state.compute(ev, *q_vol)
        ev = ev[q_vol]

        ds = xr.merge([ev, surf])

        # manual conversion
        drho = 2 * ds.rho
        dtheta = -1 / (2 * np.pi)
        dzeta = 1 / (2 * np.pi)

        out = xr.Dataset()
        out.attrs = ds.attrs
        for var in ["N_FP", "mod_B", "p"]:
            out[var] = ds[var]

        # geometry
        out["xhat"] = ds.xhat
        out["yhat"] = ds.yhat
        out["zhat"] = ds.zhat
        out["Jac"] = ds.Jac_B * (drho * dtheta * dzeta) ** (-1)
        out["Jac"].attrs = dict(
            long_name="Jacobian determinant",
            symbol=r"\mathcal{J}",
            description=r"Jacobian determinant of the Boozer straight fieldline coordinates $s,\theta,\zeta$ with $s\propto\Phi$ and $,\theta,\zeta \in [0,1)$",
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
        out["Phi"].attrs = dict(
            long_name="toroidal magnetic flux",
            symbol=r"\Phi",
        )
        out["chi"] = -2 * np.pi * ds.chi
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
        out = out.assign_coords(theta_B=(2 * np.pi - out.theta_B) % (2 * np.pi))
        out = out.sortby("theta_B")

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
        ft.attrs["creator"] = "pygvec-to-cas3d"
        ft.attrs["arguments"] = repr(
            dict(ns=ns, MN_out=MN_out, MN_booz=MN_booz, sampling=sampling)
        )
        ft.attrs["statefile"] = statefile.name
        ft.attrs["state_name"] = name
        ft.attrs["conversion_time"] = (
            datetime.datetime.now().astimezone().isoformat(timespec="seconds")
        )
        ft.attrs["fourier series"] = (
            "Assumes a fourier series of the form 'v(r, θ, ζ) = Σ v_mnc(r) cos(2π m θ - 2π n N_FP ζ) + v_mns(r) sin(2π m θ - 2π n N_FP ζ)'"
        )
        ft.attrs["stellarator_symmetry"] = str(stellsym)

        # Save to netCDF
        progress.update(1)
        progress.set_description("Saving to netCDF...")
        ft.to_netcdf(outputfile)

        if pointwise is not None:
            out["s"] = out.rho**2
            out.s.attrs = dict(
                long_name="radial coordinate, normalized toroidal flux", symbol="s"
            )
            out["theta"] = out.theta_B / (2 * np.pi)  # sign flip already done above
            out.theta.attrs = dict(
                long_name="poloidal coordinate, normalized to [0,1]", symbol=r"\theta"
            )
            out["zeta"] = out.zeta_B / (2 * np.pi)
            out.zeta.attrs = dict(
                long_name="toroidal coordinate, normalized to [0,1/N_FP] for one field period",
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
                "state_name",
                "conversion_time",
            ]:
                out.attrs[key] = ft.attrs[key]

            out.to_netcdf(pointwise)

        progress.update(1)
        progress.set_description("done")


def main(args: Sequence[str] | argparse.Namespace | None = None):
    if isinstance(args, argparse.Namespace):
        pass
    else:
        args = parser.parse_args(args)

    gvec_to_cas3d(
        args.parameterfile,
        args.statefile,
        args.outputfile,
        args.ns,
        args.MN_out,
        args.MN_booz,
        args.sampling,
        args.stellsym,
        args.pointwise,
    )


if __name__ == "__main__":
    main()
