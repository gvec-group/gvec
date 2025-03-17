# Copyright (c) 2025 GVEC Contributors, Max Planck Institute for Plasma Physics
# License: MIT
"""The pyGVEC run script for running GVEC using stages and current constraints."""

import argparse
import copy
import datetime
import logging
import re
import shutil
import time
from pathlib import Path
from typing import Mapping

import numpy as np
import xarray as xr

import gvec

# === Argument Parser === #

parser = argparse.ArgumentParser(
    prog="pygvec-run",
    formatter_class=argparse.RawDescriptionHelpFormatter,
    description="Run GVEC with a given parameterfile, optionally restarting from an existing statefile.\n\n"
    "When given an INI parameterfile, GVEC is called directly.\n"
    "With YAML and TOML parameterfiles, GVEC can be run in several stages and a current constraint with picard iterations can be performed.",
)
parser.add_argument("parameterfile", type=Path, help="input GVEC parameterfile")
parser.add_argument(
    "restartfile",
    type=Path,
    help="GVEC statefile to restart from (optional)",
    nargs="?",
)
param_type = parser.add_mutually_exclusive_group()
param_type.add_argument(
    "--ini",
    action="store_const",
    const="ini",
    dest="param_type",
    help="interpret GVEC parameterfile classicly (INI)",
)
param_type.add_argument(
    "--yaml",
    action="store_const",
    const="yaml",
    dest="param_type",
    help="interpret GVEC parameterfile as YAML",
)
param_type.add_argument(
    "--toml",
    action="store_const",
    const="toml",
    dest="param_type",
    help="interpret GVEC parameterfile as TOML",
)
verbosity = parser.add_mutually_exclusive_group()
verbosity.add_argument(
    "-v",
    "--verbose",
    action="count",
    default=0,
    help="verbosity level: -v for info, -vv for debug, -vvv for GVEC output",
)
verbosity.add_argument("-q", "--quiet", action="store_true", help="suppress output")
parser.add_argument(
    "-d",
    "--diagnostics",
    type=Path,
    default=None,
    help="output netCDF file for diagnostics",
)
parser.add_argument("-p", "--plots", action="store_true", help="plot diagnostics")

# === Script === #


def run_stages(
    parameters: Mapping,
    statefile: Path | None = None,
    progressbar: bool = False,
    redirect_gvec_stdout: bool = True,
    diagnosticfile: Path | None = None,
    plots: bool = False,
) -> tuple[Path, Path, xr.Dataset]:
    """Run GVEC with several stages (assuming hierarchical parameters)"""
    logger = logging.getLogger("pyGVEC.script")
    diagnostics: xr.Dataset | None = None
    rho = np.sqrt(np.linspace(0, 1, 101))
    rho[0] = 1e-4

    if "Itor" in parameters:
        match parameters["Itor"].get("type", "polynomial"):
            case "polynomial":
                coefs = np.array(parameters["Itor"]["coefs"][::-1])
                coefs *= parameters["Itor"].get("scale", 1.0)
                I_tor_target = np.poly1d(coefs)(rho**2)
            case "bspline":
                from scipy.interpolate import BSpline

                coefs = np.array(parameters["Itor"]["coefs"], dtype=float)
                coefs *= parameters["Itor"].get("scale", 1.0)
                knots = np.array(parameters["Itor"]["knots"], dtype=float)
                I_tor_target = BSpline(knots, coefs)(rho**2)
            case "interpolation":
                I_tor_target = np.array(parameters["Itor"]["vals"], dtype=float)
                rho = np.sqrt(np.array(parameters["Itor"]["rho2"], dtype=float))
            case _:
                raise ValueError(f"Unknown Itor type: {parameters['Itor']['type']}")

    for s, stage in enumerate(parameters.get("stages", [{}])):
        # adapt parameters for this stage
        run_params = gvec.util.CaseInsensitiveDict(copy.deepcopy(parameters))
        for key in ["stages", "Itor"]:
            if key in run_params:
                del run_params[key]
        for key, value in stage.items():
            if key in ["runs"]:
                continue
            if key in ["iota", "pres", "sgrid"]:
                if key not in run_params:
                    run_params[key] = {}
                for subkey, subvalue in value.items():
                    run_params[key][subkey] = subvalue
            if key in run_params and isinstance(value, Mapping):
                for subkey, subvalue in value.items():
                    run_params[key][subkey] = subvalue
            else:
                run_params[key] = value

        # run the stage
        runs = range(stage.get("runs", 1))
        for r in runs:
            progressstr = (
                "".join(
                    "|" + "=" * st.get("runs", 1) for st in parameters["stages"][:s]
                )
                + "|"
                + "=" * r
                + ">"
                + "." * (stage.get("runs", 1) - r - 1)
                + "|"
                + "".join(
                    "." * st.get("runs", 1) + "|"
                    for st in parameters["stages"][s + 1 :]
                )
            )
            if progressbar:
                print(f"GVEC stage {s} run {r}: {progressstr}", end="\r")
            logger.info(f"GVEC stage {s} run {r}: {progressstr}")
            start_time = time.time()
            # find previous state
            if statefile:
                logger.debug(f"Restart from statefile {statefile}")
                run_params["init_LA"] = False

            # prepare the run directory
            rundir = Path(f"{s:1d}-{r:02d}")
            if rundir.exists():
                logger.debug(f"Removing existing run directory {rundir}")
                shutil.rmtree(rundir)
            rundir.mkdir()
            logger.debug(f"Created run directory {rundir}")

            # write parameterfile & run GVEC
            gvec.util.write_parameter_file(
                gvec.util.flatten_parameters(run_params),
                rundir / "parameter.ini",
                header=f"!Auto-generated with `pygvec run` (stage {s} run {r})\n"
                "!Created at {datetime.now().isoformat()}\n"
                "!pyGVEC v{gvec.__version__}\n",
            )
            with gvec.util.chdir(rundir):
                gvec.run(
                    "parameter.ini",
                    ".." / statefile if statefile else None,
                    stdout_path="stdout.txt" if redirect_gvec_stdout else None,
                )

            # postprocessing
            statefile = sorted(rundir.glob("*State*.dat"))[-1]
            iterations = int(re.match(r".*State.*_(\d+)\.dat", statefile.name).group(1))
            max_iterations = run_params.get("maxiter")
            tolerance = run_params.get("minimize_tol")
            logger.debug(f"Postprocessing statefile {statefile}")

            with gvec.State(
                rundir / "parameter.ini",
                statefile,
                redirect_stdout=redirect_gvec_stdout,
            ) as state:
                ev = gvec.Evaluations(rho=rho, theta="int", zeta="int", state=state)
                state.compute(ev, "W_MHD", "N_FP")
                if "Itor" in parameters:
                    state.compute(ev, "iota", "iota_curr_0", "iota_0", "I_tor")

            if "Itor" in parameters:
                iota_values = ev.iota_0 + I_tor_target * ev.iota_curr_0
                run_params["iota"] = {
                    "type": "interpolation",
                    "vals": iota_values.data,
                    "rho2": (ev.rho**2).data,
                }

            # diagnostics
            # ToDo: possible early stop condition

            logger.info(f"W_MHD: {ev.W_MHD.item():.2e}")
            if "Itor" in parameters:
                iota_delta = ev.iota - iota_values
                logger.info(f"max Δiota: {np.abs(iota_delta).max().item():.2e}")
                logger.info(
                    f"rms Δiota: {np.sqrt((iota_delta**2).mean('rad')).item():.2e}"
                )
                logger.info(
                    f"max ΔItor: {np.abs(ev.I_tor - I_tor_target).max().item():.2e}"
                )

            d = xr.Dataset(
                dict(
                    W_MHD=ev.W_MHD,
                    gvec_iterations=iterations,
                    gvec_max_iterations=max_iterations,
                    gvec_tolerance=tolerance,
                )
            )
            if "Itor" in parameters:
                d["iota"] = ev.iota
                d["I_tor"] = ev.I_tor
                d["iota_delta"] = iota_delta
                d["I_tor_delta"] = ev.I_tor - I_tor_target
            d = d.drop_vars(["pol_weight", "tor_weight"])
            if diagnostics is None:
                d = d.expand_dims(dict(run=[r]))
                diagnostics = d
            else:
                d = d.expand_dims(dict(run=[diagnostics.run.size]))
                diagnostics = xr.concat([diagnostics, d], dim="run")
            if diagnosticfile:
                diagnostics.to_netcdf(diagnosticfile)

            end_time = time.time()
            logger.info(
                f"GVEC run took {end_time - start_time:5.1f} seconds for {iterations} iterations. (max {max_iterations}, tol {tolerance:.1e})"
            )
            logger.info("-" * 40)

    if plots:
        import matplotlib.pyplot as plt

        logger.debug("Plotting diagnostics...")

        if "Itor" in parameters:
            fig, axs = plt.subplots(1, 2, figsize=(10, 3), tight_layout=True)
        else:
            fig, ax = plt.subplots(1, 1, figsize=(5, 3), tight_layout=True)
            axs = [ax]
        axs[0].plot(diagnostics.run, diagnostics.W_MHD, ".-")
        axs[0].set(
            xlabel="run number",
            ylabel=f"${diagnostics.W_MHD.attrs['symbol']}$",
            title=diagnostics.W_MHD.attrs["long_name"],
        )
        if "Itor" in parameters:
            axs[1].plot(
                diagnostics.run, np.sqrt((diagnostics.iota_delta**2).mean("rad")), ".-"
            )
            axs[1].set(
                xlabel="run number",
                ylabel=r"$\sqrt{\sum \left(\Delta\iota\right)^2}$",
                title=f"Difference to target {diagnostics.iota.attrs['long_name']}\nroot mean square",
                yscale="log",
            )
        fig.savefig("iterations.png")

        if "Itor" in parameters:
            fig, axs = plt.subplots(
                2, 2, figsize=(15, 5), tight_layout=True, sharex=True
            )
            for r in diagnostics.run.data:
                if r == diagnostics.run.data[-1]:
                    kwargs = dict(marker=".", color="C0", alpha=1.0)
                else:
                    kwargs = dict(
                        color="black", alpha=0.2 + 0.3 * (r / diagnostics.run.data[-1])
                    )
                d = diagnostics.sel(run=r)
                axs[0, 0].plot(d.rho**2, d.iota, **kwargs)
                axs[1, 0].plot(d.rho**2, np.abs(d.iota_delta), **kwargs)
                axs[0, 1].plot(d.rho**2, d.I_tor, **kwargs)
                axs[1, 1].plot(d.rho**2, np.abs(d.I_tor_delta), **kwargs)
            for i, var in enumerate(["iota", "I_tor"]):
                axs[0, i].set(
                    title=diagnostics[var].attrs["long_name"],
                    ylabel=f"${diagnostics[var].attrs['symbol']}$",
                )
                axs[1, i].set(
                    title=f"Difference to target {diagnostics[var].attrs['long_name']}",
                    xlabel=r"$\rho^2$",
                    ylabel=f"$|\Delta {diagnostics[var].attrs['symbol']}|$",
                    yscale="log",
                )
            fig.savefig("profiles.png")

    logger.info("Done.")
    return rundir, statefile, diagnostics


def main(args=None):
    if args is None:
        args = parser.parse_args()
    if args.param_type is None:
        args.param_type = args.parameterfile.suffix[1:]

    if args.param_type == "ini":
        gvec.run(
            args.parameterfile,
            args.restartfile,
            stdout_path="stdout.txt" if args.quiet else None,
        )
    elif args.param_type in ["yaml", "toml"]:
        parameters = gvec.util.read_parameters(args.parameterfile)
        if "stages" not in parameters:
            parameters = gvec.util.flatten_parameters(parameters)
            parameterfile = f"{args.parameterfile.name}.ini"
            gvec.util.write_parameter_file(
                parameters,
                parameterfile,
                header=f"!Auto-generated from {args.parameterfile.name} with `pygvec run`\n!Created at {datetime.now().isoformat()}\n!pyGVEC v{gvec.__version__}\n",
            )
            gvec.run(
                parameterfile,
                args.restartfile,
                stdout_path="stdout.txt" if args.quiet else None,
            )
        else:
            logging.basicConfig(
                level=logging.WARNING
            )  # show warnings and above as normal
            logger = logging.getLogger(
                "pyGVEC.script"
            )  # show info/debug messages for this script
            logger.propagate = False
            loghandler = logging.StreamHandler()
            logformatter = logging.Formatter("{message}", style="{")
            loghandler.setFormatter(logformatter)
            logger.addHandler(loghandler)
            if args.verbose == 1:
                logger.setLevel(logging.INFO)
            elif args.verbose >= 2:
                logger.setLevel(logging.DEBUG)
            run_stages(
                parameters,
                args.restartfile,
                progressbar=not args.quiet and not args.verbose,
                redirect_gvec_stdout=args.verbose < 3,
                diagnosticfile=args.diagnostics,
                plots=args.plots,
            )
    else:
        raise ValueError("Cannot determine parameterfile type")


if __name__ == "__main__":
    main()
