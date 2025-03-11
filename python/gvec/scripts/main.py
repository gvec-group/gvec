# Copyright (c) 2025 GVEC Contributors, Max Planck Institute for Plasma Physics
# License: MIT
"""The pyGVEC executable"""

# === Imports === #

import platform
import time
from pathlib import Path
import shutil
import logging
import argparse
import re
from typing import Mapping
from datetime import datetime
import copy

import numpy as np
import xarray as xr
import yaml
import tomlkit  # also supports writing (compared to tomllib)

import gvec
from gvec.scripts import to_cas3d

# === Arguments === #

parser = argparse.ArgumentParser(
    prog="pygvec",
    description="GVEC: a 3D MHD equilibrium solver",
)
subparsers = parser.add_subparsers(
    title="mode",
    description="which mode/subcommand to run",
    dest="mode",
)
parser.add_argument(
    "-V",
    "--version",
    action="version",
    version=f"pyGVEC v{gvec.__version__} from {Path(gvec.__file__).parent} (python {platform.python_version()})",
)

# --- run --- #

run_parser = subparsers.add_parser(
    "run",
    help="run GVEC",
    formatter_class=argparse.RawDescriptionHelpFormatter,
    description="Run GVEC with a given parameterfile, optionally restarting from an existing statefile.\n\n"
    "When given an INI parameterfile, GVEC is called directly.\n"
    "With YAML and TOML parameterfiles, GVEC can be run in several stages and a current constraint with picard iterations can be performed.",
)
run_parser.add_argument("parameterfile", type=Path, help="input GVEC parameterfile")
run_parser.add_argument(
    "restartfile",
    type=Path,
    help="GVEC statefile to restart from (optional)",
    nargs="?",
)
run_parser_param_type = run_parser.add_mutually_exclusive_group()
run_parser_param_type.add_argument(
    "--ini",
    action="store_const",
    const="ini",
    dest="param_type",
    help="interpret GVEC parameterfile classicly (INI)",
)
run_parser_param_type.add_argument(
    "--yaml",
    action="store_const",
    const="yaml",
    dest="param_type",
    help="interpret GVEC parameterfile as YAML",
)
run_parser_param_type.add_argument(
    "--toml",
    action="store_const",
    const="toml",
    dest="param_type",
    help="interpret GVEC parameterfile as TOML",
)
run_parser_verbosity = run_parser.add_mutually_exclusive_group()
run_parser_verbosity.add_argument(
    "-v",
    "--verbose",
    action="count",
    default=0,
    help="verbosity level: -v for info, -vv for debug, -vvv for GVEC output",
)
run_parser_verbosity.add_argument(
    "-q", "--quiet", action="store_true", help="suppress output"
)
run_parser.add_argument(
    "-d",
    "--diagnostics",
    type=Path,
    default="GVEC-diagnostics.nc",
    help="output netCDF file for diagnostics",
)
run_parser.add_argument("-p", "--plots", action="store_true", help="plot diagnostics")

# --- convert parameterfile --- #

convert_parser = subparsers.add_parser(
    "convert-params",
    help="convert the GVEC parameterfile between different formats",
    formatter_class=argparse.RawDescriptionHelpFormatter,
    description="Convert GVEC parameterfiles between different formats.\n"
    "The INI (classical) parameter files do not support stages or the current constraint!\nAlso the formatting is lost upon conversion.",
)
convert_parser.add_argument(
    "input",
    type=Path,
    help="input GVEC parameterfile",
)
convert_parser.add_argument(
    "output",
    type=Path,
    help="output GVEC parameterfile",
)

# --- other scripts --- #

cas3d_parser = subparsers.add_parser(
    "to-cas3d",
    help="convert a GVEC state to a CAS3D compatible input file",
    formatter_class=argparse.RawDescriptionHelpFormatter,
    description="Convert a GVEC statefile to a CAS3D compatible input file.",
    parents=[to_cas3d.parser],
    add_help=False,
)

# === Script === #


def main():
    logging.basicConfig(level=logging.WARNING)  # show warnings and above as normal
    logger = logging.getLogger("pyGVEC.script")
    args = parser.parse_args()

    # --- run GVEC --- #
    if args.mode == "run":
        if args.param_type is None:
            args.param_type = args.parameterfile.suffix[1:]

        if args.param_type == "ini":
            gvec.run(
                args.parameterfile,
                args.restartfile,
                stdout_path="stdout.txt" if args.quiet else None,
            )
        elif args.param_type in ["yaml", "toml"]:
            with open(args.parameterfile, "r") as file:
                if args.param_type == "yaml":
                    parameters = yaml.safe_load(file)
                elif args.param_type == "toml":
                    parameters = tomlkit.parse(file.read()).unwrap()
            parameters = unstringify_mn_params(parameters)
            if "stages" not in parameters:
                parameters = flatten_params(parameters)
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

    # --- convert parameterfile --- #
    elif args.mode == "convert-params":
        if args.input.suffix == ".ini":
            inputs = gvec.util.read_parameter_file(args.input)
            inputs = stack_params(inputs)
        elif args.input.suffix == ".yaml":
            with open(args.input, "r") as file:
                inputs = yaml.safe_load(file)
            inputs = unstringify_mn_params(inputs)
        elif args.input.suffix == ".toml":
            with open(args.input, "r") as file:
                inputs = tomlkit.parse(file.read()).unwrap()
            inputs = unstringify_mn_params(inputs)
        else:
            raise ValueError(f"Unknown input file type: {args.input}")
        if args.output.suffix == ".ini":
            outputs = flatten_params(inputs)
            gvec.util.write_parameter_file(outputs, args.output)
        elif args.output.suffix == ".yaml":
            outputs = stringify_mn_params(inputs)
            with open(args.output, "w") as file:
                yaml.safe_dump(outputs, file)
        elif args.output.suffix == ".toml":
            outputs = stringify_mn_params(inputs)
            with open(args.output, "w") as file:
                file.write(
                    tomlkit.dumps(outputs)
                )  # ToDo: nicer output using document API
        else:
            raise ValueError(f"Unknown output file type: {args.output}")

    # --- other scripts --- #
    elif args.mode == "to-cas3d":
        to_cas3d.main(args)


def run_stages(
    parameters: Mapping,
    statefile: Path | None = None,
    progressbar: bool = False,
    redirect_gvec_stdout: bool = True,
    diagnosticfile: Path | None = None,
    plots: bool = False,
):
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
                Itor = np.poly1d(coefs)
            case _:
                raise ValueError(f"Unknown Itor type: {parameters['Itor']['type']}")

    for s, stage in enumerate(parameters["stages"]):
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
                flatten_params(run_params),
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
                Itor_values = Itor(rho**2)
                iota_values = ev.iota_0 + Itor_values * ev.iota_curr_0
                run_params["iota"] = {
                    "type": "interpolation",
                    "vals": iota_values.data,
                    "rho2": (ev.rho**2).data,
                }

            # diagnostics
            # ToDo: possible early stop condition

            logger.info(f"W_MHD: {ev.W_MHD.item():.3e}")
            if "Itor" in parameters:
                iota_delta = ev.iota - iota_values
                logger.info(f"max Δiota: {np.abs(iota_delta).max().item():.3f}")
                logger.info(
                    f"rms Δiota: {np.sqrt((iota_delta**2).mean('rad')).item():.3f}"
                )
                logger.info(
                    f"max ΔItor: {np.abs(ev.I_tor - Itor_values).max().item():.3e}"
                )

            d = xr.Dataset(
                dict(
                    W_MHD=ev.W_MHD,
                    gvec_iterations=run_params["maxiter"],
                )
            )
            if "Itor" in parameters:
                d["iota"] = ev.iota
                d["I_tor"] = ev.I_tor
                d["iota_delta"] = iota_delta
                d["I_tor_delta"] = ev.I_tor - Itor_values
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
            logger.info(f"GVEC run took {end_time-start_time:5.1f} seconds.")
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


def stack_params(parameters):
    """Stack parameters into a hierarchical dictionary"""
    output = {}
    for key, value in parameters.items():
        if "_" not in key:
            output[key] = value
            continue
        group, name = key.split("_", 1)
        if group in ["iota", "pres", "sgrid"]:
            if group not in output:
                output[group] = {}
            output[group][name] = value
        else:
            output[key] = value
    return output


def flatten_params(parameters):
    """Flatten parameters from a hierarchical dictionary"""
    output = {}
    for key, value in parameters.items():
        if isinstance(value, dict) and not re.match(r"(X1|X2|LA)_[a|b]_(sin|cos)", key):
            if key in ["stages", "Itor"]:  # not supported by fortran-GVEC
                continue
            for subkey, subvalue in value.items():
                output[f"{key}_{subkey}"] = subvalue
        else:
            output[key] = value
    return output


def stringify_mn_params(parameters):
    """Serialize parameters into a string"""
    output = {}
    for key, value in parameters.items():
        if re.match(r"(X1|X2|LA)_[a|b]_(sin|cos)", key):
            output[key] = {}
            for (m, n), val in value.items():
                output[key][f"({m}, {n:2d})"] = val
        else:
            output[key] = value
    return output


def unstringify_mn_params(parameters):
    """Deserialize parameters from a string"""
    output = {}
    for key, value in parameters.items():
        if re.match(r"(X1|X2|LA)_[a|b]_(sin|cos)", key):
            output[key] = {}
            for mn, val in value.items():
                m, n = map(int, mn.strip("()").split(","))
                output[key][(m, n)] = val
        else:
            output[key] = value
    return output


if __name__ == "__main__":
    main()
