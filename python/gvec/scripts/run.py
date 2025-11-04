# Copyright (c) 2025 GVEC Contributors, Max Planck Institute for Plasma Physics
# License: MIT
"""The pyGVEC run script for running GVEC using stages and current optimization."""

import argparse
import logging

from pathlib import Path
from collections.abc import Sequence

import numpy as np
import xarray as xr
import gvec

# === Argument Parser === #

parser = argparse.ArgumentParser(
    prog="pygvec-run",
    formatter_class=argparse.RawDescriptionHelpFormatter,
    description="Run GVEC with a given parameterfile, optionally restarting from an existing statefile.\n\n"
    "When given an INI parameterfile, GVEC is called directly.\n"
    "With YAML and TOML parameterfiles, GVEC can be run in several stages and a current optimization with picard iterations can be performed.",
)
parser.add_argument("parameterfile", type=Path, help="input GVEC parameterfile")
parser.add_argument(
    "restart_statefile",
    type=Path,
    help="GVEC statefile to restart from (optional)",
    nargs="?",
)
parser.add_argument(
    "restart_parameterfile",
    type=Path,
    help="GVEC parameterfile to restart from (optional)",
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
parser.add_argument(
    "-p",
    "--plots",
    action="count",
    help="plot diagnostics (-pp for additional plots)",
)

parser.add_argument(
    "-k",
    "--keep",
    action="count",
    default=0,
    help="keep intermediate results: -k for the last restarts of each stage , -kk for all intermediate results",
)

# === Script === #


@gvec.errors.without_traceback
def main(args: Sequence[str] | argparse.Namespace | None = None):
    if not isinstance(args, argparse.Namespace):
        args = parser.parse_args(args)
    if args.param_type is None:
        args.param_type = args.parameterfile.suffix[1:]

    if args.param_type not in ["ini", "yaml", "toml"]:
        raise ValueError("Cannot determine parameterfile type")

    parameters = gvec.util.read_parameters(args.parameterfile)

    if args.quiet:
        logging.disable()
    if args.verbose == 0:
        loglevel = logging.WARNING
    elif args.verbose == 1:
        loglevel = logging.INFO
    elif args.verbose >= 2:
        loglevel = logging.DEBUG

    if args.keep == 0:
        keep_intermediates = None
    elif args.keep == 1:
        keep_intermediates = "stages"
    elif args.keep >= 2:
        keep_intermediates = "all"

    if args.restart_parameterfile:
        restart = gvec.load_state(args.restart_parameterfile, args.restart_statefile)
    elif args.restart_statefile:
        restart = args.restart_statefile
    else:
        restart = None

    run_with_stages = gvec.run(
        parameters,
        restart,
        quiet=args.quiet,
        redirect_gvec_stdout=args.verbose < 3,
        keep_intermediates=keep_intermediates,
        loglevel=loglevel,
    )

    if args.diagnostics:
        diagnostics = xr.merge(
            [run_with_stages.diagnostics_run, run_with_stages.diagnostics_minimizer]
        )
        diagnostics.to_netcdf(args.diagnostics)
    if args.plots:
        try:
            if np.sum(np.array(run_with_stages.n_runs_in_stage) > 0) >= 2 and args.plots >= 2:
                fig_runs = run_with_stages.plot_diagnostics_run()
                fig_runs.savefig(f"{run_with_stages._state_parameters['projectName']}_runs.png")

            if run_with_stages.curr_constraint and args.plots >= 2:
                fig_profiles = run_with_stages.plot_diagnostics_current_profiles()
                fig_profiles.savefig(
                    f"{run_with_stages._state_parameters['projectName']}_profiles.png"
                )
            fig_minimization = run_with_stages.plot_diagnostics_minimization()
            fig_minimization.savefig(
                f"{run_with_stages._state_parameters['projectName']}_iterations.png"
            )
        except ImportError as e:
            logging.debug(f"Plotting failed: {e}")
            logging.error("Plotting requires matplotlib, which is not installed.")


if __name__ == "__main__":
    main()
