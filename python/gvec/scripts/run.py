# Copyright (c) 2025 GVEC Contributors, Max Planck Institute for Plasma Physics
# License: MIT
"""The pyGVEC run script for running GVEC using stages and current constraints."""

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

parser.add_argument(
    "-k",
    "--keep",
    action="count",
    default=0,
    help="keep intermediate results: -k for the last restarts of each stage , -kk for all intermediate results",
)

# === Script === #


def main(args: Sequence[str] | argparse.Namespace | None = None):
    if isinstance(args, argparse.Namespace):
        pass
    else:
        args = parser.parse_args(args)
    if args.param_type is None:
        args.param_type = args.parameterfile.suffix[1:]

    if args.param_type in ["ini", "yaml", "toml"]:
        parameters = gvec.util.read_parameters(args.parameterfile)

        logging.basicConfig(level=logging.WARNING)  # show warnings and above as normal
        logger = logging.getLogger("pyGVEC.script")  # show info/debug messages for this script
        logger.propagate = False
        loghandler = logging.StreamHandler()
        logformatter = logging.Formatter("{levelname} {message}", style="{")
        loghandler.setFormatter(logformatter)
        logger.addHandler(loghandler)
        if args.verbose == 1:
            logger.setLevel(logging.INFO)
        elif args.verbose >= 2:
            logger.setLevel(logging.DEBUG)

        if args.keep == 0:
            delete_intermediates = "all"
        elif args.keep == 1:
            delete_intermediates = "restarts"
        elif args.keep >= 2:
            delete_intermediates = None

        run_with_stages = gvec.run(
            parameters,
            args.restartfile,
            progressbar=not args.quiet and not args.verbose,
            redirect_gvec_stdout=args.verbose < 3,
            delete_intermediates=delete_intermediates,
        )

        if args.diagnostics:
            diagnostics = xr.merge(
                [run_with_stages.diagnostics_run, run_with_stages.diagnostics_minimizer]
            )
            diagnostics.to_netcdf(args.diagnostics)
        if args.plots:
            if np.sum(np.array(run_with_stages.n_runs_in_stage) > 0) >= 2:
                fig_runs = run_with_stages.plot_diagnostics_run()
                fig_runs.savefig(f"{run_with_stages.state_parameters['projectName']}_runs.png")

            if run_with_stages.curr_constraint:
                fig_profiles = run_with_stages.plot_diagnostics_current_profiles()
                fig_profiles.savefig(
                    f"{run_with_stages.state_parameters['projectName']}_profiles.png"
                )
            fig_minimization = run_with_stages.plot_diagnostics_minimization()
            fig_minimization.savefig(
                f"{run_with_stages.state_parameters['projectName']}_iterations.png"
            )
    else:
        raise ValueError("Cannot determine parameterfile type")


if __name__ == "__main__":
    main()
