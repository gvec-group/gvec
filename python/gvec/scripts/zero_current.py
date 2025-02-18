"""Run GVEC multiple times with Picard iterations on the iota profile to reach zero current.

Run GVEC multiple times with different parameters and restart from the final state of each previous run.
The new iota profile that excludes the current contribution is computed from the previous state, fitted to a polynomial and prescribed via the parameterfile to the next run.
Along the Picard iteration, radial resolution and maximum iteration count is increased, as well.

Developer's Notes:
* Providing the script as a function allows for easier integration and testing.
* The main() function is used as an entry point in pyproject.toml.
* The "__name__" check allows the script to be run standalone or using "python -m gvec.scripts.zero_current".
"""

# === Imports === #

import time
from pathlib import Path
import shutil
import logging
import argparse

import numpy as np
import xarray as xr
from tqdm import tqdm

import gvec

# === Arguments === #

parser = argparse.ArgumentParser(
    prog="gvec_zero_current",
    description="Run GVEC multiple times with Picard iterations on the iota profile to reach zero current.",
)
parser.add_argument(
    "parameterfile", type=Path, help="input (template) GVEC parameter-file"
)
parser.add_argument(
    "--picard-iterations",
    "--pi",
    type=int,
    default=5,
    help="number of Picard iterations (GVEC restarts)",
)
parser.add_argument(
    "--gvec-iterations",
    "--gi",
    type=int,
    nargs=2,
    default=(100, 1000),
    help="minimum and maximum number of GVEC iterations (min is doubled until max)",
)
parser.add_argument(
    "--gvec-nelems",
    "--ge",
    type=int,
    nargs=2,
    default=(2, 10),
    help="minimum and maximum number of elements in the sgrid (increased linearly)",
)
parser.add_argument(
    "--iota-poly-degree",
    "--ipd",
    type=int,
    default=9,
    help="degree of the polynomial fit for the iota profile",
)
parser.add_argument(
    "--reinit-LA", action="store_true", help="initialize LA at every restart"
)
parser.add_argument(
    "-o",
    "--outputfile",
    type=Path,
    default="GVEC_zero_current.nc",
    help="output netCDF file for diagnostics",
)
parser.add_argument("-p", "--plots", action="store_true", help="plot diagnostics")
parser.add_argument(
    "-v",
    "--verbose",
    action="count",
    default=0,
    help="verbosity level: -v for info, -vv for debug, -vvv for GVEC output",
)
parser.add_argument("-q", "--quiet", action="store_true", help="suppress output")

# === Script === #


def run_zero_current(
    template: str | Path,
    max_iteration: int,
    gvec_miniter: int,
    gvec_maxiter: int,
    gvec_nelems_min: int,
    gvec_nelems_max: int,
    iota_poly_degree: int,
    reinit_LA: bool,
    outputfile: str | Path,
    progressbar: bool = True,
    gvec_stdout_path: str | Path | None = "stdout.txt",
    plots: bool = False,
):
    """Run GVEC multiple times with Picard iterations on the iota profile to reach zero current.

    Parameters:
    -----------
    template : str | Path
        The template parameter file.
    max_iteration : int
        The maximum number of iterations to run.
    gvec_miniter : int
        The minimum number of iterations for GVEC.
    gvec_maxiter : int
        The maximum number of iterations for GVEC.
    gvec_nelems_min : int
        The minimum number of elements in the sgrid.
    gvec_nelems_max : int
        The maximum number of elements in the sgrid.
    iota_poly_degree : int
        The degree of the polynomial fit for the iota profile.
    reinit_LA : bool
        Whether to initialize LA at restart.
    progressbar : bool
        Whether to show a progress bar.
    gvec_stdout_path : str | Path | None
        The path to the GVEC stdout file or None to write to stdout.
    """
    rho = np.sqrt(np.linspace(0, 1, 41)[1:])

    ev: xr.Dataset = None  # Postprocssing of latest iteration
    statefile: Path = None  # Last state of latest iteration
    diagnostics: xr.Dataset = None
    logger = logging.getLogger("pyGVEC.script")

    if plots:
        import matplotlib.pyplot as plt

    iterations = range(max_iteration + 1)
    if progressbar:
        iterations = tqdm(iterations, ascii=True, desc="Picard iteration", ncols=80)
    for i in iterations:
        logger.info(
            f"Iteration {i:2d}/{max_iteration} "
            + "=" * i
            + ">"
            + "." * (max_iteration - i)
        )
        start_time = time.time()
        path = Path(f"{i:02d}")
        params = {}
        # find previous state
        if i > 0:
            # from postprocessing of previous iteration: ev, statefile
            logger.debug(f"Restart from statefile {statefile}")

            # get polynomial fit of iota - iota_curr
            iota_coefs = np.polyfit(rho**2, ev.iota_0, iota_poly_degree)
            logger.debug(f"Setting iota_coefs: {iota_coefs[::-1]}")

            # set new parameters
            params["init_LA"] = "T" if reinit_LA else "F"
            params["sign_iota"] = 1
            # currently adapt_parameter_file expects strings for iota coefficients
            params["iota_coefs"] = "(/" + ", ".join(map(str, iota_coefs[::-1])) + "/)"

        params["maxiter"] = int(np.amin([gvec_maxiter, gvec_miniter * 2**i]))
        params["sgrid_nElems"] = int(
            gvec_nelems_min + (gvec_nelems_max - gvec_nelems_min) * (i / max_iteration)
        )
        logger.debug(f"Setting maxiter to {params['maxiter']}")
        logger.debug(f"Setting sgrid_nElems to {params['sgrid_nElems']}")

        # prepare the run directory
        if path.exists():
            logger.debug(f"Removing existing run directory {path}")
            shutil.rmtree(path)
        if not path.exists():
            path.mkdir()
            logger.debug(f"Created run directory {path}")

        # copy & modify parameterfile
        gvec.util.adapt_parameter_file(template, path / "parameter.ini", **params)

        # run gvec
        with gvec.util.chdir(path):
            # ToDo: GVEC should raise an error if it fails
            gvec.run(
                "parameter.ini",
                ".." / statefile if statefile else None,
                stdout_path=gvec_stdout_path,
            )

        # postprocessing
        statefile = sorted(path.glob("*State*.dat"))[-1]
        logger.debug(f"Using statefile {statefile}")

        with gvec.State(path / "parameter.ini", statefile) as state:
            ev = gvec.Evaluations(rho=rho, theta="int", zeta="int", state=state)
            state.compute(ev, "iota", "iota_curr", "iota_0", "I_tor", "N_FP")
            # ToDo: merge into one Evaluations once MR!55 is merged
            ev_vol = gvec.Evaluations(rho="int", theta="int", zeta="int", state=state)
            state.compute(ev_vol, "W_MHD", "iota_curr")

        # get polynomial fit of iota_0 = iota - iota_curr
        iota_coefs = np.polyfit(rho**2, ev.iota_0, iota_poly_degree)

        # diagnostics
        # ToDo: possible early stop condition
        # ToDo: iota_curr seems to diverge near the axis -> check boundary conditions
        # iota_curr_rms = np.sqrt(gvec.comp.radial_integral(ev_vol.iota_curr**2))
        iota_curr_rms = np.sqrt((ev.iota_curr**2).mean("rad"))

        logger.info(f"W_MHD: {ev_vol.W_MHD.item():.3e}")
        logger.info(f"max iota_curr: {np.abs(ev.iota_curr).max().item():.3f}")
        logger.info(f"rms iota_curr: {iota_curr_rms.item():.3f}")
        logger.info(f"max I_tor: {np.abs(ev.I_tor).max().item():.3e}")

        d = xr.Dataset(
            dict(
                W_MHD=ev_vol.W_MHD,
                iota=ev.iota,
                iota_curr=ev.iota_curr,
                iota_0=ev.iota_0,
                I_tor=ev.I_tor,
                N_FP=ev.N_FP,
                gvec_iterations=params["maxiter"],
                gvec_sgrid_nElems=params["sgrid_nElems"],
            )
        ).expand_dims(dict(iteration=[i]))
        d = d.drop_vars(["pol_weight", "tor_weight"])
        if i == 0:
            diagnostics = d
        else:
            diagnostics = xr.concat([diagnostics, d], dim="iteration")
        diagnostics.to_netcdf(outputfile)

        end_time = time.time()
        logger.info(f"Iteration took {end_time-start_time:5.1f} seconds.")
        logger.info("-" * 40)

    if plots:
        logger.debug("Plotting diagnostics...")

        fig, axs = plt.subplots(1, 2, figsize=(10, 3), tight_layout=True)
        axs[0].plot(
            diagnostics.iteration, np.sqrt((diagnostics.iota_curr**2).mean("rad")), ".-"
        )
        axs[0].set(
            xlabel="picard iterations",
            ylabel=r"$\sqrt{\sum \iota_{\mathrm{curr}}^2 }$",
            title=f"{diagnostics.iota_curr.attrs['long_name']}\nroot mean square",
            yscale="log",
        )
        axs[1].plot(diagnostics.iteration, diagnostics.W_MHD, ".-")
        axs[1].set(
            xlabel="picard iterations for curr. constraint",
            ylabel=f"${diagnostics.W_MHD.attrs['symbol']}$",
            title=diagnostics.W_MHD.attrs["long_name"],
        )
        fig.savefig("iterations.png")

        fig, axs = plt.subplots(1, 3, figsize=(15, 5), tight_layout=True, sharex=True)
        for i in diagnostics.iteration.data:
            if i == max_iteration:
                kwargs = dict(marker=".", color="C0", alpha=1.0)
            else:
                kwargs = dict(color="black", alpha=0.2 + 0.3 * (i / max_iteration))
            d = diagnostics.sel(iteration=i)
            axs[0].plot(d.rho**2, d.iota, **kwargs)
            axs[1].plot(d.rho**2, np.abs(d.iota_curr), **kwargs)
            axs[2].plot(d.rho**2, np.abs(d.I_tor), **kwargs)
        for i, var in enumerate(["iota", "iota_curr", "I_tor"]):
            axs[i].set(
                title=diagnostics[var].attrs["long_name"],
                xlabel=r"$\rho^2$",
                ylabel=f"$|{diagnostics[var].attrs['symbol']}|$",
            )
        axs[1].set_yscale("log")
        axs[2].set_yscale("log")
        fig.savefig("profiles.png")

    logger.info("Done.")
    return 0


def main() -> int:
    args = parser.parse_args()
    if args.quiet and args.verbose:
        raise ValueError("Cannot be quiet and verbose at the same time.")

    logging.basicConfig(level=logging.WARNING)  # show warnings and above as normal
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

    run_zero_current(
        template=args.parameterfile,
        max_iteration=args.picard_iterations - 1,
        gvec_miniter=args.gvec_iterations[0],
        gvec_maxiter=args.gvec_iterations[1],
        gvec_nelems_min=args.gvec_nelems[0],
        gvec_nelems_max=args.gvec_nelems[1],
        iota_poly_degree=args.iota_poly_degree,
        reinit_LA=args.reinit_LA,
        outputfile=args.outputfile,
        progressbar=not args.quiet and not args.verbose,
        gvec_stdout_path=None if args.verbose >= 3 else "stdout.txt",
        plots=args.plots,
    )
    return 0


if __name__ == "__main__":
    main()
