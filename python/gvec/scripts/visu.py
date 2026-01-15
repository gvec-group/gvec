import argparse
import logging
from collections.abc import Sequence
from pathlib import Path

from gvec import State, find_state
from gvec.util import logging_setup

parser = argparse.ArgumentParser(
    prog="pygvec-plots", description="Output some default plots for diagnostics."
)

parser.add_argument("--rundir", type=Path, help="GVEC run directory", default=Path("."))

parser.add_argument("--outdir", type=Path, help="Plot output directory", default=Path("."))

parser.add_argument(
    "-v",
    "--verbose",
    action="count",
    default=0,
    help="verbosity level: -v for info, -vv for debug",
)


def write_plots_to_directory(state: State, output_directory: Path, verbose: bool = False):
    pass


def main(args: Sequence[str] | argparse.Namespace | None = None):
    if isinstance(args, argparse.Namespace):
        pass
    else:
        args = parser.parse_args(args)

    logging_setup()
    logger = logging.getLogger(__name__)
    if args.verbose == 0:
        logging.disable()
    if args.verbose >= 2:
        logger.setLevel(logging.DEBUG)
    elif args.verbose == 1:
        logger.setLevel(logging.INFO)
    logger.debug(f"parsed args: {args}")

    state = find_state(args.rundir)

    logger.info("Generating:")

    logger.info("   radial profiles.")
    f, ax = state.plot_radial_profile()
    f.savefig(args.outdir / "profiles.png")

    logger.info("   axis plots.")
    f, ax = state.plot_on_axis()
    f.savefig(args.outdir / "modB_axis.png")

    logger.info("   poloidal cuts.")
    f, ax = state.plot_poloidal_plane()
    f.savefig(args.outdir / "modB_poloidal_cuts.png")

    logger.info("   last closed flux surface.")
    f, ax = state.plot_on_flux_surface()
    f.savefig(args.outdir / "modB_lcfs.png")

    logger.info("Done.")

    pass


if __name__ == "__main__":
    main()
