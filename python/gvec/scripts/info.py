# Copyright (c) 2025 GVEC Contributors, Max Planck Institute for Plasma Physics
# License: MIT
"""info.py - Print overview information for a given GVEC State."""

import argparse
import logging
from collections.abc import Sequence
from pathlib import Path

from tabulate import tabulate
import numpy as np

from gvec import State, find_state, volume_integral
from gvec.util import logging_setup

parser = argparse.ArgumentParser(
    prog="pygvec-info", description="Print overview information for a given GVEC State."
)

parser.add_argument("--rundir", type=Path, help="GVEC run directory", default=Path("."))

verbosity = parser.add_mutually_exclusive_group()
verbosity.add_argument(
    "-v",
    "--verbose",
    action="count",
    default=0,
    help="verbosity level: -v for info, -vv for debug",
)
verbosity.add_argument("-q", "--quiet", action="store_true", help="suppress log messages")

style = parser.add_mutually_exclusive_group()
style.add_argument(
    "--raw",
    action="store_true",
    help="format output as semicolon-separated values for easier parsing",
)
style.add_argument(
    "--style",
    default="plain",
    help="table format for tabulate (e.g. 'github' or 'fancy_outline')",
)


def get_info(ev, q: str, comment: str | None = None):
    return [q, ev[q].item(), ev[q].attrs["long_name"] + (f" ({comment})" if comment else "")]


def main(args: Sequence[str] | argparse.Namespace | None = None):
    if isinstance(args, argparse.Namespace):
        pass
    else:
        args = parser.parse_args(args)

    logging_setup()
    logger = logging.getLogger("gvec")
    if args.quiet:
        logging.disable()
    if args.verbose >= 2:
        logger.setLevel(logging.DEBUG)
    elif args.verbose == 1:
        logger.setLevel(logging.INFO)
    logger.debug(f"parsed args: {args}")

    state = find_state(args.rundir)

    logger.info("Evaluating State")

    Q_geo = ["V", "L_axis", "A_surface", "r_major", "r_minor", "aspect_ratio", "elongation"]
    Q_other = ["beta_avg", "mod_F_rel_volavg", "W_MHD", "vacuum_magnetic_well_depth"]
    ev_vol = state.evaluate(
        *Q_geo,
        "iota",
        "iota_avg",
        "shear_avg",
        "mod_B",
        "Phi_edge",
        "p",
        *Q_other,
    )
    ev_lcfs = state.evaluate(
        "I_tor",
        "I_pol",
        "iota",
        "mirror_ratio",
        "L_gradB",
        rho=1.0,
    )

    output = [get_info(ev_vol, q) for q in Q_geo]
    R_avg = (
        volume_integral(np.sqrt(ev_vol.pos[0] ** 2 + ev_vol.pos[1] ** 2) * ev_vol.Jac)
        / ev_vol.V
    )
    output += [["R_avg", R_avg.item(), "average distance from Z-axis"]]
    output += [
        get_info(ev_vol, "iota_avg"),
        ["iota_max", ev_vol.iota.max().item(), "maximum rotational transform"],
        ["iota_min", ev_vol.iota.min().item(), "minimum rotational transform"],
        ["iota_edge", ev_lcfs.iota.item(), "rotational transform at the edge (rho=1.0)"],
        get_info(ev_vol, "shear_avg"),
    ]
    output += [
        ["I_tor", ev_lcfs.I_tor.item(), "total toroidal current"],
        ["I_pol", ev_lcfs.I_pol.item(), "total poloidal current"],
    ]
    B_avg = volume_integral(ev_vol.mod_B * ev_vol.Jac) / ev_vol.V
    output += [
        get_info(ev_vol, "Phi_edge"),
        ["B_avg", B_avg.item(), "average magnetic field strength"],
        ["B_max", ev_vol.mod_B.max().item(), "maximum magnetic field strength"],
        ["B_min", ev_vol.mod_B.min().item(), "minimum magnetic field strength"],
        get_info(ev_lcfs, "mirror_ratio", "on LCFS"),
    ]
    output += [
        ["p_max", ev_vol.p.max().item(), "maximum pressure"],
    ]
    output += [get_info(ev_vol, q) for q in Q_other]
    output += [
        [
            "L_gradB_min",
            ev_lcfs.L_gradB.min().item(),
            "minimum mag. gradient length scale (on LCFS)",
        ]
    ]

    if args.raw:
        print("\n".join(";".join(map(str, row)) for row in output))
    else:
        print(tabulate(output, numalign="right", tablefmt=args.style))

    logger.info("Done.")


if __name__ == "__main__":
    main()
