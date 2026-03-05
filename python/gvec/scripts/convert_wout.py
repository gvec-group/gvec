# Copyright (c) 2026 GVEC Contributors, Max Planck Institute for Plasma Physics
# License: MIT
"""convert-wout.py - convert a VMEC wout file to a GVEC parameterfile & state"""

from pathlib import Path
import logging
import argparse
from collections.abc import Sequence, Mapping
import tempfile
import shutil
import re

import numpy as np
import xarray as xr

import gvec
from gvec.util import CaseInsensitiveDict

parser = argparse.ArgumentParser(
    prog="pygvec-convert-wout",
    formatter_class=argparse.RawDescriptionHelpFormatter,
    description="Convert a VMEC wout file to a GVEC parameterfile & state.",
)
parser.add_argument(
    "wout_file",
    type=Path,
    help="input VMEC wout file",
)
parser.add_argument(
    "rundir",
    type=Path,
    default=".",
    nargs="?",
    help="optional output directory for GVEC state and parameter files (default: current directory)",
)
wout_format = parser.add_argument_group(title="wout format").add_mutually_exclusive_group()
wout_format.add_argument(
    "--nemec-ascii",
    action="store_const",
    const=1,
    dest="VMECwoutfile_format",
    help="input wout file is in nemec-ascii format (default: netCDF)",
)
wout_format.add_argument(
    "--nemec-binary",
    action="store_const",
    const=2,
    dest="VMECwoutfile_format",
    help="input wout file is in nemec-binary format (default: netCDF)",
)
parser.add_argument(
    "--param",
    type=str,
    nargs="*",
    default=[],
    help="extra parameter to be used in the conversion, e.g. --param=init_LA=True or --param=sgrid_nElems=20",
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

logger = logging.getLogger(__name__)


def convert_vmec_wout(
    wout_file: Path,
    rundir: Path | None = None,
    *,
    VMECwoutfile_format: int = 0,
    extra_parameters: Mapping = {},
) -> CaseInsensitiveDict:
    """
    Convert a VMEC wout file to a GVEC parameter & statefile.

    The boundary and profiles are extracted, to generate a GVEC parameter file corresponding to the same configuration,
    as well as a GVEC parameter and statefile describing the (interpolated) VMEC state.

    Parameters
    ----------
    wout_file
        Path to the VMEC wout file in netCDF format.
    rundir
        Optional path to a directory in which the GVEC state should be generated.
    VMECwoutfile_format
        The format of the VMEC wout file (0 for netCDF, 1 for nemec-ascii, 2 for nemec-binary).
    extra_parameters
        Extra parameters to be used in the conversion, e.g. ``init_LA`` or ``sgrid_nElems``.
    """
    if not Path(wout_file).exists():
        raise FileNotFoundError(f"wout file '{wout_file}' not found")
    wout_ds = xr.open_dataset(wout_file)

    M = int(wout_ds.xm.max())
    N = int(wout_ds.xn.max() / wout_ds.nfp)
    nfp = int(wout_ds.nfp)
    lasym = bool(wout_ds.lasym__logical__)  # stellarator asymmetric
    phiedge = float(wout_ds.phi[-1].item())
    default_parameters = gvec.util.CaseInsensitiveDict(
        which_hmap=1,
        nfp=nfp,
        X1_mn_max=[M, N],
        X2_mn_max=[M, N],
        LA_mn_max=[M, N],
        X1_sin_cos="_cos_",
        X2_sin_cos="_sin_",
        LA_sin_cos="_sin_",
        X1X2_deg=5,
        LA_deg=5,
        sgrid=dict(
            nelems=10,
        ),
        phiedge=-phiedge,
    )
    if lasym:
        for key in ["X1", "X2", "LA"]:
            default_parameters[f"{key}_sin_cos"] = "_sincos_"
    conversion_parameters = (
        gvec.util.CaseInsensitiveDict(
            whichInitEquilibrium=1,
            VMECwoutfile=wout_file.absolute(),
            VMECwoutfile_format=VMECwoutfile_format,
            init_LA=False,
        )
        | default_parameters
        | extra_parameters
    )

    if isinstance(rundir, str):
        rundir = Path(rundir)
    if rundir is not None:
        rundir.mkdir(exist_ok=True)

    with tempfile.TemporaryDirectory() as tmpdir:
        logger.info("Reading VMEC wout file...")
        state = gvec.State.new(conversion_parameters, tmpdir)
        ev = state.evaluate(
            "p", "iota", rho=np.sqrt(np.linspace(0, 1, 101)), theta=None, zeta=None
        )
        boundary_parameters = gvec.util.read_parameters(
            Path(tmpdir) / "vmec_to_gvec_boundary_and_axis.txt", format="ini"
        )
        profile_parameters = dict(
            iota=dict(
                type="interpolation",
                rho2=ev.rho.values**2,
                vals=ev.iota.values,
            ),
            pres=dict(
                type="interpolation",
                rho2=ev.rho.values**2,
                vals=ev.p.values,
            ),
        )
        parameters = (
            default_parameters | profile_parameters | boundary_parameters | extra_parameters
        )

        if rundir is not None:
            shutil.copy(sorted(Path(tmpdir).glob("*State*.dat"))[-1], rundir)
            gvec.util.write_parameters(parameters, rundir / "parameter.ini")
            gvec.util.write_parameters(parameters, rundir / "parameter.toml")
    return parameters


def parse_parameter(value: str):
    SEQUENCE = r"[\[\(](.*)[\]\)]"
    INT = r"[-+]?\d+"
    FLOAT = r"[-+]?\d*\.?\d*(?:[eE][-+]?\d+)?"

    if m := re.fullmatch(SEQUENCE, value):
        return tuple(parse_parameter(v) for v in m.group(1).split(","))
    if re.fullmatch(INT, value):
        return int(value)
    if re.fullmatch(FLOAT, value):
        return float(value)
    if value.lower() in ["true", "t"]:
        return True
    if value.lower() in ["false", "f"]:
        return False
    return value


@gvec.errors.without_traceback
def main(args: Sequence[str] | argparse.Namespace | None = None):
    if not isinstance(args, argparse.Namespace):
        args = parser.parse_args(args)

    gvec.util.logging_setup()
    logger = logging.getLogger("gvec")
    if args.quiet:
        logging.disable()
    elif args.verbose >= 2:
        logger.setLevel(logging.DEBUG)
    elif args.verbose == 1:
        logger.setLevel(logging.INFO)
    logger.debug(f"parsed args: {args}")

    if args.VMECwoutfile_format is None:
        args.VMECwoutfile_format = 0  # default to netCDF format

    params = gvec.util.CaseInsensitiveDict()
    for param in args.param:
        if "=" not in param:
            raise ValueError(f"invalid parameter format: '{param}', expected 'key=value'")
        key, value = param.split("=", 1)
        params[key] = parse_parameter(value)
    params = gvec.util.stack_parameters(params)

    convert_vmec_wout(
        wout_file=args.wout_file,
        rundir=args.rundir,
        VMECwoutfile_format=args.VMECwoutfile_format,
        extra_parameters=params,
    )


if __name__ == "__main__":
    exit(main())
