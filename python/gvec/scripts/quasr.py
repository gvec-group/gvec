# Copyright (c) 2025 GVEC Contributors, Max Planck Institute for Plasma Physics
# License: MIT
"""The pyGVEC script for converting a QUASR configuration to a G-Frame for use with GVEC.

QUASR is the A QUAsi-symmetric Stellarator Repository: https://quasr.flatironinstitute.org/

The algorithm is described in the paper Hindenlang et al. DOI: 10.1088/1361-6587/adba11 and is as follows:

### STEP 1: Evaluate surface in cartesian space:

1. using a json file from QUASR (from `https://quasr.flatironinstitute.org/`) with the simsopt interface
2. Evaluate the surface cartesian position $(x,y,z)(\vartheta_i,\zeta_j)$ at a meshgrid  on the full torus:

   $\vartheta_i=2\pi \frac{i}{n_t},i=0\dots,n_t-1,\quad \zeta_j=2\pi\frac{j}{n_z},j=0,\dots,n_z-1$

   where the angles $\vartheta,\zeta$ are just the parametrization of the given surface. (For the quasr surfaces, its a boozer angle parameterization!).

   The number $n_z$ is chosen as a multiple of the number of field periods $n_{FP}$, to be able to reduce the discrete dataset exactly to one field period.

### STEP 2: Project to surface with elliptical cross-sections

Project ${\bm x}_m(\zeta)=\frac{1}{2\pi}\int_{\vartheta=0}^{2\pi}{\bm x}(\vartheta,\zeta)\sigma_m(\vartheta)d\vartheta$ with $\sigma_0(\vartheta)=1,\sigma_{s}(\vartheta)=2\sin(\vartheta),\sigma_{c}(\vartheta)=2\cos(\vartheta)$, leading to a surface

${\bm x}_m(\vartheta,\zeta)={\bm x}_0(\zeta) + {\bm x}_s(\zeta)\sin(\vartheta)+ {\bm x}_c(\zeta)\cos(\vartheta)$

where cross-sections of $\zeta=\text{const.}$ are planar elliptical curves.


### STEP 3: Compute the plane of the ellipse cross-sections

First choice is to set the first basis unit vector $N$ from the center point of the ellipse to a point on the boundary at $\vartheta=0$ position. Then use a first guess for the second basis unit vector $B$ from the center point to $\theta=\frac{\pi}{2}$ position to span the unit normal of the plane $K=(N \times B)$, and then set the second unit vector $B=K\times N$, such that $N$ and $B$ are orthonormal and describe the plane of the ellipse.
###  STEP 4: compute fourier coefficients of the ellipse

The ellipse in a single $N,B$ plane is defined as $  $X^k(\vartheta)= x^k_c\cos(\vartheta)+x^k_s\sin(\vartheta),k=1,2$

We can deduce from the four coefficients the shift  $\vartheta_0$ and the rotation angle $\Gamma$.

### STEP 5: Final frame

The final frame is obtained by rotating the $N$ and $B$ vectors by $-(\Gamma-\vartheta_0)$, which yields constant rotation speed along $\zeta$.
Thus, in the final frame, the rotating ellipse is represented by a single poloidal and a single toroidal Fourier mode.

### STEP 6: cut the original surface with the planes of the frame

For each discrete $N,B$ plane, we compute the intersection of the all curves $\bm x(\vartheta_i,\zeta)$ and compute its position $X^1,X^2$ in the $N,B$ plane. This gives the final surface.
"""

import argparse
from pathlib import Path
import requests
import shutil
from typing import Literal
from collections.abc import Sequence
import logging

import numpy as np
import xarray as xr

from gvec.util import logging_setup
from gvec import gframe

# === Argument Parser === #

parser = argparse.ArgumentParser(
    prog="pygvec-load-quasr",
    description="Load a QUASR configuration and convert it to a G-Frame and boundary for use with GVEC.",
    usage="%(prog)s [-h] (ID | -s FILE | -f FILE) [-v | -q] [--nt NT] [--nz NZ] [--clean CLEANTOL] [--stellsym] [--cutoff=CUTOFF] [--tol=TOL] [--yaml | --toml] [--save-xyz]",
)
parser.add_argument("ID", type=int, nargs="?", help="ID of the QUASR configuration")
parser.add_argument(
    "-s",
    "--simsopt",
    type=Path,
    metavar="FILE",
    help="SIMSOPT JSON file of the boundary (e.g. QUASR configuration)",
)
parser.add_argument("-f", "--file", type=Path, help="netCDF file containing boundary data")
verbosity = parser.add_mutually_exclusive_group()
verbosity.add_argument(
    "-v",
    "--verbose",
    action="count",
    default=0,
    help="verbosity level: -v for info, -vv for debug",
)
verbosity.add_argument("-q", "--quiet", action="store_true", help="suppress output")
parser.add_argument(
    "--nt",
    type=int,
    help="number of theta points (only for ID or -s). Default is 81.",
)
parser.add_argument(
    "--nz",
    type=int,
    help="number of zeta points (only for ID or -s). Default is 81.",
)
parser.add_argument(
    "--clean",
    type=float,
    default=0.0,
    help="tolerance to reduce necessary Fourier modes (M, N) for the input surface. Default is 0., which means no cleaning.",
)
parser.add_argument(
    "--stellsym",
    action="store_true",
    help="if set, imposes stellarator symmetry for the input surface. Use this with great care!",
)
parser.add_argument(
    "--cutoff",
    type=int,
    default=-1,
    help="cutoff toroidal mode number only for G-Frame construction, reduces the number of Fourier modes in the G-Frame. Default is -1, which means no cutoff.",
)
parser.add_argument(
    "--tol",
    type=float,
    default=1e-8,
    help="tolerance for determining minimal necessary (M, N) for the output Fourier modes of X1,X2. default is 1e-8",
)
param_type = parser.add_mutually_exclusive_group()
param_type.add_argument(
    "--yaml",
    action="store_const",
    const="yaml",
    dest="param_type",
    help="write GVEC parameterfile in YAML format",
)
param_type.add_argument(
    "--toml",
    action="store_const",
    const="toml",
    dest="param_type",
    help="write GVEC parameterfile in TOML format",
)
parser.add_argument(
    "--save-xyz",
    action="store_true",
    help="save the boundary points to a netCDF file",
)
parser.add_argument(
    "--name",
    type=str,
    default="",
    help="name for outputfiles. If not set, the name of the input boundary is used.",
)


def check_args(parser, args):
    if sum([args.ID is None, args.simsopt is None, args.file is None]) != 2:
        raise parser.error("exactly one of ID, -s or -f must be provided.")
    if args.ID is not None and (args.ID < 0 or args.ID > 9999999):
        raise parser.error("ID must be between 0 and 9999999.")
    if args.simsopt is not None and not args.simsopt.exists():
        raise parser.error(f"File {args.simsopt} does not exist.")
    if args.file is None:
        if args.nt is None:
            args.nt = 81
        elif args.nt < 1:
            raise parser.error("Number of theta points must be greater than 0.")
        if args.nz is None:
            args.nz = 81
        elif args.nz < 1:
            raise parser.error("Number of zeta points must be greater than 0.")
    else:
        if not args.file.exists():
            raise parser.error(f"File {args.file} does not exist.")
        if args.nt is not None:
            raise parser.error(
                "Number of theta points cannot manually be set with a boundary file."
            )
        if args.nz is not None:
            raise parser.error(
                "Number of zeta points cannot manually be set with a boundary file."
            )
    if args.tol is not None and args.tol <= 0.0:
        raise parser.error("Tolerance must be greater than 0.")
    if args.clean is not None and args.clean < 0.0:
        raise parser.error("Cleaning tolerance must be greater than 0.")
    if args.param_type is None:
        args.param_type = "toml"


# === Functions === #


def get_json_from_quasr(configuration: int, filename: str | Path = None):
    """Retrieve a simsopt-compatible JSON for a given QUASR configuration."""

    if filename is None:
        filename = Path(f"quasr-{configuration:07d}.json")

    url = f"https://quasr.flatironinstitute.org/simsopt_serials/{configuration // 10**3:04d}/serial{configuration:07d}.json"
    with requests.get(url, stream=True) as response, open(filename, "wb") as file:
        if not response.ok:
            raise RuntimeError(
                f"Failed to download QUASR configuration {configuration}: {response.status_code} {response.reason}"
            )
        shutil.copyfileobj(response.raw, file)

    return filename


def get_surface_from_json_file(filename: Path | str):
    """Get the boundary surface as a SIMSOPT Surface object from a QUASR JSON file."""
    from simsopt._core import load

    surfaces, coils = load(filename)
    return surfaces[-1]


def get_xyz_from_surface(nt: int, nz: int, surface):
    """Sample a SIMSOPT Surface object in cartesian coordinates.

    Sample surface at nt,nz*nfp point positions on the full torus.
    Gives cartesian positions xyz[0:nz*nfp,0:nt,0:2].
    """
    nfp = surface.nfp
    # simsopt.Surface objects use [0,1] for theta & zeta
    t1d = np.linspace(0, 1, nt, endpoint=False)
    z1d = np.linspace(0, 1, nz * nfp, endpoint=False)
    t, z = np.meshgrid(t1d, z1d)

    xyz = np.zeros((nz * nfp, nt, 3))
    surface.gamma_lin(xyz, z.flatten(), t.flatten())
    return xyz


def save_xyz(xyz: np.ndarray, nfp: int, filename: Path | str, attrs: dict = {}):
    """
    Save cartesian surface points xyz[0:nz*nfp,0:nt,0:2] to a netcdf file.
    The surface cartesian position $(x,y,z)(\vartheta_i,\zeta_j)$ must be evaluated
    at a meshgrid  on the full torus, excluding the periodic endpoint:
    $\vartheta_i=2\pi \frac{i}{n_t},i=0\dots,n_t-1,\quad \zeta_j=2\pi\frac{j}{n_z},j=0,\dots,n_z-1$
    Inputs:
    xyz: cartesian surface points xyz[0:nz*nfp,0:nt,0:2]
    nfp: number of field periods
    filename: path to netcdf file (with .nc extension)

    Example usage, with a function that evaluates the surface cartesian position $\vec{x}(\vartheta,\zeta)$,
    provided the number of points `ntheta` and `nzeta` and the number of field periods `nfp`:

    ```python
    theta=np.linspace(0,2*np.pi,ntheta,endpoint=False)
    zeta=np.linspace(0,2*np.pi,nzeta*nfp,endpoint=False)
    xyz=np.zeros((nzeta,ntheta,3))
    for j in range(nzeta):
         for i in range(ntheta):
             xyz[j,i,:] = eval_surface(theta[i],zeta[j])
    save_xyz(xyz,nfp,'mycase-boundary.nc')
    ```

    """
    import datetime
    from gvec import __version__

    ds = xr.Dataset(
        data_vars=dict(
            pos=(("zeta", "theta", "xyz"), xyz),
            nfp=((), nfp),
        ),
        coords=dict(
            xyz=("xyz", ["x", "y", "z"]),
            theta=("theta", np.linspace(0, 2 * np.pi, xyz.shape[1], endpoint=False)),
            zeta=("zeta", np.linspace(0, 2 * np.pi, xyz.shape[0], endpoint=False)),
        ),
        attrs=dict(
            creator="pygvec load-quasr",
            gvec_version=__version__,
            date=str(datetime.datetime.now().date()),
        )
        | attrs,
    )
    ds.to_netcdf(filename, mode="w")


def load_xyz(filename: Path | str):
    ds = xr.load_dataset(filename)
    if "pos" not in ds or "nfp" not in ds:
        raise ValueError(
            f"File {filename} does not contain the required 'pos' and 'nfp' variables."
        )
    if set(ds.pos.dims) != {"zeta", "theta", "xyz"}:
        raise ValueError(
            f"File {filename} does not contain the required dimensions 'zeta', 'theta', and 'xyz'."
        )
    if ds.zeta.size % ds.nfp.item() != 0:
        raise ValueError(
            f"length of zeta ({ds.zeta.size}) is not compatible with nfp {ds.nfp.item()}. "
            "It must be a multiple of nfp."
        )
    xyz = ds.pos.transpose("zeta", "theta", "xyz").values
    nfp = ds.nfp.item()
    return xyz, nfp


# === Script === #


def main(args: Sequence[str] | argparse.Namespace | None = None):
    if isinstance(args, argparse.Namespace):
        pass
    else:
        args = parser.parse_args(args)
    check_args(parser, args)

    logging_setup()
    logger = logging.getLogger(__name__)
    if args.quiet:
        logging.disable()
    elif args.verbose >= 2:
        logger.setLevel(logging.DEBUG)
    elif args.verbose == 1:
        logger.setLevel(logging.INFO)
    logger.debug(f"parsed args: {args}")

    if args.ID is not None:
        logger.info("Downloading QUASR configuration")
        try:
            filename = get_json_from_quasr(args.ID)
        except RuntimeError as e:
            logger.error(e)
            return 1
    elif args.simsopt is not None:
        filename = args.simsopt

    if args.ID is not None or args.simsopt is not None:
        logger.info("Loading SIMSOPT surface")
        surface = get_surface_from_json_file(filename)
        nfp = surface.nfp
        xyz = get_xyz_from_surface(args.nt, args.nz, surface)
        name = str(filename.stem)
    else:
        filename = args.file
        logger.info(f"Reading boundary file {filename}")
        xyz, nfp = load_xyz(filename)
        if str(filename.stem).endswith("-boundary"):
            name = str(filename.stem)[:-9]
        else:
            name = str(filename.stem)

    if args.name != "":
        name = args.name

    if args.save_xyz:
        logger.info("Saving boundary points to netCDF file")
        save_xyz(xyz, nfp, f"{name}-boundary.nc", attrs={"source": str(filename)})

    gframe.construct_gframe_from_surface(
        xyz,
        nfp,
        name,
        tolerance_output=args.tol,
        format=args.param_type,
        tolerance_clean_surface=args.clean,
        impose_stell_symmetry=args.stellsym,
        cutoff_gframe=args.cutoff,
        logger=logger,
    )


if __name__ == "__main__":
    exit(main())
