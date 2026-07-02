# Copyright (c) 2025 GVEC Contributors, Max Planck Institute for Plasma Physics
# License: MIT
import contextlib
import logging
import os
from pathlib import Path

import numpy as np

logger = logging.getLogger(__name__)


def get_compile_options() -> str:
    try:
        from gvec import _compile_options as opts
    except ImportError:
        return "UNKNOWN BUILD"
    config = opts.CMAKE_BUILD_TYPE
    if opts.USE_OPENMP:
        config += ", OpenMP"
    if opts.USE_MPI:
        config += ", MPI"
    if opts.GVEC_FIX_HMAP:
        config += f", {opts.GVEC_FIX_HMAP}"
    config += f", {opts.CMAKE_Fortran_COMPILER.name}"
    if len(opts.CMAKE_HOSTNAME) > 0:
        config += f" on {opts.CMAKE_HOSTNAME}"
    return config


def version_info() -> str:
    import platform

    from gvec._version import __version__

    return f"pyGVEC v{__version__} ({get_compile_options()}, python {platform.python_version()}) from {Path(__file__).parent}"


def logging_setup():
    """Setup default logging configuration for GVEC."""
    import logging

    logging.basicConfig(
        format="{levelname:7s} {message}",
        style="{",
        level=logging.WARNING,
    )
    logging.captureWarnings(True)


@contextlib.contextmanager
def chdir(target: Path | str):
    """
    Contextmanager to change the current working directory.

    Using a context has the benefit of automatically changing back to the original directory when the context is exited, even if an exception is raised.
    """
    target = Path(target)
    source = Path.cwd()

    try:
        os.chdir(target)
        yield
    finally:
        os.chdir(source)


def compute_FD(f: np.ndarray, pos, coefs, axis=0):
    """
    1D Finite difference of a function f on equispaced n-dimnesional grid, using FD coefficients coefficients ``coefs`` and relative integer positions to the central evaluation point ``pos``, along one given axis.

    Warnings
    --------
    - if data is periodic, meaning that endpoints of periodic interval are excluded, result can be on all points.
    - If data is not periodic, the result at the **boundaries is WRONG**, for the points ``|min(pos)|`` on the left and ``max(pos)`` on the right along the given axis.

    Parameters
    ----------
    f : numpy.ndarray
        function values on equispaced n-dimnesional grid
    pos : int
        relative integer positions to the central evaluation point, as 1d list or 1d array of integers
    coefs : float
        FD coefficients for each position, as 1d list or 1d array, same size as pos!
    axis  : int, optional
        axis along which the FD is computed, default is 0

    Returns
    -------
    df    : numpy.ndarray
        Finite-Difference result, same shape as f (see warning above!)

    Examples
    --------
    - examples for first derivative of f:
        - 1st order forward FD: ``pos=[1,0]; coefs=[-1,1]/(dx)``
        - 2nd order central FD: ``pos=[-1,1]; coefs=[-1,1]/(2*dx)``
        - 4th order central FD: ``pos=[-2,-1,1,2]; coefs=[1/12,-2/3,2/3,-1/12]/(dx)``
        - 6th order central FD: ``pos=[-3,-2,-1,1,2,3]; coefs=[-1/60,3/20,-3/4,3/4,-3/20,1/60]/(dx)``
        - 8th order central FD: ``pos=[-4,-3,-2,-1,1,2,3,4]; coefs=[1/280,-4/105,1/5,-4/5,4/5,-1/5,4/105,-1/280]/(dx)``
    - examples for second derivatives of f:
        - 2nd order central FD: ``pos=[-1,0,1]; coefs=[1,-2,1]/(dx**2)``
        - 4th order central FD: ``pos=[-2,-1,1,2]; coefs=[-1/12, 4/3,-5/2, 4/3,-1/12]/(dx**2)``
        - 6th order central FD: ``pos=[-3,-2,-1,1,2,3]; coefs=[1/90,-3/20,3/2,-49/18,3/2,-3/20,1/90]/(dx**2)``
        - 8th order central FD: ``pos=[-4,-3,-2,-1,1,2,3,4]; coefs=[-1/560,8/315,-1/5,8/5,-205/72,8/5,-1/5,8/315,-1/560]/(dx**2)``

    """
    if not axis < f.ndim:
        raise ValueError(f"array does not have the requested dimension {axis}")
    if not len(pos) == len(coefs):
        raise ValueError(f"pos and coefs must have the same length, got {pos} and {coefs}")
    df = np.roll(f, -pos[0], axis=axis) * coefs[0]
    for roll, c in zip(pos[1:], coefs[1:]):
        df += np.roll(f, -roll, axis=axis) * c
    return df
