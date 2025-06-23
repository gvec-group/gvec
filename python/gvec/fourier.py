# Copyright (c) 2025 GVEC Contributors, Max Planck Institute for Plasma Physics
# License: MIT
r"""pyGVEC postprocessing - Fourier representation

This module provides functions for computing the Fourier transform in 1D and 2D.
In this context, the Fourier series is of the form :math:`x(\theta, \zeta) = \sum_{m, n} c_{m, n} \cos(m \theta - n \zeta) + s_{m, n} \sin(m \theta - n \zeta)`.
"""

# === Imports === #

from collections.abc import Iterable

import numpy as np


# === Transform functions === #


def fft1d(x: Iterable):
    """
    Compute the Fourier transform of a 1D array.

    Parameters
    ----------
    x
        Input array to transform.

    Returns
    -------
    c : ndarray
        Cosine coefficients of the Fourier series.
    s : ndarray
        Sine coefficients of the Fourier series.

    Notes
    -----
    The function uses the real-input fast Fourier transform (rfft) from numpy.
    """
    x = np.asarray(x)
    xf = np.fft.rfft(x, norm="forward")
    c = xf.real
    c[1:] *= 2
    s = -2 * xf.imag
    s[0] = 0
    return c, s


def fft2d(x: np.ndarray):
    r"""
    Compute the Fourier transform of a 2D array.

    The Fourier series is of the form :math:`x(\theta, \zeta) = \sum_{m, n} c_{m, n} \cos(m \theta - n \zeta) + s_{m, n} \sin(m \theta - n \zeta)`.
    The coefficients are given as arrays of shape (M + 1, 2 * N + 1), where M and N are the maximum poloidal and toroidal harmonics, respectively.
    The coefficients with toroidal indices :math:`n > N` are to be interpreted negatively, counted from the end of the array.

    Parameters
    ----------
    x
        Input array of shape (m, n) to transform.

    Returns
    -------
    c : ndarray
        Cosine coefficients of the double-angle Fourier series.
    s : ndarray
        Sine coefficients of the double-angle Fourier series.
    """
    x = np.asarray(x)
    xf = np.fft.rfft2(x.T, norm="forward").T

    N = (x.shape[1] - 1) // 2
    c = 2 * xf.real
    c[0, 0] /= 2  # no double counting for n = 0
    c = np.roll(c, -1, axis=1)[:, ::-1]  # invert the toroidal indices
    c[0, -N:] = 0  # zero out the negative toroidal indices

    s = -2 * xf.imag
    s = np.roll(s, -1, axis=1)[:, ::-1]  # invert the toroidal indices
    s[0, -N:] = 0  # zero out the negative toroidal indices

    if x.shape[1] % 2 == 0:
        # remove the extra toroidal harmonic if the input has an even number of points
        c = np.concatenate([c[:, : N + 1], c[:, -N:]], axis=1)
        s = np.concatenate([s[:, : N + 1], s[:, -N:]], axis=1)

    return c, s


def ifft2d(c: np.ndarray, s: np.ndarray, deriv: str | None = None, nfp: int = 1) -> np.ndarray:
    """
    Inverse Fast-Fourier-Transform of a 2D Fourier series.

    Parameters
    ----------
    c : numpy.ndarray
        Cosine coefficients of the Fourier series.
    s : numpy.ndarray
        Sine coefficients of the Fourier series.
    deriv : str, optional
        Derivative to evaluate, by default None. Specified as 'theta', 'zeta' or any string of 't' & 'z', e.g. 't', 'tz', 'ttz', ...
    nfp : int, optional
        Number of field periods, by default 1. Only used for derivatives, the data itself is always assumed to be in a single field period.

    Returns
    -------
    x : numpy.ndarray
        The values of the series at the given angles.
    """
    if c.shape != s.shape:
        raise ValueError("c and s must have the same shape")
    M = c.shape[0] - 1
    N = c.shape[1] // 2
    c = np.asarray(c)
    s = np.asarray(s)

    if deriv is not None:
        mg, ng = fft2d_modes(c.shape[0] - 1, c.shape[1] // 2, grid=True)
        ng *= nfp
        if set(deriv) <= {"t", "z"}:
            ts, zs = deriv.count("t"), deriv.count("z")
            for _ in range(ts):
                c, s = mg * s, -mg * c
            for _ in range(zs):
                c, s = -ng * s, ng * c
        elif deriv == "theta":
            c, s = mg * s, -mg * c
        elif deriv == "zeta":
            c, s = -ng * s, ng * c
        else:
            raise ValueError(
                f"Invalid derivative specification, got '{deriv}', expected 'theta', 'zeta', 't', 'z', 'tt', 'tz', ..."
            )

    c = np.roll(c[:, ::-1], 1, axis=1)
    c[0, :] *= 2
    c = c / 2

    s = np.roll(s[:, ::-1], 1, axis=1)
    s[0, :] *= 2
    s = -s / 2

    xf = c + 1j * s
    # always use an odd number of points in both directions
    x = np.fft.irfft2(xf.T, s=(2 * N + 1, 2 * M + 1), norm="forward").T
    return x


def fft2d_modes(M: int, N: int, grid: bool = False):
    """
    Generate the modenumbers for a 2D FFT, as performed by `fft2d`.

    Parameters
    ----------
    M : int
        The maximum poloidal modenumber.
    N : int
        The maximum toroidal modenumber.

    Returns
    -------
    m : numpy.ndarray
        The poloidal modenumbers.
    n : numpy.ndarray
        The toroidal modenumbers.
    """
    m = np.arange(M + 1)
    n = np.concatenate([np.arange(N + 1), np.arange(-N, 0)])
    if grid:
        m, n = np.meshgrid(m, n, indexing="ij")
    return m, n


def scale_modes2d(c, M, N):
    """
    Scale the coefficients of a 2D Fourier series to a new maximum poloidal and toroidal harmonics.

    Parameters
    ----------
    c : numpy.ndarray
        The coefficients of the original Fourier series.
    M : int
        The new maximum poloidal harmonic.
    N : int
        The new maximum toroidal harmonic.

    Returns
    -------
    c2 : numpy.ndarray
        The coefficients of the scaled Fourier series.
    """
    if c.shape[1] % 2 != 1:
        raise ValueError("Expects an odd number of toroidal harmonics: [0 ... N, -N ... -1]")
    M1, N1 = c.shape[0] - 1, c.shape[1] // 2
    m1, n1 = fft2d_modes(M1, N1, grid=True)
    m2, n2 = fft2d_modes(M, N, grid=True)
    Mmin, Nmin = min(M1, M), min(N1, N)

    c2 = np.zeros((M + 1, 2 * N + 1), dtype=c.dtype)
    c2[(m2 <= Mmin) & (np.abs(n2) <= Nmin)] = c[(m1 <= Mmin) & (np.abs(n1) <= Nmin)]
    return c2


def eval2d(
    c: np.ndarray,
    s: np.ndarray,
    theta: np.ndarray,
    zeta: np.ndarray,
    deriv: str | None = None,
    nfp: int = 1,
):
    """
    Evaluate a 2D Fourier series at given poloidal and toroidal angles.

    Parameters
    ----------
    c : numpy.ndarray
        Cosine coefficients of the Fourier series.
    s : numpy.ndarray
        Sine coefficients of the Fourier series.
    theta : numpy.ndarray
        Poloidal angles at which to evaluate the series.
    zeta : numpy.ndarray
        Toroidal angles at which to evaluate the series.
    deriv : str, optional
        Derivative to evaluate, by default None. Specified as 'theta', 'zeta' or any string of 't' & 'z', e.g. 't', 'tz', 'ttz', ...
    nfp : int, optional
        Number of field periods, by default 1.

    Returns
    -------
    x : numpy.ndarray
        The values of the series at the given angles.
    """
    if theta.shape != zeta.shape:
        raise ValueError("theta and zeta must have the same shape")

    shape = theta.shape
    theta, zeta = theta.ravel(), zeta.ravel()

    if deriv is not None:
        mg, ng = fft2d_modes(c.shape[0] - 1, c.shape[1] // 2, grid=True)
        ng *= nfp
        if set(deriv) <= {"t", "z"}:
            ts, zs = deriv.count("t"), deriv.count("z")
            for _ in range(ts):
                c, s = mg * s, -mg * c
            for _ in range(zs):
                c, s = -ng * s, ng * c
        elif deriv == "theta":
            c, s = mg * s, -mg * c
        elif deriv == "zeta":
            c, s = -ng * s, ng * c
        else:
            raise ValueError(
                f"Invalid derivative specification, got '{deriv}', expected 'theta', 'zeta', 't', 'z', 'tt', 'tz', ..."
            )

    ms, ns = fft2d_modes(c.shape[0] - 1, c.shape[1] // 2)
    x = np.zeros_like(theta)
    for m in ms:
        for n in ns:
            # this python double loop is NOT slower than numpy array operations
            x += c[m, n] * np.cos(m * theta - n * nfp * zeta)
            x += s[m, n] * np.sin(m * theta - n * nfp * zeta)
    return x.reshape(shape)
