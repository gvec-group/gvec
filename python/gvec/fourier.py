# ============================================================================================================================== #
# Copyright (C) 2024 Robert Babin <robert.babin@ipp.mpg.de>
#
# This file is part of GVEC. GVEC is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3
# of the License, or (at your option) any later version.
#
# GVEC is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
# of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License v3.0 for more details.
#
# You should have received a copy of the GNU General Public License along with GVEC. If not, see <http://www.gnu.org/licenses/>.
# ============================================================================================================================== #
r"""pyGVEC postprocessing - Fourier representation

This module provides functions for computing the Fourier transform in 1D and 2D.
In this context, the Fourier series is of the form :math:`x(\theta, \zeta) = \sum_{m, n} c_{m, n} \cos(m \theta - n \zeta) + s_{m, n} \sin(m \theta - n \zeta)`.
"""

# === Imports === #

from typing import Iterable
import logging

import numpy as np
import xarray as xr


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
        raise ValueError(
            "Expects an odd number of toroidal harmonics: [0 ... N, -N ... -1]"
        )
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
    theta, zeta = np.broadcast_arrays(theta, zeta)
    x = np.zeros_like(theta)
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
    for m in ms:
        for n in ns:
            x += c[m, n] * np.cos(m * theta - n * nfp * zeta)
            x += s[m, n] * np.sin(m * theta - n * nfp * zeta)
    return x


def ev2ft(ev, quiet=False):
    m, n = None, None
    data = {}

    if "N_FP" not in ev.data_vars and not quiet:
        logging.warning("recommended quantity 'N_FP' not found in the provided dataset")

    for var in ev.data_vars:
        if ev[var].dims == ():  # scalar
            data[var] = ((), ev[var].data.item(), ev[var].attrs)

        elif ev[var].dims == ("rad",):  # profile
            data[var] = ("rad", ev[var].data, ev[var].attrs)

        elif {"pol", "tor"} <= set(ev[var].dims) <= {"rad", "pol", "tor"}:
            if "rad" in ev[var].dims:
                vft = []
                for r in ev.rad:
                    vft.append(fft2d(ev[var].sel(rad=r).transpose("pol", "tor").data))
                vcos, vsin = map(np.array, zip(*vft))
                dims = ("rad", "m", "n")
            else:
                vcos, vsin = fft2d(ev[var].transpose("pol", "tor").data)
                dims = ("m", "n")

            if m is None:
                m, n = fft2d_modes(vcos.shape[-2] - 1, vcos.shape[-1] // 2, grid=False)

            attrs = {
                k: v for k, v in ev[var].attrs.items() if k in {"long_name", "symbol"}
            }
            data[f"{var}_mnc"] = (
                dims,
                vcos,
                dict(
                    long_name=f"{ev[var].long_name}, cosine coefficient",
                    symbol=f"{{{ev[var].symbol}}}_{{mn}}^c",
                )
                | attrs,
            )
            data[f"{var}_mns"] = (
                dims,
                vsin,
                dict(
                    long_name=f"{ev[var].long_name}, sine coefficient",
                    symbol=f"{{{ev[var].symbol}}}_{{mn}}^s",
                )
                | attrs,
            )

        elif "xyz" in ev[var].dims and not quiet:
            logging.info(f"skipping quantity '{var}' with cartesian components")

        elif not quiet:
            logging.info(f"skipping quantity '{var}' with dims {ev[var].dims}")

    coords = dict(
        rho=(ev.rho.dims, ev.rho.data, ev.rho.attrs),
        m=(
            "m",
            m if m is not None else [],
            dict(long_name="poloidal mode number", symbol="m"),
        ),
        n=(
            "n",
            n if n is not None else [],
            dict(long_name="toroidal mode number", symbol="n"),
        ),
    )

    ft = xr.Dataset(data, coords=coords)
    ft.attrs["fourier series"] = (
        "Assumes a fourier series of the form 'v(r, θ, ζ) = Σ v^c_mn(r) cos(m θ - n N_FP ζ) + v^s_mn(r) sin(m θ - n N_FP ζ)'"
    )
    return ft
