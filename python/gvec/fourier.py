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


def fft1d(x: Iterable, angle0=0.0):
    """
    Compute the Fourier transform of a 1D array.

    Parameters
    ----------
    x
        Input array to transform, assumed to be sampled on `angle=angle0+np.linspace(0,2*np.pi,len(x),endpoint=False)`.
    angle0: starting value of angle, where x was sampled, in [0,2pi] , defaults to 0.

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
    if angle0 != 0.0:
        ks = np.arange(xf.shape[0])
        xf *= np.exp(-1j * angle0 * ks)

    c = xf.real
    c[1:] *= 2
    s = -2 * xf.imag
    s[0] = 0
    return c, s


def shift_1d(y: np.ndarray, x0, axis, newshape=None):
    """
    shift periodic data along one given axis (and upsample):
    from
    $y(x_i)$ with $x_i=x_0+2\pi i/N$, $i=0,..N-1,$ `N=len(y)`
    to
    $y(x_j)$ with $x_j=0+2\pi j/M$, $j=0,..M-1$

    with `M>=N`.
    can be used to upsample, by setting `x0=0`.
    Inputs:
        y: periodic data
        x0: origin position where y[0] was evaluated
        axis: axis along which to shift
        newshape: output shape along the shifted axis >= input shape. Defaults to input shape
    """
    c = np.fft.rfft(y, norm="forward", axis=axis)
    ks = np.expand_dims(
        np.arange(c.shape[axis]), axis=tuple(i for i in np.arange(y.ndim) if i != axis)
    )
    if newshape is None:
        newshape = y.shape[axis]
    else:
        assert newshape >= y.shape[axis], f"new shape must be >= y.shape[axis]={y.shape[axis]}"
    cshft = c * np.exp(-1j * ks * x0)
    yshft = np.fft.irfft(cshft, newshape, norm="forward", axis=axis)
    return yshft


def fft2d(x: np.ndarray, theta0=0.0, zeta0=0.0):
    r"""
    Compute the Fourier transform of a 2D array.

    The Fourier series is of the form :math:`x(\theta, \zeta) = \sum_{m, n} c_{m, n} \cos(m \theta - n \zeta) + s_{m, n} \sin(m \theta - n \zeta)`.
    The coefficients are given as arrays of shape (M + 1, 2 * N + 1), where M and N are the maximum poloidal and toroidal harmonics, respectively.
    The coefficients with toroidal indices :math:`n > N` are to be interpreted negatively, counted from the end of the array.

    Parameters
    ----------
    x
        Input array of shape (ntheta, nzeta) to transform,
        assumed to be sampled on `theta=theta0+np.linspace(0,2*np.pi,ntheta,endpoint=False)`
        and `zeta=zeta0+np.linspace(0,2*np.pi,nzeta*nfp,endpoint=False)`
    theta0 : starting value of theta, where x was sampled , defaults to 0.
    zeta0 : starting value of zeta, where x was sampled, defaults to 0.
    Returns
    -------
    c : ndarray
        Cosine coefficients of the double-angle Fourier series.
    s : ndarray
        Sine coefficients of the double-angle Fourier series.
    """
    x = np.asarray(x)
    xf = np.fft.rfft2(x.T, norm="forward").T
    if theta0 != 0.0:
        kt = np.arange(
            xf.shape[0]
        )  # theta is second dimension in rfft2, has only positive modes!
        # shift frequencies to get the same modes as if sampled with theta0=0
        xf *= np.exp(-1j * theta0 * kt[:, None])
    if zeta0 != 0.0:
        nzeta = x.shape[1]
        kz = nzeta * np.fft.fftfreq(
            nzeta
        )  # zeta is first dimension in rfft2, has positive and negative modes
        # shift frequencies to get the same modes as if sampled with zeta0=0
        xf *= np.exp(-1j * zeta0 * kz[None, :])

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
        The values of the series evaluated at `theta=np.linspace(0,2*np.pi,2*M+1,endpoint=False)` and `zeta=np.linspace(0,2*np.pi,N*nfp+1,endpoint=False)`.

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
    Scale/Cutoff the coefficients of a 2D Fourier series to a new maximum poloidal and toroidal harmonics.

    Parameters
    ----------
    c : numpy.ndarray
        The coefficients of the original Fourier series, with poloidal mode numbers in the first dimension and toroidal mode numbers in the second. second dimension must be odd.
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


def real_dft_mat(x_in, x_out, nfp=1, modes=None, deriv=0):
    """
    Flexible Direct Fourier Transform for real data
    takes an input array of equidistant points in [0,2pi/nfp[ (exclude endpoint!),
    evaluate the discrete fourier transform with the given 1d mode vector (all >=0) using the input points x_in, then evaluate the inverse transform (or its derivative deriv>0) on the output points x_out anywhere...
    len(x_in) must be > 2*max(modes)
    output is the matrix that transforms real function to real function [derivative]:
     f^deriv(x_out) = Mat f(x_in) (can then be used to do 2d transforms with matmul!)

    nfp is the number of field periods, default 1 (int), all modes are multiples of nfp

    """
    if modes is None:
        modes = np.arange((len(x_in) - 1) // 2 + 1)  # all modes up to Nyquist
    assert np.allclose(x_in[-1] + (x_in[1] - x_in[0]) - x_in[0], 2 * np.pi / nfp)
    assert np.all(modes >= 0), "modes must be positive"
    zeromode = np.where(modes == 0)
    assert len(zeromode) <= 1, "only one zero mode allowed"
    maxmode = np.amax(modes)
    assert len(x_in) > 2 * maxmode, (
        f"number of sampling points ({len(x_in)}) > 2*maxmodenumber ({maxmode})"
    )
    # matrix for forward transform
    Fmat = np.exp(1j * nfp * (modes[:, None] * x_in[None, :]))
    mass_re = Fmat.real @ Fmat.real.T
    mass_im = Fmat.imag @ Fmat.imag.T
    diag_re = np.copy(np.diag(mass_re))
    diag_im = np.copy(np.diag(mass_im))

    assert np.all(np.abs(mass_re - np.diag(diag_re)) < 1.0e-8), "massre must be diagonal"
    assert np.all(np.abs(mass_im - np.diag(diag_im)) < 1.0e-8), "massim must be diagonal"
    diag_im[zeromode] = 1  # imag (=sin) is zero at zero mode
    assert np.all(diag_re > 0.0)
    assert np.all(diag_im > 0.0)

    # inverse mass matrix applied (for real and imag)
    Fmat_mod = np.diag(1 / diag_re) @ Fmat.real + np.diag(1j / diag_im) @ Fmat.imag
    Bmat = get_B_dft(x_out=x_out, nfp=nfp, modes=modes, deriv=deriv)
    Mat = (Bmat @ Fmat_mod).real
    return {
        "F": Fmat_mod,
        "B": Bmat,
        "BF": Mat,
        "modes": modes,
        "x_in": x_in,
        "x_out": x_out,
        "deriv": deriv,
        "nfp": nfp,
    }


def get_B_dft(x_out, deriv, nfp, modes):
    """
    get the matrix B for Fourier transform from modes -> points, with derivative
    """
    modes_back = np.exp(-1j * nfp * (modes[None, :] * x_out[:, None]))
    if deriv > 0:
        modes_back *= (-1j * nfp * modes[None, :]) ** deriv
    return modes_back
