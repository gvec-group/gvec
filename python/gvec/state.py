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
"""pygvec postprocessing"""

from ._fgvec import modgvec_py_post as _post

from pathlib import Path
from typing import Mapping, Callable, Iterable, Literal
import re
import inspect
import functools
import tempfile

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable


def _assert_init(func):
    @functools.wraps(func)
    def wrapped(self, *args, **kwargs):
        if not self.initialized:
            raise RuntimeError("State is not initialized.")
        if not _post.initialized:
            raise RuntimeError("State is initialized, but GVEC libaray is not!")
        return func(self, *args, **kwargs)

    return wrapped


def _evaluate_1D_factory(
    func: callable, argnames: Iterable[str], n_out: int, vector_out: bool = False
):
    params = [inspect.Parameter("self", inspect.Parameter.POSITIONAL_OR_KEYWORD)] + [
        inspect.Parameter(
            name, inspect.Parameter.POSITIONAL_OR_KEYWORD, annotation=np.ndarray
        )
        for name in argnames
    ]
    returns = tuple[tuple(np.ndarray for _ in range(n_out))]
    sig = inspect.Signature(params, return_annotation=returns)

    @_assert_init
    def wrapper(self, *args, **kwargs):
        bound_args = sig.bind(self, *args, **kwargs)
        inputs = [
            np.asfortranarray(value, dtype=np.float64)
            for key, value in bound_args.arguments.items()
            if key != "self"
        ]
        n = inputs[0].size
        for value in inputs:
            if value.shape != (n,):
                raise ValueError("All arguments must be 1D arrays of the same size.")

        if vector_out:
            outputs = [
                np.zeros((3, n), dtype=np.float64, order="F") for _ in range(n_out)
            ]
        else:
            outputs = [np.zeros(n, dtype=np.float64) for _ in range(n_out)]
        func(n, *inputs, *outputs)
        return outputs

    wrapper.__signature__ = sig
    wrapper.__name__ = func.__name__
    return wrapper


class State:
    # === Constructor & Destructor === #

    def __init__(
        self,
        parameterfile: str | Path,
        statefile: str | Path,
        redirect_stdout: bool = True,
    ):
        self.initialized: bool = False
        self.parameterfile: Path | None = None
        self.statefile: Path | None = None
        self.GB_available: bool = False

        if _post.initialized:
            raise NotImplementedError("Only one instance of State is allowed.")
        if not Path(parameterfile).exists():
            raise FileNotFoundError(f"Parameter file {parameterfile} does not exist.")
        if not Path(statefile).exists():
            raise FileNotFoundError(f"State file {statefile} does not exist.")

        if redirect_stdout:
            self._stdout = tempfile.NamedTemporaryFile(mode="r", prefix="gvec-stdout-")
            _post.redirect_stdout(self._stdout.name)
        self.parameterfile = Path(parameterfile)
        _post.init(parameterfile)
        self.statefile = Path(statefile)
        _post.readstate(statefile)
        self.initialized = True

    @_assert_init
    def finalize(self):
        _post.finalize()
        self.initialized = False

    def __del__(self):
        if hasattr(self, "_stdout"):
            self._stdout.close()
        # silently ignore non-initialized states
        if self.initialized:
            self.finalize()

    # === Context Manager === #

    @_assert_init
    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        # silently ignore non-initialized states
        if self.initialized:
            self.finalize()

    # === Debug Information === #

    def __repr__(self):
        return (
            "<pygvec.State("
            + ",".join(
                [
                    "initialized" if self.initialized else "finalized",
                    self.parameterfile.name,
                    self.statefile.name,
                ]
            )
            + ")>"
        )

    @property
    def stdout(self):
        if not hasattr(self, "_stdout"):
            return None
        self._stdout.seek(0)
        return self._stdout.read()

    # === Evaluation Methods === #

    @property
    @_assert_init
    def nfp(self):
        return _post.nfp

    @_assert_init
    def get_integration_points(self, quantity: str = "LA"):
        if not isinstance(quantity, str):
            raise ValueError("Quantity must be a string.")
        elif quantity not in ["X1", "X2", "LA", "GB"]:
            raise ValueError(
                f"Unknown quantity: {quantity}, expected one of 'X1', 'X2', 'LA', 'GB'."
            )
        elif quantity == "GB" and not self.GB_available:
            raise AttributeError(
                "Boozer potential GB is not available. Run `construct_GB` first."
            )

        r_n, t_n, z_n = _post.get_integration_points_num(quantity)
        r_GP, r_w = (np.zeros(r_n, dtype=np.float64) for _ in range(2))
        t_w, z_w = _post.get_integration_points(quantity, r_GP, r_w)
        return r_GP, r_w, t_n, t_w, z_n, z_w

    @_assert_init
    def evaluate_base_tens(
        self,
        quantity: str,
        derivs: str | None,
        rho: np.ndarray,
        theta: np.ndarray,
        zeta: np.ndarray,
    ):
        if not isinstance(quantity, str):
            raise ValueError("Quantity must be a string.")
        elif quantity not in ["X1", "X2", "LA", "GB"]:
            raise ValueError(
                f"Unknown quantity: {quantity}, expected one of 'X1', 'X2', 'LA', 'GB'."
            )
        elif quantity == "GB" and not self.GB_available:
            raise AttributeError(
                "Boozer potential GB is not available. Run `construct_GB` first."
            )
        if derivs is not None:
            if not isinstance(derivs, str):
                raise ValueError("Derivatives must be a string.")
            if m := re.match(r"(r{0,2})(t{0,2}|z{0,2}|tz)$", derivs):
                sel_derivs = m.groups()
            else:
                raise ValueError(f"Unknown derivative: {derivs}")
        else:
            sel_derivs = ("", "")

        rho = np.asfortranarray(rho, dtype=np.float64)
        theta = np.asfortranarray(theta, dtype=np.float64)
        zeta = np.asfortranarray(zeta, dtype=np.float64)
        if rho.ndim != 1 or theta.ndim != 1 or zeta.ndim != 1:
            raise ValueError("rho, theta, and zeta must be 1D arrays.")
        if rho.max() > 1.0 or rho.min() < 0.0:
            raise ValueError("rho must be in the range [0, 1].")

        result = np.zeros(
            (rho.size, theta.size, zeta.size), dtype=np.float64, order="F"
        )
        _post.evaluate_base_tens(rho, theta, zeta, quantity, *sel_derivs, result)
        return result

    @_assert_init
    def evaluate_base_list_tz(
        self, quantity: str, derivs: str | None, rho: np.ndarray, thetazeta: np.ndarray
    ):
        if not isinstance(quantity, str):
            raise ValueError("Quantity must be a string.")
        elif quantity not in ["X1", "X2", "LA", "GB"]:
            raise ValueError(
                f"Unknown quantity: {quantity}, expected one of 'X1', 'X2', 'LA', 'GB'."
            )
        elif quantity == "GB" and not self.GB_available:
            raise AttributeError(
                "Boozer potential GB is not available. Run `construct_GB` first."
            )
        if derivs is not None:
            if not isinstance(derivs, str):
                raise ValueError("Derivatives must be a string.")
            if m := re.match(r"(r{0,2})(t{0,2}|z{0,2}|tz)$", derivs):
                sel_derivs = m.groups()
            else:
                raise ValueError(f"Unknown derivative: {derivs}")
        else:
            sel_derivs = ("", "")

        rho = np.asfortranarray(rho, dtype=np.float64)
        thetazeta = np.asfortranarray(thetazeta, dtype=np.float64)
        if rho.ndim != 1:
            raise ValueError("rho must be a 1D array.")
        if thetazeta.ndim != 2 or thetazeta.shape[0] != 2:
            raise ValueError("thetazeta must be a 2D array with shape (2, n).")
        if rho.max() > 1.0 or rho.min() < 0.0:
            raise ValueError("rho must be in the range [0, 1].")

        result = np.zeros((rho.size, thetazeta.shape[1]), dtype=np.float64, order="F")
        _post.evaluate_base_list_tz(
            rho.size, thetazeta.shape[1], rho, thetazeta, quantity, *sel_derivs, result
        )
        return result

    @_assert_init
    def evaluate_base_list_tz_all(
        self, quantity: str, rho: np.ndarray, thetazeta: np.ndarray
    ):
        if not isinstance(quantity, str):
            raise ValueError("Quantity must be a string.")
        elif quantity not in ["X1", "X2", "LA", "GB"]:
            raise ValueError(
                f"Unknown quantity: {quantity}, expected one of 'X1', 'X2', 'LA', 'GB'."
            )
        elif quantity == "GB" and not self.GB_available:
            raise AttributeError(
                "Boozer potential GB is not available. Run `construct_GB` first."
            )

        rho = np.asfortranarray(rho, dtype=np.float64)
        thetazeta = np.asfortranarray(thetazeta, dtype=np.float64)
        if rho.ndim != 1:
            raise ValueError("rho must be a 1D array.")
        if thetazeta.ndim != 2 or thetazeta.shape[0] != 2:
            raise ValueError("thetazeta must be a 2D array with shape (2, n).")
        if rho.max() > 1.0 or rho.min() < 0.0:
            raise ValueError("rho must be in the range [0, 1].")

        # Q, dQ_drho, dQ_dtheta, dQ_dzeta, dQ_drr, dQ_drt, dQ_drz, dQ_dtt, dQ_dtz, dQ_dzz
        outputs = [
            np.zeros((rho.size, thetazeta.shape[1]), dtype=np.float64, order="F")
            for _ in range(10)
        ]

        _post.evaluate_base_list_tz_all(
            rho.size, thetazeta.shape[1], rho, thetazeta, quantity, *outputs
        )
        return outputs

    @_assert_init
    def evaluate_base_tens_all(
        self, quantity: str, rho: np.ndarray, theta: np.ndarray, zeta: np.ndarray
    ):
        if not isinstance(quantity, str):
            raise ValueError("Quantity must be a string.")
        elif quantity not in ["X1", "X2", "LA", "GB"]:
            raise ValueError(
                f"Unknown quantity: {quantity}, expected one of 'X1', 'X2', 'LA', 'GB'."
            )
        elif quantity == "GB" and not self.GB_available:
            raise AttributeError(
                "Boozer potential GB is not available. Run `construct_GB` first."
            )

        rho = np.asfortranarray(rho, dtype=np.float64)
        theta = np.asfortranarray(theta, dtype=np.float64)
        zeta = np.asfortranarray(zeta, dtype=np.float64)
        if rho.ndim != 1 or theta.ndim != 1 or zeta.ndim != 1:
            raise ValueError("rho, theta, and zeta must be 1D arrays.")
        if rho.max() > 1.0 or rho.min() < 0.0:
            raise ValueError("rho must be in the range [0, 1].")

        # Q, dQ_drho, dQ_dtheta, dQ_dzeta, dQ_drr, dQ_drt, dQ_drz, dQ_dtt, dQ_dtz, dQ_dzz
        outputs = [
            np.zeros((rho.size, theta.size, zeta.size), dtype=np.float64, order="F")
            for _ in range(10)
        ]

        _post.evaluate_base_tens_all(
            rho.size, theta.size, zeta.size, rho, theta, zeta, quantity, *outputs
        )
        return outputs

    evaluate_hmap = _evaluate_1D_factory(
        _post.evaluate_hmap,
        [
            "X1",
            "X2",
            "zeta",
        ]
        + [f"d{Q}_d{i}" for i in "rtz" for Q in ["X1", "X2"]],
        4,
        True,
    )  # -> pos, e_rho, e_theta, e_zeta

    evaluate_hmap_only = _evaluate_1D_factory(
        _post.evaluate_hmap_only, ["X1", "X2", "zeta"], 4, True
    )  # -> pos, e_X1, e_X2, e_zeta3

    evaluate_metric = _evaluate_1D_factory(
        _post.evaluate_metric,
        [
            "X1",
            "X2",
            "zeta",
        ]
        + [
            f"d{Q}_d{i}"
            for i in "r t z rr rt rz tt tz zz".split()
            for Q in ["X1", "X2"]
        ],
        24,
    )  # -> g_rr, g_rt ... g_zz, dg_rr_dr, dg_rt_dr ... dg_zz_dz

    evaluate_jacobian = _evaluate_1D_factory(
        _post.evaluate_jacobian,
        [
            "X1",
            "X2",
            "zeta",
            "dX1_dr",
            "dX2_dr",
            "dX1_dt",
            "dX2_dt",
            "dX1_dz",
            "dX2_dz",
        ],
        4,
    )  # -> Jac_h, dJac_h_dr, dJac_h_dt, dJac_h_dz

    @_assert_init
    def evaluate_profile(self, quantity: str, rho: np.ndarray):
        if not isinstance(quantity, str):
            raise ValueError("Quantity must be a string.")
        elif quantity not in [
            "iota",
            "iota_prime",
            "p",
            "p_prime",
            "chi",
            "chi_prime",
            "Phi",
            "Phi_prime",
            "Phi_2prime",
            "PhiNorm",
            "PhiNorm_prime",
        ]:
            raise ValueError(f"Unknown quantity: {quantity}")

        rho = np.asfortranarray(rho, dtype=np.float64)
        if rho.ndim != 1:
            raise ValueError("rho must be a 1D array.")
        if rho.max() > 1.0 or rho.min() < 0.0:
            raise ValueError("rho must be in the range [0, 1].")

        result = np.zeros(rho.size, dtype=np.float64, order="F")
        _post.evaluate_profile(rho.size, rho, quantity, result)
        return result

    # === Boozer transform === #

    @_assert_init
    def construct_GB(
        self,
        degree: int,
        m: int,
        n: int,
        fourier: Literal["sin", "cos", "both", "sincos"] = "sin",
        continuity: int | Literal["smooth"] = "smooth",
        method: Literal["interpolation", "projection"] = "projection",
        n_theta: int | None = None,
        n_zeta: int | None = None,
    ):
        """Construct the Boozer potential $G_B$ ('GB') by projection onto a basis.
        Recomputes $\\lambda$ with the new basis to satisfy the integrability condition.

        Parameters
        ----------
        degree
            The polynomial degree of the radial spline basis.
        m
            The maximum poloidal mode number.
        n
            The maximum toroidal mode number.
        fourier
            The type of Fourier series to use for the poloidal and toroidal direction. 'sin', 'cos', 'sincos' / 'both'. Default is 'both'.
        continuity
            The continuity of the basis. Only 'degree - 1' / 'smooth' and '-1' is supported for now. Default is 'smooth'.
        method
            The method for the construction of the radial DoFs: 'interpolation' on the IPs or 'projection' using Gauss points. Default is 'projection'.
        n_theta
            The number of poloidal integration points. Set to 2 * m + 1 per default.
        n_zeta
            The number of toroidal integration points. Set to 2 * n + 1 per default.

        Raises
        ------
        ValueError
            If the type of Fourier series is not one of 'sin', 'cos', 'sincos', or 'both'.
            If the continuity is not either degree - 1 (smooth) or -1.
            If m or n is negative.
            If the degree is less than 1.
        """
        pass
        # --- argument handling --- #
        if fourier == "both":
            fourier = "sincos"
        if fourier not in ["sin", "cos", "sincos"]:
            raise ValueError(
                f"Type of fourier series must be one of 'sin', 'cos', 'sincos' or 'both', got {fourier}."
            )
        fourier = f"{f'_{fourier}_':8s}"
        if continuity == "smooth":
            continuity = degree - 1
        if continuity not in [degree - 1, -1]:
            raise ValueError(
                f"Continuity must be either degree - 1 (smooth) or -1, got {continuity}."
            )
        if m < 0 or n < 0:
            raise ValueError("m and n must be non-negative, got {m} and {n}.")
        if degree < 1:
            raise ValueError("Degree must be larger than 0, got {degree}.")
        if method not in ["interpolation", "projection"]:
            raise ValueError(
                f"Method must be one of 'interpolation' or 'projection', got {method}."
            )
        mn_max = [(m, n)] + [
            _post.get_mn_max(var)
            for var in ["X1", "X2"] + (["GB"] if self.GB_available else [])
        ]
        m_max, n_max = max(list(zip(*mn_max))[0]), max(list(zip(*mn_max))[1])
        if n_theta is not None and n_theta < 2 * m_max + 1:
            raise ValueError(
                f"n_theta must be at least 2 * m_max ({m_max}) + 1, got {n_theta}."
            )
        if n_zeta is not None and n_zeta < 2 * n_max + 1:
            raise ValueError(
                f"n_zeta must be at least 2 * n_max ({n_max}) + 1, got {n_zeta}."
            )

        _, n_theta0, n_zeta0 = _post.get_integration_points_num("LA")
        if n_theta is None:
            n_theta = 2 * m + 1
            if n_theta < n_theta0:
                n_theta = n_theta0
        if n_zeta is None:
            n_zeta = 2 * n + 1
            if n_zeta < n_zeta0:
                n_zeta = n_zeta0

        # --- prepare lambda --- #
        _post.set_integration_points_tz(n_theta, n_zeta)
        _post.change_base_la(degree, continuity, m, n, fourier)
        _post.recompute_la()

        # --- initialize new base --- #
        _post.init_base_gb(degree, continuity, m, n, fourier)

        # --- construction --- #
        if method == "interpolation":
            n_IP = _post.get_s_nbase("GB")
            s_IP = np.zeros(n_IP, dtype=np.float64)
            _post.get_s_ip("GB", s_IP)

            # ev = Evaluations(self, s_IP, "int", "int")
        # elif method == "projection":
        # ev = Evaluations(self, "int", "int", "int")
        ev = NotImplemented

        ev.compute("dGB_dt_def", "dGB_dz_def")
        GBt = np.asfortranarray(
            ev.dGB_dt_def.transpose("rho", "zeta", "theta").values.reshape(
                ev.rho.size, -1
            )
        )
        GBz = np.asfortranarray(
            ev.dGB_dz_def.transpose("rho", "zeta", "theta").values.reshape(
                ev.rho.size, -1
            )
        )

        if method == "interpolation":
            _post.gb_project_f_interpolate_s(GBt, GBz)
        elif method == "projection":
            _post.gb_project_f_project_s(GBt, GBz)
        self.GB_available = True

    # === Integration with computable quantities === #

    def compute(self, ds: xr.Dataset, *quantities):
        from .comp import compute

        return compute(ds, *quantities, state=self)

    # === Plotting Methods === #

    def plot_surfaces(
        self,
        quantities: (
            str
            | Callable[[xr.Dataset], xr.DataArray | tuple[xr.DataArray, str]]
            | Iterable[
                str | Callable[[xr.Dataset], xr.DataArray | tuple[xr.DataArray, str]]
            ]
        ),
        rho: int | tuple[float, float, int] | np.ndarray = 3,
        theta: int | tuple[float, float, int] | np.ndarray = 101,
        zeta: int | tuple[float, float, int] | np.ndarray = 81,
        n_xtick: int = 5,
        n_ytick: int = 5,
        share_colorbar: bool = True,
        figkwargs: Mapping = {},
        ckwargs: Mapping = {},
    ):
        # --- argument handling --- #
        if isinstance(quantities, str) or isinstance(quantities, Callable):
            quantities = [quantities]
        theta_linear = isinstance(theta, int)
        zeta_linear = isinstance(zeta, int)
        # ds = Evaluations(self, rho=rho, theta=theta, zeta=zeta)
        ds = NotImplemented
        # --- plotting --- #
        fig, axs = plt.subplots(
            len(quantities), ds.rho.size, tight_layout=True, **figkwargs
        )
        axs = axs.reshape((len(quantities), ds.rho.size))
        for q, quantity in enumerate(quantities):
            if isinstance(quantity, str):
                ds.compute(quantity)
                values = ds[quantity]
            else:
                values = quantity(ds)
                match values:
                    case (xr.DataArray() as v, str() as s):
                        v.attrs["symbol"] = s
                        values = v
            for r, rhoi in enumerate(ds.rho):
                v = values.isel(rho=r).transpose("zeta", "theta").values
                if share_colorbar:
                    mesh = axs[q, r].contourf(
                        ds.theta,
                        ds.zeta,
                        v,
                        vmin=values.min(),
                        vmax=values.max(),
                        **ckwargs,
                    )
                    if r == ds.rho.size - 1:
                        divider = make_axes_locatable(axs[q, r])
                        cax = divider.append_axes("right", size="5%", pad=0.05)
                        fig.colorbar(mesh, cax=cax)
                else:
                    mesh = axs[q, r].contourf(ds.theta, ds.zeta, v, **ckwargs)
                    fig.colorbar(mesh, ax=axs[q, r])
                if theta_linear:
                    axs[q, r].set_xticks(np.linspace(0, 2 * np.pi, n_xtick))
                if zeta_linear:
                    axs[q, r].set_yticks(np.linspace(0, 2 * np.pi / self.nfp, n_ytick))
                axs[q, r].set(
                    xlim=(ds.theta[0], ds.theta[-1]),
                    ylim=(ds.zeta[0], ds.zeta[-1]),
                )
                if q != len(quantities) - 1:
                    axs[q, r].set_xticklabels(
                        np.full_like(axs[q, r].get_xticks(), "", dtype="U0")
                    )
                if r != 0:
                    axs[q, r].set_yticklabels(
                        np.full_like(axs[q, r].get_yticks(), "", dtype="U0")
                    )
            symbol = f"${values.attrs['symbol']}$" if "symbol" in values.attrs else "?"
            axs[q, 0].set_ylabel(f"{symbol}\n$\\zeta$")
            if zeta_linear:
                axs[q, 0].set_yticklabels(
                    [
                        rf"$\frac{{{2*i} \pi}}{{{self.nfp * (n_ytick - 1)}}}$"
                        for i in range(n_ytick)
                    ]
                )
        for r, rhoi in enumerate(ds.rho):
            axs[0, r].set(
                title=rf"$\rho = {rhoi:.2f}$",
            )
            if theta_linear:
                axs[-1, r].set_xticklabels(
                    [rf"$\frac{{{2*i} \pi}}{{{n_xtick - 1}}}$" for i in range(n_xtick)]
                )
            axs[-1, r].set_xlabel(r"$\theta$")
        fig.suptitle(self.statefile.name)
        return fig, axs

    def plot_poloidal(
        self,
        quantities: (
            str
            | Callable[[xr.Dataset], xr.DataArray | tuple[xr.DataArray, str]]
            | Iterable[
                str | Callable[[xr.Dataset], xr.DataArray | tuple[xr.DataArray, str]]
            ]
        ),
        rho: int | tuple[float, float, int] | np.ndarray = 81,
        theta: int | tuple[float, float, int] | np.ndarray = 101,
        zeta: int | tuple[float, float, int] | np.ndarray = 3,
        share_colorbar: bool = True,
        figkwargs: Mapping = {},
        ckwargs: Mapping = {},
    ):
        # --- argument handling --- #
        if isinstance(quantities, str) or isinstance(quantities, Callable):
            quantities = [quantities]
        # ds = Evaluations(self, rho=rho, theta=theta, zeta=zeta)
        ds = NotImplemented
        ds.compute("X1", "X2")
        # --- plotting --- #
        fig, axs = plt.subplots(
            len(quantities),
            ds.zeta.size,
            tight_layout=True,
            sharex=True,
            sharey=True,
            **figkwargs,
        )
        axs = np.asarray(axs).reshape((len(quantities), ds.zeta.size))
        for q, quantity in enumerate(quantities):
            if isinstance(quantity, str):
                ds.compute(quantity)
                values = ds[quantity]
            else:
                values = quantity(ds)
                match values:
                    case (xr.DataArray() as v, str() as s):
                        v.attrs["symbol"] = s
                        values = v
            for z, zetai in enumerate(ds.zeta):
                dsz = ds.isel(zeta=z)
                v = values.isel(zeta=z)
                if share_colorbar:
                    mesh = axs[q, z].contourf(
                        dsz.X1,
                        dsz.X2,
                        v,
                        vmin=values.min(),
                        vmax=values.max(),
                        **ckwargs,
                    )
                    if z == ds.zeta.size - 1:
                        divider = make_axes_locatable(axs[q, z])
                        cax = divider.append_axes("right", size="5%", pad=0.05)
                        fig.colorbar(mesh, cax=cax)
                else:
                    mesh = axs[q, z].contourf(dsz.X1, dsz.X2, v, **ckwargs)
                    fig.colorbar(mesh, ax=axs[q, z])
                axs[q, z].set(
                    aspect="equal",
                )
            symbol = f"${values.attrs['symbol']}$" if "symbol" in values.attrs else "?"
            axs[q, 0].set_ylabel(f"{symbol}\n$X^2$")
        for z, zetai in enumerate(ds.zeta):
            axs[0, z].set(
                title=rf"$\zeta = {zetai:.2f}$",
            )
            axs[-1, z].set_xlabel(r"$X^1$")
        fig.suptitle(self.statefile.name)
        return fig, axs
