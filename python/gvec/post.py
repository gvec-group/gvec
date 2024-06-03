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
from typing import Mapping, Callable, Hashable, Iterable, Literal
import re
import inspect
import functools
import tempfile
import logging

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

    # === Constructor & Desctructor === #

    def __init__(
        self,
        parameterfile: str | Path,
        statefile: str | Path,
        redirect_stdout: bool = True,
    ):
        self.initialized = False
        self.parameterfile = None
        self.statefile = None

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
            f"<pygvec.State("
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
    def get_integration_points(self, quantity: str):
        if not isinstance(quantity, str):
            raise ValueError("Quantity must be a string.")
        elif quantity not in ["X1", "X2", "LA"]:
            raise ValueError(f"Unknown quantity: {quantity}")

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
        elif quantity not in ["X1", "X2", "LA"]:
            raise ValueError(f"Unknown quantity: {quantity}")
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
        elif quantity not in ["X1", "X2", "LA"]:
            raise ValueError(f"Unknown quantity: {quantity}")
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
        thetazeta = np.asfortranarray(thetazeta, dtype=np.float64, order="F")
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
    def evaluate_base_tens_all(
        self, quantity: str, rho: np.ndarray, theta: np.ndarray, zeta: np.ndarray
    ):
        if not isinstance(quantity, str):
            raise ValueError("Quantity must be a string.")
        elif quantity not in ["X1", "X2", "LA"]:
            raise ValueError(f"Unknown quantity: {quantity}")

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
        ds = Evaluations(self, rho=rho, theta=theta, zeta=zeta)
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
        ds = Evaluations(self, rho=rho, theta=theta, zeta=zeta)
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


class Evaluations(xr.Dataset):
    __slots__ = ["state"]
    _quantities = {}

    def __init__(
        self,
        state: State,
        rho: (
            int | Literal["int"] | tuple[float, float, int] | np.ndarray | None
        ) = "int",
        theta: (
            int | Literal["int"] | tuple[float, float, int] | np.ndarray | None
        ) = "int",
        zeta: (
            int | Literal["int"] | tuple[float, float, int] | np.ndarray | None
        ) = "int",
        **kwargs,
    ):
        self.state = state
        coords = kwargs["coords"] if "coords" in kwargs else {}
        data_vars = kwargs["data_vars"] if "data_vars" in kwargs else {}
        # --- get integration points --- #
        intp = [state.get_integration_points(q) for q in ["X1", "X2", "LA"]]
        match rho:
            case "int":
                if any(
                    [
                        not np.allclose(intp[0][j], intp[i][j])
                        for i in (1, 2)
                        for j in (0, 1)
                    ]
                ):
                    raise ValueError(
                        "Integration points for rho do not align for X1, X2 and LA."
                    )
                coords["rho"] = intp[0][0]
                data_vars["rho_weights"] = ("rho", intp[0][1])
            case int() as num:
                coords["rho"] = np.linspace(0, 1, num)
                coords["rho"][0] = 0.01
            case (start, stop):
                coords["rho"] = np.linspace(start, stop)
            case (start, stop, num):
                coords["rho"] = np.linspace(start, stop, num)
            case None:
                pass
            case _:
                coords["rho"] = rho
        match theta:
            case "int":
                if any(
                    [
                        not np.allclose(intp[0][j], intp[i][j])
                        for i in (1, 2)
                        for j in (2, 3)
                    ]
                ):
                    raise ValueError(
                        "Integration points for theta do not align for X1, X2 and LA."
                    )
                coords["theta"] = np.linspace(0, 2 * np.pi, intp[0][2], endpoint=False)
                data_vars["theta_weight"] = intp[0][3]
            case int() as num:
                coords["theta"] = np.linspace(0, 2 * np.pi, num)
            case (start, stop):
                coords["theta"] = np.linspace(start, stop)
            case (start, stop, num):
                coords["theta"] = np.linspace(start, stop, num)
            case None:
                pass
            case _:
                coords["theta"] = theta
        match zeta:
            case "int":
                if any(
                    [
                        not np.allclose(intp[2][j], intp[i][j])
                        for i in (1, 2)
                        for j in (4, 5)
                    ]
                ):
                    raise ValueError(
                        "Integration points for zeta do not align for X1, X2 and LA."
                    )
                coords["zeta"] = np.linspace(
                    0, 2 * np.pi / state.nfp, intp[0][4], endpoint=False
                )
                data_vars["zeta_weight"] = intp[0][5]
            case int() as num:
                coords["zeta"] = np.linspace(0, 2 * np.pi / state.nfp, num)
            case (start, stop):
                coords["zeta"] = np.linspace(start, stop)
            case (start, stop, num):
                coords["zeta"] = np.linspace(start, stop, num)
            case None:
                pass
            case _:
                coords["zeta"] = zeta

        # --- init Dataset --- #
        kwargs["coords"] = coords
        kwargs["data_vars"] = data_vars
        super().__init__(**kwargs)

        # --- set attributes --- #
        if "rho" in self:
            self.rho.attrs["long_name"] = "Logical radial coordinate"
            self.rho.attrs["symbol"] = r"\rho"
            self.rho.attrs["integration_points"] = rho == "int"
        if "theta" in self:
            self.theta.attrs["long_name"] = "Logical poloidal angle"
            self.theta.attrs["symbol"] = r"\theta"
            self.theta.attrs["integration_points"] = theta == "int"
        if "zeta" in self:
            self.zeta.attrs["long_name"] = "Logical toroidal angle"
            self.zeta.attrs["symbol"] = r"\zeta"
            self.zeta.attrs["integration_points"] = zeta == "int"

    @classmethod
    def register_compute_func(
        cls,
        quantities: None | str | Iterable[str] = None,
        requirements: Iterable[str] = (),
        integration: Iterable[str] = (),
    ):
        """Decorator to register equilibrium quantities.

        Will select the dependencies specified in `requirements` recursively and compute them if necessary.
        The required quantities will be selected based on their specified `dims`, allowing different computations for different dimensional layouts
         * tensorproduct of rho, theta, zeta
         * tensorproduct of rho, theta, zeta at integration points
        """

        def _register(
            func: (
                Callable[[Evaluations], Evaluations]
                | Callable[[Evaluations, State], Evaluations]
            )
        ):
            nonlocal quantities
            if quantities is None:
                quantities = [func.__name__]
            if isinstance(quantities, str):
                quantities = [quantities]
            func.quantities = quantities
            func.requirements = requirements
            func.integration = integration
            for q in quantities:
                if q in cls._quantities:
                    logging.warning(
                        f"A quantity `{q}` is already registered with {cls}."
                    )
                cls._quantities[q] = func
            return func

        return _register

    def compute(self, *quantities: Iterable[str]):
        """Compute the target equilibrium quantity.

        This method will compute required parameters recursively.
        """
        for quantity in quantities:
            # --- get the compute function --- #
            if quantity in self:
                continue  # already computed
            if quantity not in self._quantities:
                raise KeyError(f"The quantity `{quantity}` is not registered.")
            func = self._quantities[quantity]
            # --- handle requirements --- #
            for i in func.integration:
                if i not in self:
                    raise ValueError(f"Cannot compute `{quantity}` without `{i}`.")
                if (
                    "integration_points" not in self[i].attrs
                    or self[i].attrs["integration_points"] != True
                ):
                    raise ValueError(
                        f"Computation of `{quantity}` requires integration points for `{i}`."
                    )
            self.compute(*func.requirements)
            # --- compute the quantity --- #
            if "state" in inspect.signature(func).parameters:
                func(self, self.state)
            else:
                func(self)

    def __getitem__(self, key: Mapping | Hashable | Iterable[Hashable]):
        """Access variables or coordinates of this dataset as a `xarray.DataArray` or a subset of variables or a indexed dataset.

        Indexing with a list of names will return a new `Dataset` object.
        Automatically computes a registered quantity.
        """
        # --- compute quantities if possible --- #
        if not isinstance(key, str) and isinstance(key, Iterable):
            for k in key:
                if k in self._quantities:
                    self.compute(k)
        else:
            if key in self._quantities:
                self.compute(key)
        # --- pass on to xarray --- #
        return super().__getitem__(key)

    def radial_integral(self, quantity: str | xr.DataArray):
        """Compute the radial average of the given quantity."""
        if isinstance(quantity, str):
            self.compute(quantity)
            quantity = self[quantity]
        if not self.rho.attrs["integration_points"]:
            raise ValueError("Radial average requires integration points for rho.")
        return (quantity * self.rho_weights).sum("rho")

    def fluxsurface_integral(self, quantity: str | xr.DataArray):
        """Compute the flux surface average of the given quantity."""
        if isinstance(quantity, str):
            self.compute(quantity)
            quantity = self[quantity]
        if (
            not self.theta.attrs["integration_points"]
            and not self.zeta.attrs["integration_points"]
        ):
            raise ValueError(
                "Flux surface average requires integration points for theta and zeta."
            )
        return quantity.sum(("theta", "zeta")) * self.theta_weight * self.zeta_weight

    def fluxsurface_average(self, quantity: str | xr.DataArray):
        """Compute the flux surface average of the given quantity."""
        return self.fluxsurface_integral(quantity) / 4 / np.pi**2

    def volume_integral(self, quantity: str | xr.DataArray):
        """Compute the volume integral of the given quantity."""
        if isinstance(quantity, str):
            self.compute(quantity)
            quantity = self[quantity]
        if (
            not self.rho.attrs["integration_points"]
            and not self.theta.attrs["integration_points"]
            and not self.zeta.attrs["integration_points"]
        ):
            raise ValueError(
                "Volume integral requires integration points for rho, theta and zeta."
            )
        return (
            (quantity * self.rho_weights).sum(("rho", "theta", "zeta"))
            * self.theta_weight
            * self.zeta_weight
        )

    def volume_average(self, quantity: str | xr.DataArray):
        """Compute the volume average of the given quantity."""
        return self.volume_integral(quantity) / 4 / np.pi**2
