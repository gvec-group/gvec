import numpy as np
import xarray as xr
from numpy.typing import ArrayLike
from gvec.lib import modgvec_biotsavart as _BS
from warnings import warn
from collections.abc import Iterable
from pathlib import Path
from scipy.integrate import solve_ivp
from tqdm import tqdm
import os
from gvec.util import logging_setup
from gvec.core.state import State, CoordinateSpec
from gvec.core.compute import volume_integral
from logging import getLogger
from typing import Literal

try:
    from joblib import Parallel, delayed, parallel_config

    has_joblib = True

    # wrapper to enable progressbar for parallel fieldline tracing
    class ProgressParallel(Parallel):
        def __init__(self, use_tqdm=True, total=None, *args, **kwargs):
            self._use_tqdm = use_tqdm
            self._total = total
            super().__init__(*args, **kwargs)

        def __call__(self, *args, **kwargs):
            with tqdm(disable=not self._use_tqdm, total=self._total) as self._pbar:
                return Parallel.__call__(self, *args, **kwargs)

        def print_progress(self):
            if self._total is None:
                self._pbar.total = self.n_dispatched_tasks
            self._pbar.n = self.n_completed_tasks
            self._pbar.refresh()

except ImportError:
    has_joblib = False

mu_0 = 4 * np.pi * 0.99999999987 * 1e-7


def _stack_for_BiotSavart(ds: xr.Dataset | xr.DataArray):
    if "points" in ds.dims:
        return ds
    else:
        stack_dims = [dim for dim in ds.dims if dim != "xyz"]
        ds = ds.stack(points=stack_dims)
        return ds


class Coil:
    def __init__(self, coil_points: ArrayLike, coil_current: float):
        """Filament coil discretized via straight line segments.

        Parameters
        ----------
        coil_points : ArrayLike
            Discretization points of the coil in Cartesian coordinates. Expected shape is (3,n_points).
        coil_current : float
            Current flowing in the coil in [A].
        """
        if not isinstance(coil_points, xr.DataArray):
            coil_points = np.atleast_2d(coil_points)
            self.coil_points = xr.DataArray(
                coil_points, dims=("xyz", "points"), coords={"xyz": ("xyz", ["x", "y", "z"])}
            )
        else:
            self.coil_points = coil_points
        self.coil_current = coil_current
        self.n_points = self.coil_points.shape[1]
        self.e_vec = self.coil_points[:, 1:] - self.coil_points[:, :-1]
        self.L = self._mod_vector(self.e_vec)
        self.e_hat_vec = self.e_vec / self.L
        self.prefactor = self.coil_current * mu_0 / (4 * np.pi)
        delta = abs(self.coil_points.isel(points=0) - self.coil_points.isel(points=-1))
        if any(delta > 1e-12):
            raise ValueError(
                "First and last point of Coil is not the same! Please close the Coil."
            )

    def __call__(self, pos: ArrayLike):
        return self.eval_B(pos)

    def __repr__(self):
        return f"Coil(points={self.n_points}, current={self.coil_current:.2e} A)"

    def _mod_vector(self, xyz):
        return np.sqrt(xr.dot(xyz, xyz, dim="xyz"))

    def eval_B(self, pos: ArrayLike):
        """Evaluate the magnetic field due to currents flowing in the coil at given positions.

        Parameters
        ----------
        pos : ArrayLike
            Cartesian coordinates of the evaluation points with shape (3,n_positions).

        Returns
        -------
        xr.Dataset
            Dataset containing the magnetic field at the evaluation positions.
        """
        if not isinstance(pos, xr.DataArray):
            pos = np.asanyarray(pos)
            if pos.ndim == 1:
                pos = pos.reshape(-1, 1)
            pos = xr.DataArray(
                pos,
                dims=("xyz", "points"),
                coords={"xyz": ("xyz", ["x", "y", "z"])},
                attrs=dict(long_name="position vector", symbol=r"\mathbf{x}"),
            )
        n_positions = pos.shape[1]

        B = np.asfortranarray(np.zeros(pos.shape))

        _BS.biotsavart(
            n_positions=n_positions,
            xyz=pos,
            n_points=self.n_points,
            coil_points=self.coil_points,
            prefactor=self.prefactor,
            b=B,
        )

        ds = xr.Dataset(
            data_vars=dict(B=(("xyz", "points"), B)),
            coords={"xyz": ("xyz", ["x", "y", "z"])},
        )
        ds["pos"] = pos

        ds.B.attrs = dict(long_name="magnetic coil field", symbol=r"\mathbf{B}_C")
        return ds

    def eval_mod_B(self, pos: ArrayLike):
        """Evaluate the modulus of the magnetic field due to currents flowing in the coil at given positions.

        Parameters
        ----------
        pos : ArrayLike
            Cartesian coordinates of the evaluation points with shape (3,n_positions).

        Returns
        -------
        ds : xr.Dataset
            Dataset containing the magnetic field and its modulus at the evaluation positions.
        """
        ds = self.eval_B(pos)
        ds["mod_B"] = self._mod_vector(ds.B)
        return ds

    def get_as_dataobject(self):
        """Return the coil as an xarray data-set.

        Returns
        -------
        Dataset
            Dataset containing the coil-points and the current.
        """
        ds = xr.Dataset()
        ds["pos"] = self.coil_points
        ds["current"] = self.coil_current
        ds.pos.attrs = {
            "long_name": "Positions of the coil discretization points in Cartesian coordinates."
        }
        ds.current.attrs = {"long_name": "Coil current in [A]."}
        return ds

    def save(self, path: str | Path, **kwargs):
        """Safe the coils discretization points and current in a netcdf file.

        Parameters
        ----------
        path : str | Path
            Defines where to save te file.
        """
        ds = self.get_as_dataobject()
        ds.to_netcdf(path, **kwargs)

    @classmethod
    def load(cls, path: str | Path, **kwargs):
        ds = xr.open_dataset(path, **kwargs)
        return cls(ds.pos, ds.current)

    def plot(self, ax=None, show=True, name: str = None, **kwargs):
        """Visualize the coil in 3D via plotly.

        Parameters
        ----------
        ax : plotly.graph_objects.Figure, optional
            Existing plotly figure to add the coil to, by default None
        show : bool, optional
            Whether to show the plot, by default True
        name : str, optional
            Label for the coil.

        Returns
        -------
        plotly.graph_objects.Figure
            3D plotly figure.
        """
        try:
            from plotly import graph_objects as go
        except ImportError:
            warn("plotly not installed, cannot plot coil.")
        if ax is None:
            ax = go.Figure()
        if name is None:
            name = self.__repr__()
        ax.add_trace(
            go.Scatter3d(
                x=self.coil_points[0, :],
                y=self.coil_points[1, :],
                z=self.coil_points[2, :],
                mode="lines",
                name=name,
                **kwargs,
            )
        )
        if show:
            ax.show()
        return ax


class CoilSet(Coil):
    def __init__(self, coils: list[Coil], coil_names: list[str] = None):
        """Initialize a coil set from a list of coils.

        Parameters
        ----------
        coils : list[Coil]
            Coils to be added to the set.
        coil_names : list[str], optional
            List of names for the coils, by default None
        """
        if coil_names is None:
            coil_names = [f"coil_{i}" for i in range(len(coils))]
        self.coils = {coil_name: coil for coil_name, coil in zip(coil_names, coils)}

    def __repr__(self):
        return_str = ""
        for coil_name, coil in self.coils.items():
            return_str += f"{coil_name}: {coil}\n"
        return return_str

    def __getitem__(self, key: str):
        return self.coils[key]

    def get_as_dataobject(self):
        dt = xr.DataTree()
        for coil in self.coils:
            dt[coil] = self.coils[coil].get_as_dataobject()
        return dt

    @classmethod
    def load(cls, path: str | Path, **kwargs):
        dt = xr.open_datatree(path, **kwargs)
        coil_list = [Coil(dt[name].pos, dt[name].current) for name in dt]
        coil_names = [name for name in dt]
        return cls(coil_list, coil_names)

    def eval_B(self, pos: ArrayLike):
        """Evaluate the magnetic field due to currents flowing in the coil at given positions.

        Parameters
        ----------
        pos : ArrayLike
            Evaluation positions in xyz with shape (3,n_positions).
            If a DataArray is provided, "xyz" is an expected dimensions.

        Returns
        -------
        xr.Dataset
            Magnetic field at the evaluation positions.
        """
        stacked = False
        if not isinstance(pos, xr.DataArray):
            pos = np.asanyarray(pos)
            if pos.ndim == 1:
                pos = pos.reshape(-1, 1)
            pos = xr.DataArray(
                pos,
                dims=("xyz", "points"),
                coords={"xyz": ("xyz", ["x", "y", "z"])},
                attrs=dict(long_name="position vector", symbol=r"\mathbf{x}"),
            )
        else:
            pos = _stack_for_BiotSavart(pos)
            stacked = True
        n_positions = pos.shape[1]
        B_aux = np.asfortranarray(np.zeros(pos.shape))
        ds = xr.Dataset(
            data_vars=dict(B=(("xyz", "points"), np.zeros(pos.shape))),
            coords={"xyz": ("xyz", ["x", "y", "z"])},
        )

        ds["pos"] = pos
        for coil_name in self.coils:
            coil = self.coils[coil_name]
            _BS.biotsavart(
                n_positions=n_positions,
                xyz=pos,
                n_points=coil.n_points,
                coil_points=coil.coil_points,
                prefactor=coil.prefactor,
                b=B_aux,
            )
            ds.B[:, :] += B_aux
            B_aux *= 0.0
        ds.B.attrs = dict(long_name="magnetic coil field", symbol=r"\mathbf{B}_C")
        if stacked:
            ds = ds.unstack("points")
            for dimname, actual_dimname in zip(["rho", "theta", "zeta"], ["rad", "pol", "tor"]):
                if dimname in ds.dims:
                    ds = ds.rename_dims({dimname: actual_dimname})

        return ds

    def _eval_B_direct(self, pos: ArrayLike):
        """Evaluate the magnetic field due to currents flowing in the coil at given positions.
        This routine skips the xarray wrapping of the output.
        Useful when the B-field evaluation has to be called repeatedly for few positions (e.g. field line tracing).

        Parameters
        ----------
        pos : ArrayLike
            Evaluation positions in xyz with shape (3,n_positions).

        Returns
        -------
        np.ndarray
            Magnetic field at the evaluation positions in the shape (3,n_positions).
        """
        pos = np.asanyarray(pos)
        if pos.ndim == 1:
            pos = pos.reshape(-1, 1)
        n_positions = pos.shape[1]
        B = np.zeros(pos.shape)
        B_aux = np.asfortranarray(np.zeros(pos.shape))
        for coil_name in self.coils:
            coil = self.coils[coil_name]
            _BS.biotsavart(
                n_positions=n_positions,
                xyz=pos,
                n_points=coil.n_points,
                coil_points=coil.coil_points,
                prefactor=coil.prefactor,
                b=B_aux,
            )
            B += B_aux
            B_aux *= 0.0
        return B

    def plot(self, ax=None, show=True, **kwargs):
        """Visualize the coil set in 3D via plotly.

        Parameters
        ----------
        ax : plotly.graph_objects.Figure, optional
            Existing plotly figure to add the coils to, by default None
        show : bool, optional
            Whether to show the plot, by default True

        Returns
        -------
        plotly.graph_objects.Figure
            3D plotly figure.
        """
        for coil_name in self.coils:
            ax = self.coils[coil_name].plot(ax=ax, show=False, name=coil_name, **kwargs)
        if show:
            ax.show()
        return ax


def check_plane_factory(surf_normal: np.ndarray, surf_point: np.ndarray):
    """Function factory for checking intersection with a plane defined by a normal vector.
    Used for the poincare plots in fieldline tracing.

    Parameters
    ----------
    surf_normal : np.ndarray
        Normal vector of the plane.
    """

    def check_plane(t, R):
        R = np.array(R)
        dot = np.sum((R - surf_point) * surf_normal)
        return dot

    return check_plane


def intersection_planes_from_state(
    state: State, zetas: CoordinateSpec = 3, min_box_size: float = 1e-3
):
    """Create IntersectionPlane objects from a State object at specified zeta positions.

    Parameters
    ----------
    state : State
        GVEC state object to construct the planes from.
    zetas : CoordinateSpec, optional
        Zeta positions of the desired planes.
        If an integer is given linearly spaced zetas will be used. By default 3
    min_box_size : float, optional
        Minimum length of the bounding box, by default 1e-3.

    Returns
    -------
    List[IntersectionPlane]
        List of IntersectionPlane objects.
    """

    boundary = state.evaluate("pos", zeta=zetas, theta=128, rho=1.0).pos

    planes = []
    for zeta in boundary.zeta:
        surf_point, e_q1, e_q2, e_q3 = state.evaluate_hmap_only(X1=0, X2=0, zeta=zeta)
        e_q3 = e_q3[:, 0]
        surf_normal_length = np.sqrt(np.dot(e_q3, e_q3))
        if surf_normal_length <= 1e-12:  # RZ case
            major_radius = state.evaluate("major_radius", zeta=zeta).major_radius.data
            e_q1, e_q2, e_q3 = state.evaluate_hmap_only(X1=major_radius, X2=0, zeta=zeta)[1:]
            e_q3 = e_q3[:, 0]
            surf_normal_length = np.sqrt(np.dot(e_q3, e_q3))
        surf_normal = e_q3 / surf_normal_length

        surf_point = surf_point[:, 0]
        e_q1 = e_q1[:, 0]
        e_q2 = e_q2[:, 0]

        boundary_x = boundary.sel(rho=1.0, xyz="x", zeta=zeta)
        boundary_y = boundary.sel(rho=1.0, xyz="y", zeta=zeta)
        boundary_z = boundary.sel(rho=1.0, xyz="z", zeta=zeta)

        xlim = np.array([boundary_x.min().data, boundary_x.max().data])
        ylim = np.array([boundary_y.min().data, boundary_y.max().data])
        zlim = np.array([boundary_z.min().data, boundary_z.max().data])
        dx = xlim[1] - xlim[0]
        dy = ylim[1] - ylim[0]
        dz = zlim[1] - zlim[0]

        if abs(dx) < min_box_size:
            xlim[0] -= min_box_size / 2
            xlim[1] += min_box_size / 2

        if abs(dy) < min_box_size:
            ylim[0] -= min_box_size / 2
            ylim[1] += min_box_size / 2

        if abs(dz) < min_box_size:
            zlim[0] -= min_box_size / 2
            zlim[1] += min_box_size / 2

        planes.append(
            IntersectionPlane(
                surf_normal=surf_normal,
                surf_point=surf_point,
                xlim=xlim,
                ylim=ylim,
                zlim=zlim,
                e_q1=e_q1,
                e_q2=e_q2,
            )
        )

    return planes


class IntersectionPlane:
    def __init__(
        self,
        surf_normal: np.ndarray,
        surf_point: np.ndarray,
        xlim: ArrayLike = np.inf,
        ylim: ArrayLike = np.inf,
        zlim: ArrayLike = np.inf,
        e_q1: ArrayLike = None,
        e_q2: ArrayLike = None,
    ):
        """Define a plane for Poincaré plots with a bounding box.

        Parameters
        ----------
        surf_normal : np.ndarray
            Surface normal of the plane.
        surf_point : np.ndarray
            Point on the plane. If e_q1 and e_q2 are specified, this is assumed to be the origin.
        xlim : ArrayLike, optional
            X limits of the bounding box, by default np.inf
        ylim : ArrayLike, optional
            Y limits of the bounding box, by default np.inf
        zlim : ArrayLike, optional
            Z limits of the bounding box, by default np.inf
        """

        self.surf_normal = surf_normal
        self.surf_point = surf_point
        self.xlim = xlim
        self.ylim = ylim
        self.zlim = zlim
        self.e_q1 = e_q1
        self.e_q2 = e_q2

    def transform_events_to_X1_X2(self, events):
        if self.e_q1 is None or self.e_q2 is None:
            warn("Quantities not set, cannot transform events to X1 and X2.")
            return
        X1 = np.dot((events - self.surf_point), self.e_q1)
        X2 = np.dot((events - self.surf_point), self.e_q2)
        return (X1, X2)

    def __call__(self, t, R):
        R = np.array(R)
        dot = np.sum((R - self.surf_point) * self.surf_normal)
        return dot


def trace_fieldlines(
    starts: np.ndarray | xr.Dataset,
    coils: CoilSet,
    t: float,
    surf_normals: list[np.ndarray] = None,
    surf_points: list[np.ndarray] = None,
    n_jobs: int = 1,
    return_solves: bool = False,
    verbosity: Literal["WARNING", "INFO", "DEBUG"] = "INFO",
    **kwargs,
):
    """
    Trace magnetic field lines of a coil field.

    Parameters
    ----------
    starts : np.ndarray | xr.Dataset
        Initial positions of the magnetic field lines in Carthesian coordinates. Expected shape (3,n_fieldlines).
    coils : CoilSet
        Coils set used for evaluating the magnetic field.
    t : float
        Time for which to trace the field lines so that t_span = [0,t].
    surf_normals : list[np.ndarray], optional
        List of normal vectors of the planes for which intersections of the field lines should be checked.
        This can then be used for Poincare plots. The default is None
    surf_points : list[np.ndarray], optional
        List of points on the planes for which intersections of the field lines should be checked.
        If None, the origin is used. However, this does not work for the G-frame.
    n_jobs : int, optional
        Number of jobs for parallelization over fieldlines, by default 1.
    return_solves : bool, optional
        Whether to return a list of solve_ivp objects or an xarray datatree.
    verbosity: str
        Level of the logger.

    """
    logging_setup()
    logger = getLogger("gvec_fieldlines")
    logger.setLevel(verbosity)
    if n_jobs > 1 and not has_joblib:
        logger.warning(
            "n_jobs > 1 but joblib not installed, parallelization over fieldlines is not possible. Falling back to n_jobs=1."
        )
        n_jobs = 1

    if isinstance(starts, xr.DataArray):
        starts = _stack_for_BiotSavart(starts)
        starts = starts.transpose("xyz", "points")

    def _push_cart(t, R):
        B = coils._eval_B_direct(R)
        B = B / np.sqrt(np.sum(B**2, axis=0))
        return B[:, 0]

    check_planes = []

    if surf_normals is not None:
        if surf_points is None:
            surf_points = [np.zeros(3) for _ in range(len(surf_normals))]
        for surf_normal, surf_point in zip(surf_normals, surf_points):
            check_planes.append(check_plane_factory(surf_normal, surf_point))

    if "events" not in kwargs:
        kwargs["events"] = check_planes
    else:
        kwargs["events"] = check_planes + kwargs["events"]

    solves = []
    logger.info("Tracing fieldlines using n_jobs=%d", n_jobs)
    if n_jobs > 1:  # embarrassingly parallel over fieldlines
        os.environ["OMP_NUM_THREADS"] = "1"
        with parallel_config(n_jobs=n_jobs):
            solves = ProgressParallel(total=starts.shape[1])(
                delayed(solve_ivp)(
                    fun=_push_cart,
                    y0=starts[:, i],
                    t_span=[0, t],
                    **kwargs,
                )
                for i in range(starts.shape[1])
            )
    else:  # serial over fieldlines but still OMP parallelization for magnetic field evaluation
        for i in tqdm(range(starts.shape[1])):
            solve = solve_ivp(
                fun=_push_cart,
                y0=starts[:, i],
                t_span=[0, t],
                **kwargs,
            )
            solves.append(solve)
    if return_solves:
        return solves
    else:
        dt = xr.DataTree()
        for fieldline, output in enumerate(solves):
            coords = {
                "t": output.t,
                "xyz": ["x", "y", "z"],
            }
            for i, t_event in enumerate(output.t_events):
                coords[f"intersect_time_plane_{i}"] = t_event
            ds = xr.Dataset(
                {
                    "pos": (("xyz", "t"), output.y),
                },
                coords=coords,
            )
            try:
                for i, event_aux in enumerate(kwargs["events"]):
                    ds[f"event_{i}"] = (
                        (f"time_event_{i}", "xyz"),
                        output.y_events[i],
                    )
                    if type(event_aux) is IntersectionPlane:
                        ds_event = ds[f"event_{i}"]
                        mask_x = (event_aux.xlim[0] <= ds_event.sel(xyz="x")) & (
                            ds_event.sel(xyz="x") <= event_aux.xlim[1]
                        )
                        mask_y = (event_aux.ylim[0] <= ds_event.sel(xyz="y")) & (
                            ds_event.sel(xyz="y") <= event_aux.ylim[1]
                        )
                        mask_z = (event_aux.zlim[0] <= ds_event.sel(xyz="z")) & (
                            ds_event.sel(xyz="z") <= event_aux.zlim[1]
                        )
                        ds[f"event_{i}"] = ds_event.where(mask_x & mask_y & mask_z)
                        if hasattr(event_aux, "e_q1") and hasattr(event_aux, "e_q2"):
                            X1, X2 = event_aux.transform_events_to_X1_X2(ds[f"event_{i}"])
                            ds[f"event_{i}_X1"] = ((f"time_event_{i}"), X1)
                            ds[f"event_{i}_X2"] = ((f"time_event_{i}"), X2)
            except ValueError:
                logger.debug(
                    f"Error when extracting event-data for event {i} of fieldline {fieldline}. Maybe no intersection occurred."
                )
                pass
            dt[f"fieldline_{fieldline}"] = ds
        return dt


def get_phi_edge_from_coils(state: State, coil_set: CoilSet | Coil):
    """Calculate the total toroidal flux from the geometry of a state object and the magnetic field of a coil-set.

    Parameters
    ----------
    state : State
        State object used for evaluating the geometry.
    coil_set : CoilSet | Coil
        Coil-set used for evaluating the magnetic field.

    Returns
    -------
    float
        Averaged total toroidal flux.
    """
    ev = state.evaluate("pos", "B", "grad_zeta")
    B_coils = coil_set.eval_mod_B(ev.pos)

    Bn = xr.dot(B_coils.B, ev.grad_zeta, dim="xyz")
    BnJac = Bn * ev.Jac
    Psi_vol = volume_integral(BnJac) / (np.pi * 2)
    return Psi_vol.item()
