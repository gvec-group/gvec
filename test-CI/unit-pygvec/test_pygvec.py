import pytest
import numpy as np
import xarray as xr
import os
from pathlib import Path

import helpers


# === Fixtures === #


@pytest.fixture(scope="session")
def pygvec():
    import pygvec

    return pygvec


@pytest.fixture()
def ellipstell(tmpdir):
    """prepare the ellipstell parameters"""
    refdir = Path("test-CI/run/example/ellipstell_lowres").absolute()
    paramfile = "parameter.ini"
    statefile = "ELLIPSTELL_LOWRES_State_0000_00000000.dat"
    with helpers.chdir(tmpdir):
        os.symlink(refdir / paramfile, paramfile)
        os.symlink(refdir / statefile, statefile)
        yield


@pytest.fixture()
def ellipstell_state(pygvec, ellipstell):
    paramfile = "parameter.ini"
    statefile = "ELLIPSTELL_LOWRES_State_0000_00000000.dat"
    with pygvec.post.State(paramfile, statefile) as state:
        yield state


@pytest.fixture()
def evals_r(pygvec, ellipstell_state):
    rho = np.linspace(0, 1, 6)
    ds = pygvec.post.Evaluations(state=ellipstell_state, coords={"rho": rho})
    return ds


@pytest.fixture()
def evals_rtz(pygvec, ellipstell_state):
    rho = np.linspace(0, 1, 6)
    theta = np.linspace(0, 2 * np.pi, 32, endpoint=False)
    zeta = np.linspace(0, 2 * np.pi, 10, endpoint=False)
    ds = pygvec.post.Evaluations(
        state=ellipstell_state, coords={"rho": rho, "theta": theta, "zeta": zeta}
    )
    return ds


# === Tests === #


def test_version(pygvec):
    import pkg_resources

    assert isinstance(pygvec.__version__, str)
    assert pygvec.__version_tuple__ >= (0, 2, 1)
    assert pkg_resources.get_distribution("pygvec").version == pygvec.__version__


def test_state(pygvec, ellipstell):
    paramfile = "parameter.ini"
    statefile = "ELLIPSTELL_LOWRES_State_0000_00000000.dat"

    with pygvec.post.State(paramfile, statefile) as state:
        assert isinstance(state, pygvec.post.State)


def test_state_args(pygvec, ellipstell):
    paramfile = "parameter.ini"
    statefile = "ELLIPSTELL_LOWRES_State_0000_00000000.dat"

    with pytest.raises(FileNotFoundError):
        state = pygvec.post.State(paramfile, "nonexistent.dat")

    with pytest.raises(FileNotFoundError):
        state = pygvec.post.State("nonexistent.ini", statefile)


def test_state_explicit(pygvec, ellipstell):
    paramfile = "parameter.ini"
    statefile = "ELLIPSTELL_LOWRES_State_0000_00000000.dat"

    state = pygvec.post.State(paramfile, statefile)
    assert isinstance(state, pygvec.post.State)
    assert state.initialized
    state.finalize()
    assert not state.initialized

    with pytest.raises(RuntimeError):
        state.evaluate_base_tens("X1", None, [0.5], [0.5], [0.5])


def test_state_twice(pygvec, ellipstell):
    paramfile = "parameter.ini"
    statefile = "ELLIPSTELL_LOWRES_State_0000_00000000.dat"

    # double context
    with pygvec.post.State(paramfile, statefile) as state:
        pass
    with pygvec.post.State(paramfile, statefile) as state:
        pass

    # double explicit
    state = pygvec.post.State(paramfile, statefile)
    state.finalize()
    state = pygvec.post.State(paramfile, statefile)
    state.finalize()

    # simultaneous context (not implemented yet)
    with pytest.raises(NotImplementedError):
        with pygvec.post.State(paramfile, statefile) as state1:
            with pygvec.post.State(paramfile, statefile) as state2:
                pass


def test_evaluate_base_tens(ellipstell_state):
    rho = np.linspace(0, 1, 6)
    theta = np.linspace(0, 2 * np.pi, 32, endpoint=False)
    zeta = np.linspace(0, 2 * np.pi, 10, endpoint=False)

    # base evaluation
    X1 = ellipstell_state.evaluate_base_tens("X1", None, rho, theta, zeta)
    X2 = ellipstell_state.evaluate_base_tens("X2", None, rho, theta, zeta)
    LA = ellipstell_state.evaluate_base_tens("LA", None, rho, theta, zeta)
    dX2_drtz = ellipstell_state.evaluate_base_tens("X2", "rtz", rho, theta, zeta)

    assert X1.shape == (6, 32, 10)
    assert not np.any(np.isnan(X1))
    assert not np.any(np.isnan(X2))
    assert not np.any(np.isnan(LA))
    # magnetic axis collapses to a line
    assert np.allclose(np.std(X1[0, ...], axis=0), 0.0)
    assert np.allclose(np.std(X2[0, ...], axis=0), 0.0)


def test_evaluate_base_tens_vectorize(ellipstell_state):
    """Test if the vectorization works properly (correct order of indices)"""
    X1_222 = ellipstell_state.evaluate_base_tens(
        "X1", None, [0.5, 0.6], [0, 0.1], [0, 0.1]
    )
    X1_122 = ellipstell_state.evaluate_base_tens("X1", None, [0.5], [0, 0.1], [0, 0.1])
    X1_212 = ellipstell_state.evaluate_base_tens(
        "X1", None, [0.5, 0.6], [0.1], [0, 0.1]
    )
    X1_221 = ellipstell_state.evaluate_base_tens(
        "X1", None, [0.5, 0.6], [0, 0.1], [0.1]
    )
    X1_112 = ellipstell_state.evaluate_base_tens("X1", None, [0.5], [0.1], [0, 0.1])
    assert np.allclose(X1_222[0, :, :], X1_122[0, :, :])
    assert np.allclose(X1_222[:, 1, :], X1_212[:, 0, :])
    assert np.allclose(X1_222[:, :, 1], X1_221[:, :, 0])
    assert np.allclose(X1_222[0, 1, :], X1_112[0, 0, :])


def test_evaluate_base_tens_bounds(ellipstell_state):
    rho = np.linspace(0, 1, 6)
    theta = np.linspace(0, 4 * np.pi, 8, endpoint=False)
    zeta = np.linspace(-2 * np.pi, 0, 10, endpoint=False)

    ellipstell_state.evaluate_base_tens("X1", None, rho, theta, zeta)
    with pytest.raises(ValueError):
        ellipstell_state.evaluate_base_tens("X3", None, rho, theta, zeta)
    rho = np.linspace(-1, 1, 6)
    with pytest.raises(ValueError):
        ellipstell_state.evaluate_base_tens("X1", None, rho, theta, zeta)
    rho = np.linspace(0, 2, 6)
    with pytest.raises(ValueError):
        ellipstell_state.evaluate_base_tens("X1", None, rho, theta, zeta)


def test_evaluate_base_tens_all(ellipstell_state):
    rho = np.linspace(0, 1, 6)
    theta = np.linspace(0, 2 * np.pi, 32, endpoint=False)
    zeta = np.linspace(0, 2 * np.pi, 10, endpoint=False)

    results = ellipstell_state.evaluate_base_tens_all("X1", rho, theta, zeta)
    assert len(results) == 4
    assert all(result.shape == (6, 32, 10) for result in results)
    with pytest.raises(ValueError):
        ellipstell_state.evaluate_base_tens_all("X3", rho, theta, zeta)


def test_compute_base(evals_rtz):
    ds = evals_rtz

    ds.compute("X1")
    ds.compute("X2")
    ds.compute("LA")

    assert ds.X1.shape == (6, 32, 10)
    assert "dX1_dr" in ds
    assert "dX2_dt" in ds
    assert "dLA_dz" in ds
    assert "dX1_drt" in ds
    assert not np.any(np.isnan(ds.X1))


def test_evaluate_hmap(ellipstell_state):
    rho = np.linspace(0, 1, 6)
    theta = np.linspace(0, 2 * np.pi, 32, endpoint=False)
    zeta = np.linspace(0, 2 * np.pi, 10, endpoint=False)
    R, T, Z = np.meshgrid(rho, theta, zeta, indexing="ij")
    # X1, X2, dX1_dr, dX2_dr, dX1_dt, dX2_dt, dX1_dz, dX2_dz
    X1 = ellipstell_state.evaluate_base_tens_all("X1", rho, theta, zeta)
    X2 = ellipstell_state.evaluate_base_tens_all("X2", rho, theta, zeta)
    inputs = sum([[x1, x2] for x1, x2 in zip(X1[:4], X2[:4])], [])
    inputs = inputs[:2] + [Z] + inputs[2:]
    inputs = [i.flatten() for i in inputs]

    outputs = ellipstell_state.evaluate_hmap(*inputs)
    assert len(outputs) == 4
    assert all(output.shape == (3, 6 * 32 * 10) for output in outputs)


def test_compute_hmap(evals_rtz):
    ds = evals_rtz
    ds.compute("pos")

    assert ds.pos.shape == (3, 6, 32, 10)
    assert "vector" in ds.coords
    assert "e_X1" in ds
    assert "e_X2" in ds
    assert "e_zeta3" in ds
    assert not np.any(np.isnan(ds.pos))

    ds.compute("e_rho")
    assert ds.e_rho.shape == (3, 6, 32, 10)


def test_compute_metric(evals_rtz):
    ds = evals_rtz

    ds.compute("g_tt")
    for ij in ("rr", "rt", "rz", "tt", "tz", "zz"):
        key = f"g_{ij}"
        assert key in ds
        assert set(ds[key].coords) == {"rho", "theta", "zeta"}
        for k in "rtz":
            assert f"dg_{ij}_d{k}" in ds

    ds.compute("e_rho")
    ds.compute("e_theta")
    ds.compute("e_zeta")
    idxs = {"r": "rho", "t": "theta", "z": "zeta"}
    for ij in "rr rt rz tt tz zz".split():
        key = f"g_{ij}"
        assert np.allclose(
            ds[f"g_{ij}"],
            xr.dot(ds[f"e_{idxs[ij[0]]}"], ds[f"e_{idxs[ij[1]]}"], dim="vector"),
        )


def test_evaluate_profile(ellipstell_state):
    rho = np.linspace(0, 1, 6)

    iota = ellipstell_state.evaluate_profile("iota", rho)
    assert iota.shape == rho.shape
    dp_dr = ellipstell_state.evaluate_profile("p_prime", rho)


@pytest.mark.parametrize(
    "quantity",
    [
        "iota",
        "diota_dr",
        "p",
        "dp_dr",
        "Phi",
        "dPhi_dr",
        "dPhi_drr",
        "chi",
        "dchi_dr",
        "Phi_n",
        "dPhi_n_dr",
    ],
)
def test_compute_profile(evals_r, quantity):
    ds = evals_r
    ds.compute(quantity)
    assert quantity in ds
    assert set(ds[quantity].coords) == {"rho"}
    assert ds[quantity].data.size == ds.rho.size
    assert not np.any(np.isnan(ds[quantity]))


@pytest.mark.parametrize(
    "quantity",
    [
        "Jac",
        "B",
    ],
)
def test_compute_derived(evals_rtz, quantity):
    ds = evals_rtz
    ds.compute(quantity)
    assert quantity in ds
    assert not np.any(np.isnan(ds[quantity]))


def test_compute_J(evals_rtz):
    ds = evals_rtz
    ds.compute("J")
    assert "J" in ds
    assert not np.any(np.isnan(ds["J"]))
