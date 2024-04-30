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


def test_post(pygvec, ellipstell):
    paramfile = "parameter.ini"
    statefile = "ELLIPSTELL_LOWRES_State_0000_00000000.dat"

    with pygvec.post.State(paramfile, statefile) as state:
        assert isinstance(state, pygvec.post.State)


@pytest.mark.xfail()
def test_post_twice(pygvec, ellipstell):
    paramfile = "parameter.ini"
    statefile = "ELLIPSTELL_LOWRES_State_0000_00000000.dat"

    with pygvec.post.State(paramfile, statefile) as state:
        pass

    with pygvec.post.State(paramfile, statefile) as state:
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


def test_compute_base(pygvec, ellipstell_state):
    rho = np.linspace(0, 1, 6)
    theta = np.linspace(0, 2 * np.pi, 32, endpoint=False)
    zeta = np.linspace(0, 2 * np.pi, 10, endpoint=False)
    ds = pygvec.post.Evaluations(coords={"rho": rho, "theta": theta, "zeta": zeta})

    ds.compute("X1", ellipstell_state)
    ds.compute("X2", ellipstell_state)
    ds.compute("LA", ellipstell_state)

    assert ds.X1.shape == (6, 32, 10)
    assert "dX1_drho" in ds
    assert "dX2_dtheta" in ds
    assert "dLA_dzeta" in ds
    assert not np.any(np.isnan(ds.X1))


def test_evaluate_hmap(ellipstell_state):
    rho = np.linspace(0, 1, 6)
    theta = np.linspace(0, 2 * np.pi, 32, endpoint=False)
    zeta = np.linspace(0, 2 * np.pi, 10, endpoint=False)
    R, T, Z = np.meshgrid(rho, theta, zeta, indexing="ij")
    # X1, X2, dX1_drho, dX2_drho, dX1_dtheta, dX2_dtheta, dX1_dzeta, dX2_dzeta
    X1 = ellipstell_state.evaluate_base_tens_all("X1", rho, theta, zeta)
    X2 = ellipstell_state.evaluate_base_tens_all("X2", rho, theta, zeta)
    inputs = sum([[x1, x2] for x1, x2 in zip(X1, X2)], [])
    inputs = inputs[:2] + [T] + inputs[2:]
    inputs = [i.flatten() for i in inputs]

    outputs = ellipstell_state.evaluate_hmap(*inputs)
    assert len(outputs) == 4
    assert all(output.shape == (3, 6 * 32 * 10) for output in outputs)


def test_compute_hmap(pygvec, ellipstell_state):
    rho = np.linspace(0, 1, 6)
    theta = np.linspace(0, 2 * np.pi, 32, endpoint=False)
    zeta = np.linspace(0, 2 * np.pi, 10, endpoint=False)
    ds = pygvec.post.Evaluations(coords={"rho": rho, "theta": theta, "zeta": zeta})

    ds.compute("pos", ellipstell_state)

    assert ds.pos.shape == (3, 6, 32, 10)
    assert "vector" in ds.coords
    assert "e_rho" in ds
    assert "e_theta" in ds
    assert "e_zeta" in ds
    assert "dX1_drho" in ds
    assert not np.any(np.isnan(ds.pos))


def test_evaluate_profile(ellipstell_state):
    rho = np.linspace(0, 1, 6)
    ds = xr.Dataset(coords={"rho": rho})

    iota = ellipstell_state.evaluate_profile("iota", rho)
    assert iota.shape == rho.shape
    dp_drho = ellipstell_state.evaluate_profile("p_prime", rho)


@pytest.mark.parametrize(
    "quantity",
    [
        "iota",
        "diota_drho",
        "p",
        "dp_drho",
        "Phi",
        "dPhi_drho",
        "d2Phi_drho2",
        "chi",
        "dchi_drho",
        "Phi_n",
        "dPhi_n_drho",
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
