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


def test_evaluate_base(ellipstell_state):
    rho = np.linspace(0, 1, 6)
    theta = np.linspace(0, 2 * np.pi, 32, endpoint=False)
    zeta = np.linspace(0, 2 * np.pi, 10, endpoint=False)

    # base evaluation
    X1 = ellipstell_state.evaluate_base("X1", rho, theta, zeta)
    X2 = ellipstell_state.evaluate_base("X2", rho, theta, zeta)
    LA = ellipstell_state.evaluate_base("LA", rho, theta, zeta)
    D_rtzX2 = ellipstell_state.evaluate_base("D_rtz X2", rho, theta, zeta)

    assert X1.shape == (6, 32, 10)
    assert not np.any(np.isnan(X1))
    assert not np.any(np.isnan(X2))
    assert not np.any(np.isnan(LA))
    # magnetic axis collapses to a line
    assert np.allclose(np.std(X1[0, ...], axis=0), 0.0)
    assert np.allclose(np.std(X2[0, ...], axis=0), 0.0)


def test_evaluate_base_vectorize(ellipstell_state):
    X1_222 = ellipstell_state.evaluate_base("X1", [0.5, 0.6], [0, 0.1], [0, 0.1])
    X1_122 = ellipstell_state.evaluate_base("X1", [0.5], [0, 0.1], [0, 0.1])
    X1_212 = ellipstell_state.evaluate_base("X1", [0.5, 0.6], [0.1], [0, 0.1])
    X1_221 = ellipstell_state.evaluate_base("X1", [0.5, 0.6], [0, 0.1], [0.1])
    X1_112 = ellipstell_state.evaluate_base("X1", [0.5], [0.1], [0, 0.1])
    assert np.allclose(X1_222[0, :, :], X1_122[0, :, :])
    assert np.allclose(X1_222[:, 1, :], X1_212[:, 0, :])
    assert np.allclose(X1_222[:, :, 1], X1_221[:, :, 0])
    assert np.allclose(X1_222[0, 1, :], X1_112[0, 0, :])


def test_evaluate_base_bounds(ellipstell_state):
    rho = np.linspace(0, 1, 6)
    theta = np.linspace(0, 4 * np.pi, 8, endpoint=False)
    zeta = np.linspace(-2 * np.pi, 0, 10, endpoint=False)

    ellipstell_state.evaluate_base("X1", rho, theta, zeta)
    rho = np.linspace(-1, 1, 6)
    with pytest.raises(ValueError):
        ellipstell_state.evaluate_base("X1", rho, theta, zeta)
    rho = np.linspace(0, 2, 6)
    with pytest.raises(ValueError):
        ellipstell_state.evaluate_base("X1", rho, theta, zeta)


def test_evaluate_base_all(ellipstell_state):
    rho = np.linspace(0, 1, 6)
    theta = np.linspace(0, 2 * np.pi, 32, endpoint=False)
    zeta = np.linspace(0, 2 * np.pi, 10, endpoint=False)

    results = ellipstell_state.evaluate_base("all", rho, theta, zeta)
    assert len(results) == 8
    assert all(result.shape == (6, 32, 10) for result in results)
    with pytest.raises(ValueError):
        ellipstell_state.evaluate_base("all", rho, theta, zeta, selection="X1")


def test_evaluate_base_args(ellipstell_state):
    rho = np.linspace(0, 1, 6)
    theta = np.linspace(0, 2 * np.pi, 32, endpoint=False)
    zeta = np.linspace(0, 2 * np.pi, 10, endpoint=False)
    Z, T = np.meshgrid(zeta, theta)
    thetazeta = np.stack([T.flatten(), Z.flatten()], axis=0)

    # tensorproduct - positional
    X1p = ellipstell_state.evaluate_base("X1", rho, theta, zeta)
    # tensorproduct - keyword
    X1k = ellipstell_state.evaluate_base("X1", rho=rho, theta=theta, zeta=zeta)
    assert np.allclose(X1p, X1k)
    # list-tz - positional
    X1ptz = ellipstell_state.evaluate_base("X1", rho, thetazeta)
    assert np.allclose(X1p, X1ptz)
    # list-tz - keyword
    X1ktz = ellipstell_state.evaluate_base("X1", rho=rho, thetazeta=thetazeta)
    assert np.allclose(X1p, X1ktz)
    # wrong errors
    with pytest.raises(ValueError):
        ellipstell_state.evaluate_base("X1", rho, theta)
    with pytest.raises(ValueError):
        ellipstell_state.evaluate_base("X1", rho, zeta, zeta=zeta)
    with pytest.raises(ValueError):
        ellipstell_state.evaluate_base("X1", rho, thetazeta=thetazeta)
    with pytest.raises(ValueError):
        ellipstell_state.evaluate_base(
            "X1", rho=rho, theta=theta, zeta=zeta, thetazeta=thetazeta
        )
    with pytest.raises(ValueError):
        ellipstell_state.evaluate_base("X1", np.zeros((3, 32)))
    with pytest.raises(ValueError):
        ellipstell_state.evaluate_base("all", rho, theta, zeta)


def test_evaluate_base_xr(ellipstell_state):
    rho = np.linspace(0, 1, 6)
    theta = np.linspace(0, 2 * np.pi, 32, endpoint=False)
    zeta = np.linspace(0, 2 * np.pi, 10, endpoint=False)
    ds = xr.Dataset(coords={"rho": rho, "theta": theta, "zeta": zeta})

    dsall = ellipstell_state.evaluate_base("all", ds)
    assert all(
        [
            var in dsall
            for var in [
                "X1",
                "X2",
                "dX1_drho",
                "dX2_drho",
                "dX1_dtheta",
                "dX2_dtheta",
                "dX1_dzeta",
                "dX2_dzeta",
            ]
        ]
    )
    dsX1 = ellipstell_state.evaluate_base("X1", ds)
    assert "X1" in dsX1


def test_evaluate_hmap(ellipstell_state):
    rho = np.linspace(0, 1, 6)
    theta = np.linspace(0, 2 * np.pi, 32, endpoint=False)
    zeta = np.linspace(0, 2 * np.pi, 10, endpoint=False)
    R, T, Z = np.meshgrid(rho, theta, zeta, indexing="ij")
    # X1, X2, dX1_drho, dX2_drho, dX1_dtheta, dX2_dtheta, dX1_dzeta, dX2_dzeta
    inputs = ellipstell_state.evaluate_base("all", rho, theta, zeta)
    inputs = inputs[:2] + [T] + inputs[2:]

    outputs = ellipstell_state.evaluate_hmap(*inputs)
    assert len(outputs) == 4
    assert all(output.shape == (3, 6 * 32 * 10) for output in outputs)

def test_evaluate_profile(ellipstell_state):
    rho = np.linspace(0, 1, 6)
    ds = xr.Dataset(coords={"rho": rho})

    iota = ellipstell_state.evaluate_profile("iota", rho)
    assert iota.shape == rho.shape
    dp_drho = ellipstell_state.evaluate_profile("D_s p", rho)
    ds = ellipstell_state.evaluate_profile("Phi", ds)
    assert "Phi" in ds and ds.Phi.data.size == rho.size
