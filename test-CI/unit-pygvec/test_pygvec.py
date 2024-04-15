import pytest
import numpy as np
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
    X1 = ellipstell_state.evaluate_base(rho, theta, zeta, "X1")
    X2 = ellipstell_state.evaluate_base(rho, theta, zeta, "X2")
    LA = ellipstell_state.evaluate_base(rho, theta, zeta, "LA")
    D_rtzX2 = ellipstell_state.evaluate_base(rho, theta, zeta, "D_rtz X2")

    assert X1.shape == (6, 32, 10)
    assert not np.any(np.isnan(X1))
    assert not np.any(np.isnan(X2))
    assert not np.any(np.isnan(LA))
    # magnetic axis collapses to a line
    assert np.allclose(np.std(X1[0, ...], axis=0), 0.0)
    assert np.allclose(np.std(X2[0, ...], axis=0), 0.0)


def test_evaluate_base_vectorize(ellipstell_state):
    X1_222 = ellipstell_state.evaluate_base([.5, .6], [0, .1], [0, .1], "X1")
    X1_122 = ellipstell_state.evaluate_base([.5    ], [0, .1], [0, .1], "X1")
    X1_212 = ellipstell_state.evaluate_base([.5, .6], [   .1], [0, .1], "X1")
    X1_221 = ellipstell_state.evaluate_base([.5, .6], [0, .1], [   .1], "X1")
    X1_112 = ellipstell_state.evaluate_base([.5    ], [   .1], [0, .1], "X1")
    assert np.allclose(X1_222[0,:,:], X1_122[0,:,:])
    assert np.allclose(X1_222[:,1,:], X1_212[:,0,:])
    assert np.allclose(X1_222[:,:,1], X1_221[:,:,0])
    assert np.allclose(X1_222[0,1,:], X1_112[0,0,:])


def test_evaluate_base_bounds(ellipstell_state):
    s = np.linspace(0, 1, 6)
    theta = np.linspace(0, 4 * np.pi, 8, endpoint=False)
    zeta = np.linspace(-2 * np.pi, 0, 10, endpoint=False)

    X1 = ellipstell_state.evaluate_base(s, theta, zeta, "X1")

    s = np.linspace(-1, 1, 6)
    with pytest.raises(ValueError):
        ellipstell_state.evaluate_base(s, theta, zeta, "X1")
    
    s = np.linspace(0, 2, 6)
    with pytest.raises(ValueError):
        ellipstell_state.evaluate_base(s, theta, zeta, "X1")
