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


def test_evaluate(ellipstell_state):
    s = np.linspace(0.1, 0.9, 6)
    theta = np.linspace(0, 2 * np.pi, 8, endpoint=False)
    zeta = np.linspace(0, 2 * np.pi, 10, endpoint=False)

    result = ellipstell_state.evaluate_base(s, theta, zeta, "X1")
    assert result.shape == (6, 8, 10)
    print(result)
    
