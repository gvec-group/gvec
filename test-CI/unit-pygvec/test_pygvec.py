import pytest
import importlib
import numpy as np
import os
from pathlib import Path

import helpers

# === Fixtures === #


@pytest.fixture()
def pygvec():
    try:
        pygvec = importlib.reload(pygvec)
    except NameError:
        import pygvec
    return pygvec

@pytest.fixture()
def ellipstell_tmpdir(tmpdir):
    """prepare the ellipstell parameters"""
    refdir = Path("test-CI/run/example/ellipstell_lowres").absolute()
    paramfile = "parameter.ini"
    statefile = "ELLIPSTELL_LOWRES_State_0000_00000000.dat"
    with helpers.chdir(tmpdir):
        os.symlink(refdir / paramfile, paramfile)
        os.symlink(refdir / statefile, statefile)
        yield

@pytest.fixture()
def ellipstell_params(pygvec, ellipstell_tmpdir):
    paramfile = "parameter.ini"
    # pygvec.main.ReadParameterFile(paramfile)
    yield

# === Tests === #


def test_version(pygvec):
    import pkg_resources

    assert isinstance(pygvec.__version__, str)
    assert pygvec.__version__ == pygvec.version
    assert pygvec.__version_tuple__ >= (0, 2, 1)
    assert pkg_resources.get_distribution("pygvec").version == pygvec.__version__
