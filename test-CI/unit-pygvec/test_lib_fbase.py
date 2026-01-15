import pytest

try:
    import numpy as np

    import gvec.lib
    from gvec.lib import modgvec_fbase as mod
except ImportError:
    pytest.skip("Import Error", allow_module_level=True)


# === Fixtures === #


@pytest.fixture(scope="module", autouse=True)
def redirect_abort():
    """Redirect abort to raise an exception instead."""
    gvec.lib.modgvec_py_binding.redirect_abort()


@pytest.fixture()
def fbase():
    """Fixture to create a fBase object for testing."""
    return mod.t_fBase((3, 4), (7, 9), 2, "_sin_", True)


# === Tests === #


def test_lib_fbase_init():
    fbase = mod.t_fBase((3, 4), (7, 9), 2, "_sin_", True)
    assert fbase._alloc
    del fbase


def test_lib_fbase_properties():
    fbase = mod.t_fBase((3, 4), (7, 9), 2, "_sin_", True)

    assert all(fbase.mn_max == (3, 4))
    assert all(fbase.mn_nyq == (7, 9))
    assert fbase.nfp == 2
    assert fbase.sin_cos == 1
    assert fbase.exclude_mn_zero

    assert fbase.modes == 31
    assert all(fbase.cos_range == (0, 0))
    assert all(fbase.sin_range == (0, 31))

    fbase.d_thet
    fbase.d_zeta
    fbase.thet_ip
    fbase.zeta_ip
    fbase.xmn


def test_lib_fbase_eval(fbase):
    assert fbase.eval(0, np.array([0.0, 1.0])).shape == (fbase.modes,)


def test_lib_fbase_eval_xn(fbase):
    t = np.linspace(0, 2 * np.pi, 3)
    z = np.linspace(0, 2 * np.pi, 4)
    tgrid, zgrid = np.meshgrid(t, z, indexing="ij")
    tz = np.asfortranarray(np.vstack([tgrid.ravel(), zgrid.ravel()]))
    assert fbase.eval_xn(0, tz.shape[1], tz).shape == (12, fbase.modes)
