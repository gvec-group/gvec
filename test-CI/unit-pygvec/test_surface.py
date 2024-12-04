import pytest

try:
    import numpy as np
    from numpy.random import random
    import xarray as xr

    from gvec import State, Evaluations, compute
    from gvec import surface
except ImportError:
    pytest.skip("Import Error", allow_module_level=True)


# === FIXTURES === #


@pytest.fixture()
def state(testfiles):
    with State(*testfiles) as state:
        yield state


@pytest.fixture()
def ev(state):
    rho = [0.1, 0.5, 0.9]
    ds = Evaluations(rho=rho, theta=20, zeta=20, state=state)
    compute(ds, "pos", "N_FP", state=state)
    return ds


@pytest.fixture()
def surfs(ev):
    x = ev.pos.loc["x"].data
    y = ev.pos.loc["y"].data
    z = ev.pos.loc["z"].data
    nfp = ev.N_FP.item()
    return surface.init_surface(x, y, z, nfp)


# === TESTS === #


def test_init_surface_single(ev):
    evs = ev.isel(rad=0)
    x = evs.pos.loc["x"].data
    y = evs.pos.loc["y"].data
    z = evs.pos.loc["z"].data
    nfp = evs.N_FP.item()

    surf = surface.init_surface(x, y, z, nfp)
    assert isinstance(surf, xr.Dataset)
    assert {"pol", "tor"} == set(surf.dims)
    assert {"xhat", "yhat", "zhat", "dxhat_dt", "dzhat_dtz"} < set(surf.data_vars)


def test_init_surface_multiple(ev):
    x = ev.pos.loc["x"].data
    y = ev.pos.loc["y"].data
    z = ev.pos.loc["z"].data
    nfp = ev.N_FP.item()

    surf = surface.init_surface(x, y, z, nfp)
    assert {"rad", "pol", "tor"} == set(surf.dims)
    assert {"xhat", "yhat", "zhat", "dxhat_dt", "dzhat_dtz"} < set(surf.data_vars)


@pytest.mark.parametrize(
    "Q",
    [
        "pos",
        "e_theta",
        "e_zeta",
        "k_tt",
        "k_tz",
        "k_zz",
        "g_tt",
        "g_tz",
        "g_zz",
        "normal",
        "II_tt",
        "II_tz",
        "II_zz",
    ],
)
def test_compute(surfs, Q):
    surface.compute(surfs, Q)
    assert Q in surfs.data_vars
    assert np.isnan(surfs[Q].data).sum() == 0


@pytest.mark.parametrize("Q", ["pos", "e_theta", "e_zeta", "g_tt", "g_tz", "g_zz"])
def test_compare(ev, state, surfs, Q):
    surface.compute(surfs, Q)
    compute(ev, Q, state=state)
    np.testing.assert_allclose(ev[Q].data, surfs[Q].data, atol=1e-12)
