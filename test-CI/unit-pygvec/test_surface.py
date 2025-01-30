import pytest

try:
    import numpy as np
    from numpy.random import random
    import xarray as xr

    from gvec import State, EvaluationsBoozer, compute
    from gvec import surface
except ImportError:
    pass  # tests will be skipped via the `check_import` fixture


# === FIXTURES === #


@pytest.fixture()
def state(testfiles):
    with State(*testfiles) as state:
        yield state


@pytest.fixture()
def ev(state):
    rho = [0.1, 0.5, 0.9]
    ds = EvaluationsBoozer(rho=rho, n_theta=20, n_zeta=50, state=state)
    compute(ds, "pos", "N_FP", state=state)
    return ds


@pytest.fixture()
def surfs(ev):
    return surface.init_surface(ev.pos, ev.N_FP)


# === TESTS === #


def test_init_surface_single(ev):
    evs = ev.isel(rad=0)
    surf = surface.init_surface(evs.pos, evs.N_FP)
    assert isinstance(surf, xr.Dataset)
    assert {"pol", "tor"} <= set(surf.dims)
    assert {"xhat", "yhat", "zhat", "dxhat_dt", "dzhat_dtz"} < set(surf.data_vars)


def test_init_surface_multiple(ev):
    surf = surface.init_surface(ev.pos, ev.N_FP)
    assert {"rad", "pol", "tor"} <= set(surf.dims)
    assert {"xhat", "yhat", "zhat", "dxhat_dt", "dzhat_dtz"} < set(surf.data_vars)


@pytest.mark.parametrize(
    "Q",
    [
        "pos",
        "e_theta_B",
        "e_zeta_B",
        "k_tt_B",
        "k_tz_B",
        "k_zz_B",
        "g_tt_B",
        "g_tz_B",
        "g_zz_B",
        "normal",
        "II_tt_B",
        "II_tz_B",
        "II_zz_B",
    ],
)
def test_compute(surfs, Q):
    surface.compute(surfs, Q)
    assert Q in surfs.data_vars
    assert np.isnan(surfs[Q].data).sum() == 0
