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


@pytest.fixture(scope="session")
def testcase():
    return "ellipstell_lowres"


@pytest.fixture(scope="session")
def state(testfiles):
    return State(*testfiles)


@pytest.fixture(params=["eval", "fft"])
def ift(request):
    return request.param


@pytest.fixture(scope="session")
def ev_boozer(state):
    rho = [0.1, 0.5, 0.9]
    ds = state.evaluate_sfl("pos", "N_FP", sfl="boozer", rho=rho, theta=21, zeta=31, MNfactor=2)
    return ds


@pytest.fixture()
def ev(ev_boozer):
    # only for testing! MNfactor should be ~4 or above for most applications
    return ev_boozer.copy(deep=True)


@pytest.fixture(scope="session")
def ev_boozer_highres(state):
    rho = [0.1, 0.5, 0.9]
    ds = state.evaluate_sfl(
        "pos", "N_FP", sfl="boozer", rho=rho, theta=101, zeta=101, MNfactor=2
    )
    return ds


@pytest.fixture()
def ev_highres(ev_boozer_highres):
    # only for testing! MNfactor should be ~4 or above for most applications
    return ev_boozer_highres.copy(deep=True)


@pytest.fixture(params=[0, 1, -1])
def winding(request, state):
    return 1 + request.param * state.nfp


@pytest.fixture()
def surfs(ev, ift, winding):
    return surface.init_surface(ev.pos, ev.N_FP, ift=ift, winding=winding)


@pytest.fixture()
def surfs_highres(ev_highres, ift, winding):
    return surface.init_surface(ev_highres.pos, ev_highres.N_FP, ift=ift, winding=winding)


# === TESTS === #


def test_init_surface_single(ev, ift):
    evs = ev.isel(rad=0)
    surf = surface.init_surface(evs.pos, evs.N_FP, ift=ift)
    assert isinstance(surf, xr.Dataset)
    assert {"pol", "tor"} <= set(surf.dims)
    assert {"xhat", "yhat", "zhat", "dxhat_dt", "dzhat_dtz"} < set(surf.data_vars)


def test_init_surface_multiple(ev, ift):
    surf = surface.init_surface(ev.pos, ev.N_FP, ift=ift)
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


@pytest.mark.parametrize("ift", ["fft"])
@pytest.mark.parametrize(
    "Q",
    [
        "pos",
        "e_theta_B",
        "e_zeta_B",
        "g_tt_B",
        "g_tz_B",
        "g_zz_B",
        "normal",
        "k_tt_B",
        "k_tz_B",
        "k_zz_B",
        "II_tt_B",
        "II_tz_B",
        "II_zz_B",
    ],
)
def test_compare(state, ev_highres, surfs_highres, Q, ift):
    """Test that the computed surface quantities match the expected values."""
    surface.compute(surfs_highres, Q)
    compute(ev_highres, Q, state=state)

    surfQ, evQ = xr.broadcast(surfs_highres[Q], ev_highres[Q])
    np.testing.assert_allclose(surfQ.data, evQ.data, rtol=1e-2, atol=1e-3)
