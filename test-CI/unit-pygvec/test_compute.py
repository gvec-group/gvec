"""
test the gvec.compute & gvec.quantities modules

In general this test suite assumes that test_state.py has passed.
"""

import pytest

try:
    import numpy as np
    import xarray as xr

    import gvec
    from gvec.state import State
    from gvec.comp import (
        Evaluations,
        compute,
        volume_integral,
    )
except ImportError:
    pytest.skip("Import Error", allow_module_level=True)

# === FIXTURES === #


@pytest.fixture()
def teststate(testfiles):
    with State(*testfiles) as state:
        yield state


@pytest.fixture()
def evals_r():
    rho = np.linspace(0, 1, 6)
    ds = Evaluations(rho=rho, theta=None, zeta=None)
    return ds


@pytest.fixture()
def evals_rtz():
    rho = np.linspace(0, 1, 6)
    theta = np.linspace(0, 2 * np.pi, 32, endpoint=False)
    zeta = np.linspace(0, 2 * np.pi, 10, endpoint=False)
    ds = Evaluations(rho=rho, theta=theta, zeta=zeta)
    return ds


@pytest.fixture()
def evals_rtz_int(teststate):
    ds = Evaluations(
        state=teststate,
        rho="int",
        theta="int",
        zeta="int",
    )
    return ds


@pytest.fixture()
def evals_r_int_tz(teststate):
    ds = Evaluations(
        state=teststate,
        rho="int",
        theta=100,
        zeta=100,
    )
    return ds


@pytest.fixture()
def evals_r_tz_list(teststate):
    rho = np.linspace(0, 1, 6)
    theta = np.linspace(0, 2 * np.pi, 32, endpoint=False)
    zeta = np.linspace(0, 2 * np.pi, 10, endpoint=False)
    R, T, Z = np.meshgrid(rho, theta, zeta, indexing="ij")
    dims = ("rad", "pol", "tor")
    ds = xr.Dataset(
        coords=dict(rho=("rad", rho), theta=(dims, T), zeta=(dims, Z)),
    )
    return ds


@pytest.fixture(
    params=["rtz", "r-tz-list"],
)
def evals_rtz_and_list(request, evals_rtz, evals_r_tz_list):
    if request.param == "rtz":
        return evals_rtz
    elif request.param == "r-tz-list":
        return evals_r_tz_list


# === TESTS === #


def test_evaluations_init(teststate):
    ds = Evaluations(rho=[0.5, 0.6], theta=None, zeta=None)
    assert np.allclose(ds.rho, [0.5, 0.6])
    assert {"rho"} == set(ds.coords)

    with pytest.raises(ValueError):
        ds = Evaluations(rho="int")
    with pytest.raises(ValueError):
        ds = Evaluations(theta="int")
    with pytest.raises(ValueError):
        ds = Evaluations(zeta="int")

    ds = Evaluations(rho=np.array([0.5, 0.6]), theta="int", zeta="int", state=teststate)
    assert np.allclose(ds.rho, [0.5, 0.6])
    assert {"rho", "theta", "zeta", "pol_weight", "tor_weight"} == set(ds.coords)
    assert ds.rho.attrs["integration_points"] is False

    ds = Evaluations(rho=[0.5], theta="int", zeta="int", state=teststate)
    assert np.allclose(ds.rho, [0.5])
    assert {"rho", "theta", "zeta", "pol_weight", "tor_weight"} == set(ds.coords)
    assert ds.rho.attrs["integration_points"] is False

    ds = Evaluations(rho="int", theta=None, zeta=None, state=teststate)
    assert {"rho", "rad_weight"} == set(ds.coords)
    assert ds.rho.attrs["integration_points"] is True

    ds = Evaluations("int", "int", "int", state=teststate)
    for c in ["rho", "theta", "zeta"]:
        assert c in ds.coords
        assert ds[c].attrs["integration_points"] is True


def test_compute_base(teststate, evals_rtz_and_list):
    ds = evals_rtz_and_list

    compute(ds, "X1", state=teststate)
    assert "X1" in ds
    assert ds.X1.shape == (6, 32, 10)
    assert "dX1_dr" in ds
    assert not np.any(np.isnan(ds.X1))

    compute(ds, "X2", "LA", state=teststate)
    assert "X1" in ds and "X2" in ds and "LA" in ds
    assert "dX2_dt" in ds
    assert "dLA_dz" in ds


def test_compute_hmap(teststate, evals_rtz):
    ds = evals_rtz
    compute(ds, "pos", state=teststate)

    assert ds.pos.shape == (3, 6, 32, 10)
    assert "xyz" in ds.dims
    assert "xyz" in ds.coords
    assert "e_X1" in ds
    assert "e_X2" in ds
    assert "e_zeta3" in ds
    assert not np.any(np.isnan(ds.pos))

    compute(ds, "e_rho", state=teststate)
    assert ds.e_rho.shape == (3, 6, 32, 10)


def test_compute_metric(teststate, evals_rtz):
    ds = evals_rtz

    compute(ds, "g_tt", state=teststate)
    for ij in ("rr", "rt", "rz", "tt", "tz", "zz"):
        key = f"g_{ij}"
        assert key in ds
        assert set(ds[key].coords) == {"rho", "theta", "zeta"}
        for k in "rtz":
            assert f"dg_{ij}_d{k}" in ds

    compute(ds, "e_rho", "e_theta", "e_zeta", state=teststate)
    idxs = {"r": "rho", "t": "theta", "z": "zeta"}
    for ij in "rr rt rz tt tz zz".split():
        key = f"g_{ij}"
        assert np.allclose(
            ds[f"g_{ij}"],
            xr.dot(ds[f"e_{idxs[ij[0]]}"], ds[f"e_{idxs[ij[1]]}"], dim="xyz"),
        )


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
def test_compute_profile(teststate, evals_r, quantity):
    ds = evals_r
    compute(ds, quantity, state=teststate)
    assert quantity in ds
    assert set(ds[quantity].coords) == {"rho"}
    assert ds[quantity].data.size == ds.rho.size
    assert not np.any(np.isnan(ds[quantity]))


@pytest.mark.parametrize(
    "quantity,on_axis",
    [
        ("Jac", True),
        ("B", False),
        ("J", False),
        ("mod_B", False),
        ("mod_J", False),
    ],
)
def test_compute_derived(teststate, evals_rtz, quantity, on_axis):
    ds = evals_rtz
    compute(ds, quantity, state=teststate)
    assert quantity in ds
    if on_axis:
        assert not np.any(np.isnan(ds[quantity]))
    else:
        assert np.all(np.isnan(ds[quantity].isel(rad=0)))
        assert not np.any(np.isnan(ds[quantity].isel(rad=slice(1, None))))


def test_compute_basis(teststate, evals_rtz):
    ds = evals_rtz
    for coord in ["rho", "theta", "zeta"]:
        compute(ds, f"e_{coord}", f"grad_{coord}", state=teststate)
    ds = ds.isel(rad=slice(1, None))
    for coord in ["rho", "theta", "zeta"]:
        assert np.allclose(
            xr.dot(ds[f"e_{coord}"], ds[f"grad_{coord}"], dim="xyz"), 1.0
        )
        for coord2 in ["rho", "theta", "zeta"]:
            if coord2 == coord:
                continue
            compute(ds, f"e_{coord2}", f"grad_{coord2}", state=teststate)
            assert np.allclose(
                xr.dot(ds[f"e_{coord}"], ds[f"grad_{coord2}"], dim="xyz"), 0.0
            )
            assert np.allclose(
                xr.dot(ds[f"grad_{coord}"], ds[f"e_{coord2}"], dim="xyz"), 0.0
            )


def test_volume_integral(teststate, evals_rtz_int, evals_rtz):
    ds = evals_rtz_int
    compute(ds, "Jac", state=teststate)
    volume_integral(ds.Jac)

    with pytest.raises(ValueError):
        compute(evals_rtz, "Jac", state=teststate)
        volume_integral(evals_rtz.Jac)


@pytest.mark.parametrize(
    "quantity",
    [
        "V",
        "dV_dPhi_n",
        "dV_dPhi_n2",
        "minor_radius",
        "major_radius",
        "iota_mean",
        "iota_tor",
        "I_tor",
        "I_pol",
        "B_theta_avg",
    ],
)
def test_integral_quantities(teststate, evals_rtz_int, quantity):
    ds = evals_rtz_int
    compute(ds, quantity, state=teststate)
    # --- check that metadata is preserved --- #
    assert ds.rho.attrs["integration_points"] is True


@pytest.mark.xfail()
@pytest.mark.parametrize(
    "quantity",
    [
        "V",
        "dV_dPhi_n",
        "dV_dPhi_n2",
        "minor_radius",
        "major_radius",
        "iota_mean",
        "iota_tor",
        "I_tor",
        "I_pol",
        "B_theta_avg",
    ],
)
def test_integral_quantities_aux(teststate, evals_r_int_tz, quantity):
    # automatic change to integration points if necessary
    ds = evals_r_int_tz
    compute(ds, quantity, state=teststate)
    assert ds.rho.attrs["integration_points"] is True
    assert ds.theta.attrs["integration_points"] is False
    assert ds.zeta.attrs["integration_points"] is False
    assert "theta_weight" not in ds and "zeta_weight" not in ds
    assert "rho_weights" in ds


def test_table_of_quantities():
    s = gvec.comp.table_of_quantities()
    assert len(s.split("\n")) == 3 + len(
        gvec.comp.QUANTITIES
    )  # 2 lines header, 1 trailing \n
