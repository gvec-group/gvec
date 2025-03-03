import pytest
from pathlib import Path

try:
    import numpy as np

    import xarray as xr
    import gvec
    from gvec.state import State
except ImportError:
    pytest.skip("Import Error", allow_module_level=True)


tol = 1e-12


# Utility functions for transforming polynomials into B-Splines.
# ============================================================== #
def binc(n, k):
    """Binomial coefficient implementation to avoid using scipy.

    Args:
        n (int): top value
        k (int): bottom value

    Returns:
        float: binomial coefficient (n k)
    """
    if k > k:
        AssertionError("error in binomial coefficient (n k) calculation: k>n")
    if abs(k) <= tol:
        coef = 1.0
    elif k > n / 2:
        coef = binc(n, n - k)
    else:
        coef = n * binc(n - 1, k - 1) / k
    return coef


def dual_fatcor(r: int, d: int, j: int, t):
    """Prefactor for the transformation of polynomial coefficients to B-spline coefficients.
    Only works on the interval s=rho^2 in [0,1].

    Args:
        r (int): polynomial degree.
        d (int): max polynomial degree.
        j (int): index of the B-spline.
        t (np.ArrayLike): B-spline knots.

    Returns:
        float: prefactor from the dual-polynomial.
    """
    n = int(np.sum(t[j : j + d + 1]))
    if r > n:
        return 0
    c_jr = binc(n, r) / binc(d, r)
    return c_jr


def poly2bspl_coeff(c, j, t):
    """

    Args:
        c (int): polynomial coefficients
        j (int): index of the B-spline.
        t (int): B-spline knots.

    Returns:
        float: B-spline coefficient.
    """
    d = len(c) - 1
    c_spl = 0
    for r in range(d + 1):
        c_jr = dual_fatcor(r, d, j, t)
        c_spl += c[r] * c_jr
    return c_spl


# === Fixtures === #


@pytest.fixture(scope="module")
def c_poly():
    """Reference polynomial to test against.

    Returns:
        list: polynomial coefficients
    """
    return [1, 2, 3, -1, -2, -3]


@pytest.fixture(scope="module", params=["poly", "bspl"])
def testfile_aux(request, testcaserundir, c_poly):
    """prepare the testcase parameters"""
    params_gvec = {"sign_iota": 1, "pres_scale": 1}
    deg = len(c_poly) - 1
    if request.param == "poly":
        paramfile = "parameter_poly.ini"
        params_gvec["pres_coefs"] = gvec.util.np2gvec(c_poly)
        params_gvec["pres_type"] = "polynomial"
        params_gvec["iota_coefs"] = gvec.util.np2gvec(c_poly)
        params_gvec["iota_type"] = "polynomial"
    elif request.param == "bspl":
        paramfile = "parameter_bspl.ini"
        c_bspl = np.zeros(deg + 1)
        knots = np.concatenate([np.zeros(deg + 1), np.ones(deg + 1)])
        for j in range(deg + 1):
            c_bspl[j] = poly2bspl_coeff(c_poly, j, knots)
        params_gvec = gvec.util.bspl2gvec(
            "pres", knots=knots, coefs=c_bspl, params=params_gvec
        )
        params_gvec = gvec.util.bspl2gvec(
            "iota", knots=knots, coefs=c_bspl, params=params_gvec
        )
    gvec.util.adapt_parameter_file(
        testcaserundir / "parameter.ini", testcaserundir / paramfile, **params_gvec
    )
    yield (testcaserundir / paramfile)


@pytest.fixture()
def teststate(testfile_aux):
    with State(testfile_aux) as state:
        yield state


@pytest.fixture(scope="module")
def vmecfiles():
    testcase = "w7x_from_vmec_initLA_F"
    path_example = Path(__file__).parent / "../examples/" / testcase
    return (path_example / "parameter.ini", path_example / "wout_d23p4_tm.nc")


# === Tests === #
@pytest.mark.parametrize("type", ["p", "iota"])
def test_eval_profile(teststate, type, c_poly):
    rho = np.linspace(0, 1, 100)
    ref_poly = np.polynomial.Polynomial(c_poly)
    ref_poly_vals = ref_poly(rho**2)
    gvec_profile = teststate.evaluate_profile(type, rho=rho)
    np.testing.assert_allclose(gvec_profile, ref_poly_vals)


@pytest.mark.parametrize("type", ["p", "iota"])
def test_eval_profile_prime(teststate, type, c_poly):
    rho = np.linspace(0, 1, 100)
    ref_poly = np.polynomial.Polynomial(c_poly)
    ref_poly_vals = ref_poly.deriv(m=1)(rho**2) * 2 * rho
    gvec_profile = teststate.evaluate_profile(type, rho=rho, deriv=1)
    np.testing.assert_allclose(gvec_profile, ref_poly_vals)


@pytest.mark.parametrize("type", ["p", "iota"])
def test_eval_profile_n_deriv(teststate, type, c_poly):
    rho = np.linspace(0, 1, 100)
    ref_poly = np.polynomial.Polynomial(c_poly)

    # derivatives of s=rho**2 with respect to rho
    ds1 = 2 * rho
    ds2 = 2

    # derivatives of the test-polynomial with respect to s=rho**2
    drefs1 = ref_poly.deriv(m=1)(rho**2)
    drefs2 = ref_poly.deriv(m=2)(rho**2)
    drefs3 = ref_poly.deriv(m=3)(rho**2)
    drefs4 = ref_poly.deriv(m=4)(rho**2)

    # derivatives of the test-polynomial with respect to rho
    dref1 = drefs1 * ds1
    dref2 = ds1**2 * drefs2 + drefs1 * ds2
    dref3 = 3 * ds1 * drefs2 * ds2 + ds1**3 * drefs3
    dref4 = ds1**4 * drefs4 + 6 * ds2 * ds1**2 * drefs3 + 3 * ds2**2 * drefs2

    for i, dref in enumerate([dref1, dref2, dref3, dref4]):
        gvec_profile_di = teststate.evaluate_profile(type, rho=rho, deriv=i + 1)
        np.testing.assert_allclose(dref, gvec_profile_di)


def test_eval_profile_iota_vs_phi_and_chi(teststate, c_poly):
    """
    since chi and phi are computed using an anti-derivative from dPhi/ds and dchi/ds,
    this is a test for the anti-derivative.
    """
    rho = np.linspace(1.0e-16, 1.0, 131)
    ref_poly = np.polynomial.Polynomial(c_poly)
    ref_poly_vals = ref_poly(rho**2)
    eval_iota = teststate.evaluate_profile("iota", rho=rho)
    np.testing.assert_allclose(eval_iota, ref_poly_vals, atol=tol)
    eval_phiprime = teststate.evaluate_profile("Phi", rho=rho, deriv=1)
    eval_chiprime = teststate.evaluate_profile("chi", rho=rho, deriv=1)
    np.testing.assert_allclose(eval_chiprime / eval_phiprime, ref_poly_vals, atol=tol)


def test_eval_rho2_profile_iota_vs_phi_and_chi(teststate, c_poly):
    """
    since chi and phi are computed using an anti-derivative from dPhi/ds and dchi/ds,
    this is a test for the anti-derivative.
    """
    rho2 = np.linspace(1.0e-16, 1.0, 129)
    ref_poly = np.polynomial.Polynomial(c_poly)
    ref_poly_vals = ref_poly(rho2)
    eval_iota = teststate.evaluate_rho2_profile("iota", rho2=rho2)
    np.testing.assert_allclose(eval_iota, ref_poly_vals, atol=tol)
    eval_dphids = teststate.evaluate_rho2_profile("Phi", rho2=rho2, deriv=1)
    eval_dchids = teststate.evaluate_rho2_profile("chi", rho2=rho2, deriv=1)
    np.testing.assert_allclose(eval_dchids / eval_dphids, ref_poly_vals, atol=tol)


@pytest.mark.parametrize("BC_type_axis", ["not_a_knot", "1st_deriv", "2nd_deriv"])
@pytest.mark.parametrize("BC_type_edge", ["not_a_knot", "1st_deriv", "2nd_deriv"])
@pytest.mark.parametrize("n_points", [4, 5, 6, 11])
def test_interpolation(testcaserundir, c_poly, BC_type_axis, BC_type_edge, n_points):
    pres_scale = 1500
    params_gvec = {"sign_iota": 1, "pres_scale": pres_scale}
    cubic_poly_c = c_poly[:4]
    paramfile = "parameter_interpolation.ini"

    cubic_poly = np.polynomial.Polynomial(cubic_poly_c)
    rho2_vals = np.linspace(0, 1, n_points)
    P_vals = cubic_poly(rho2_vals)
    iota_vals = cubic_poly(rho2_vals)

    params_gvec["pres_type"] = "interpolation"
    params_gvec["pres_rho2"] = gvec.util.np2gvec(rho2_vals)
    params_gvec["pres_vals"] = gvec.util.np2gvec(P_vals)
    params_gvec["pres_BC_type_axis"] = BC_type_axis
    params_gvec["pres_BC_type_edge"] = BC_type_edge

    params_gvec["iota_type"] = "interpolation"
    params_gvec["iota_rho2"] = gvec.util.np2gvec(rho2_vals)
    params_gvec["iota_vals"] = gvec.util.np2gvec(iota_vals)
    params_gvec["iota_BC_type_axis"] = BC_type_axis
    params_gvec["iota_BC_type_edge"] = BC_type_edge

    rho2 = np.linspace(0, 1, 200)
    P = pres_scale * cubic_poly(rho2)
    iota = cubic_poly(rho2)

    gvec.util.adapt_parameter_file(
        testcaserundir / "parameter.ini", testcaserundir / paramfile, **params_gvec
    )
    with State(testcaserundir / paramfile) as state:
        P_p_interpol = state.evaluate_rho2_profile("p", rho2_vals, deriv=0)
        iota_p_interpol = state.evaluate_rho2_profile("iota", rho2_vals, deriv=0)

        P_interpol = state.evaluate_rho2_profile("p", rho2, deriv=0)
        iota_interpol = state.evaluate_rho2_profile("iota", rho2, deriv=0)

        dP_interpol_ds = state.evaluate_rho2_profile("p", rho2, deriv=1)
        diota_interpol_ds = state.evaluate_rho2_profile("iota", rho2, deriv=1)

        dP_interpol_dss = state.evaluate_rho2_profile("p", rho2, deriv=2)
        diota_interpol_dss = state.evaluate_rho2_profile("iota", rho2, deriv=2)

    # check that the profiles are the same at all provided interpolation points
    np.testing.assert_allclose(P_p_interpol / pres_scale, P_vals)
    np.testing.assert_allclose(iota_p_interpol, iota_vals)

    # check that the profiles are the same if not-a-knot BC are used.
    if BC_type_axis == "not_a_knot" and BC_type_edge == "not_a_knot":
        np.testing.assert_allclose(P_interpol, P)
        np.testing.assert_allclose(iota_interpol, iota)

    # check that the BC-are enforced.
    if BC_type_axis == "1st_deriv":
        np.testing.assert_allclose(dP_interpol_ds[0] / pres_scale, 0, atol=tol)
        np.testing.assert_allclose(diota_interpol_ds[0], 0, atol=tol)

    if BC_type_axis == "2nd_deriv":
        np.testing.assert_allclose(dP_interpol_dss[0] / pres_scale, 0, atol=tol)
        np.testing.assert_allclose(diota_interpol_dss[0], 0, atol=tol)

    if BC_type_edge == "1st_deriv":
        np.testing.assert_allclose(dP_interpol_ds[-1] / pres_scale, 0, atol=tol)
        np.testing.assert_allclose(diota_interpol_ds[-1], 0, atol=tol)

    if BC_type_edge == "2nd_deriv":
        np.testing.assert_allclose(dP_interpol_dss[-1] / pres_scale, 0, atol=tol)
        np.testing.assert_allclose(diota_interpol_dss[-1], 0, atol=tol)


def test_vmec_profile_init(vmecfiles):
    wout = xr.open_dataset(vmecfiles[1])
    phi = wout.phi
    P_vmec = wout.presf
    iota_vmec = -1 * wout.iotaf
    rho_wout = np.sqrt((phi - phi[0]) / (phi[-1] - phi[0]))  # Normalized toroidal flux

    with State(vmecfiles[0]) as state:
        P_gvec = state.evaluate_profile("p", rho_wout, deriv=0)
        dP_dr_gvec = state.evaluate_profile("p", 0, deriv=1)
        iota_gvec = state.evaluate_profile("iota", rho_wout, deriv=0)
        diota_dr_gvec = state.evaluate_profile("iota", 0, deriv=1)

    np.testing.assert_allclose(P_gvec, P_vmec)
    np.testing.assert_allclose(iota_gvec, iota_vmec)
    np.testing.assert_allclose(diota_dr_gvec, 0, atol=tol)
    np.testing.assert_allclose(dP_dr_gvec, 0, atol=tol)


@pytest.mark.parametrize("profile_type", ["pres", "iota"])
def test_vmec_with_custom_profile(testcaserundir, vmecfiles, c_poly, profile_type):
    wout = xr.open_dataset(vmecfiles[1])
    phi = wout.phi
    rho_wout = np.sqrt((phi - phi[0]) / (phi[-1] - phi[0]))  # Normalized toroidal flux

    params = {"vmecwoutFile": vmecfiles[1]}
    if profile_type == "pres":
        var = "p"
        scale = 1
        params["init_with_profile_pressure"] = "T"
        params["pres_scale"] = scale
    else:
        var = "iota"
        scale = -1
        params["init_with_profile_iota"] = "T"
        params["iota_sign"] = scale

    params[f"{profile_type}_coefs"] = gvec.util.np2gvec(c_poly)
    params[f"{profile_type}_type"] = "polynomial"
    param_file = testcaserundir / f"parameter_vmec_custom_{profile_type}.ini"
    gvec.util.adapt_parameter_file(vmecfiles[0], param_file, **params)
    rho = np.linspace(0, 1, 100)
    ref_poly = np.polynomial.Polynomial(c_poly)
    profile_true = scale * ref_poly(rho**2)
    with State(param_file) as state:
        profile_gvec = state.evaluate_profile(var, rho, deriv=0)
        profile_gvec_at_vmec = state.evaluate_profile(var, rho_wout, deriv=0)

    with np.testing.assert_raises(AssertionError):
        np.testing.assert_allclose(profile_gvec_at_vmec, wout[f"{profile_type}f"])
    np.testing.assert_allclose(profile_gvec, profile_true)
