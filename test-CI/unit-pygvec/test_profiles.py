import pytest
from pathlib import Path

try:
    import numpy as np
    import xarray as xr
    from numpy.polynomial import Polynomial
    from scipy.interpolate import BSpline

    import gvec
    from gvec import State
    from gvec.core.minimize import MinimizationState
except ImportError as e:
    pytest.skip(f"Import Error: {e}", allow_module_level=True)

EXAMPLES = Path(__file__).parent / "../examples"
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
    if k > n:
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


@pytest.fixture(scope="session")
def testcase():
    return "ellipstell_lowres"


@pytest.fixture()
def woutfile(testcase):
    return EXAMPLES / testcase / "wout_d23p4_tm.nc"


@pytest.fixture()
def parameters(testcase, woutfile):
    params = gvec.util.read_parameters(EXAMPLES / testcase / "parameter.ini")
    if woutfile.exists():
        params["vmecwoutFile"] = woutfile.absolute()
    return params


@pytest.fixture(params=["poly", "bspl"])
def teststate(request, c_poly, testcase, tmp_path):
    """prepare the testcase parameters"""
    parameters = gvec.util.read_parameters(EXAMPLES / testcase / "parameter.ini")
    deg = len(c_poly) - 1
    if request.param == "poly":
        parameters["pres"] = {"type": "polynomial", "coefs": c_poly}
        parameters["iota"] = {"type": "polynomial", "coefs": c_poly}
    elif request.param == "bspl":
        c_bspl = np.zeros(deg + 1)
        knots = np.concatenate([np.zeros(deg + 1), np.ones(deg + 1)])
        for j in range(deg + 1):
            c_bspl[j] = poly2bspl_coeff(c_poly, j, knots)
        parameters["pres"] = {"type": "bspline", "knots": knots, "coefs": c_bspl}
        parameters["iota"] = {"type": "bspline", "knots": knots, "coefs": c_bspl}
    parameters["iota"]["scale"] = 1
    parameters["pres"]["scale"] = 1
    return State.new(parameters, tmp_path)


@pytest.fixture()
def teststate_min(teststate):
    return MinimizationState(teststate.parameterfile, teststate.statefile)


# === Tests === #
@pytest.mark.parametrize("type", ["p", "iota"])
def test_eval_profile(teststate, type, c_poly):
    rho = np.linspace(0, 1, 100)
    ref_poly = Polynomial(c_poly)
    ref_poly_vals = ref_poly(rho**2)
    gvec_profile = teststate.evaluate_profile(type, rho=rho)
    np.testing.assert_allclose(gvec_profile, ref_poly_vals)


@pytest.mark.parametrize("type", ["p", "iota"])
def test_eval_profile_prime(teststate, type, c_poly):
    rho = np.linspace(0, 1, 100)
    ref_poly = Polynomial(c_poly)
    ref_poly_vals = ref_poly.deriv(m=1)(rho**2) * 2 * rho
    gvec_profile = teststate.evaluate_profile(type, rho=rho, deriv=1)
    np.testing.assert_allclose(gvec_profile, ref_poly_vals)


@pytest.mark.parametrize("type", ["p", "iota"])
def test_eval_profile_n_deriv(teststate, type, c_poly):
    rho = np.linspace(0, 1, 100)
    ref_poly = Polynomial(c_poly)

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
    ref_poly = Polynomial(c_poly)
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
    ref_poly = Polynomial(c_poly)
    ref_poly_vals = ref_poly(rho2)
    eval_iota = teststate.evaluate_rho2_profile("iota", rho2=rho2)
    np.testing.assert_allclose(eval_iota, ref_poly_vals, atol=tol)
    eval_dphids = teststate.evaluate_rho2_profile("Phi", rho2=rho2, deriv=1)
    eval_dchids = teststate.evaluate_rho2_profile("chi", rho2=rho2, deriv=1)
    np.testing.assert_allclose(eval_dchids / eval_dphids, ref_poly_vals, atol=tol)


# NOTE: the get_iota_fobj method is not actually used and has issues with ifx in DEBUG mode
# def test_get_iota_fobj(teststate_min):
#     iota_fobj = teststate_min._get_iota_fobj()
#     # f90wrap apparently doesn't know which implementation the returned object is
#     assert isinstance(iota_fobj, gvec.lib.c_rProfile)
#     r = 0.1
#     v1 = iota_fobj.eval_at_rho(r)
#     v2 = teststate_min.evaluate_profile("iota", rho=r)
#     np.testing.assert_allclose(v1, v2, atol=1e-12)


def test_get_iota(teststate_min):
    derivs = [1, 2, 3]
    iota_profile = teststate_min.iota
    typ = teststate_min.parameters["iota"]["type"]
    assert typ in ["polynomial", "bspline"]
    if typ == "polynomial":
        assert isinstance(iota_profile, Polynomial)
        derivatives = [iota_profile.deriv(m=d) for d in derivs]
    elif typ == "bspline":
        assert isinstance(iota_profile, BSpline)
        derivatives = [iota_profile.derivative(nu=d) for d in derivs]

    rho2 = np.linspace(0, 1, 101)
    python_values = iota_profile(rho2)
    fortran_values = teststate_min.evaluate_rho2_profile("iota", rho2)
    np.testing.assert_allclose(python_values, fortran_values, atol=1e-16)
    for d, deriv in zip(derivs, derivatives):
        python_values = deriv(rho2)
        fortran_values = teststate_min.evaluate_rho2_profile("iota", rho2, deriv=d)
        np.testing.assert_allclose(python_values, fortran_values, atol=1e-16)


def test_set_iota(teststate_min):
    # set a new iota profile and check that it is updated correctly
    new_iota = Polynomial([1, 2, 3])
    teststate_min.iota = new_iota
    assert teststate_min.iota == new_iota

    rho2 = np.linspace(0, 1, 101)
    python_values = new_iota(rho2)
    fortran_values = teststate_min.evaluate_rho2_profile("iota", rho2)
    np.testing.assert_allclose(python_values, fortran_values, atol=1e-16)


@pytest.mark.parametrize("BC_type_axis", ["not_a_knot", "1st_deriv", "2nd_deriv"])
@pytest.mark.parametrize("BC_type_edge", ["not_a_knot", "1st_deriv", "2nd_deriv"])
@pytest.mark.parametrize("n_points", [4, 5, 6, 11])
def test_interpolation(parameters, c_poly, BC_type_axis, BC_type_edge, n_points, tmp_path):
    pres_scale = 1500
    cubic_poly_c = c_poly[:4]

    cubic_poly = Polynomial(cubic_poly_c)
    rho2_vals = np.linspace(0, 1, n_points)
    P_vals = cubic_poly(rho2_vals)
    iota_vals = cubic_poly(rho2_vals)

    parameters["pres"] = dict(
        type="interpolation",
        rho2=rho2_vals,
        vals=P_vals,
        scale=pres_scale,
        BC_type_axis=BC_type_axis,
        BC_type_edge=BC_type_edge,
    )
    parameters["iota"] = dict(
        type="interpolation",
        rho2=rho2_vals,
        vals=iota_vals,
        scale=1,
        BC_type_axis=BC_type_axis,
        BC_type_edge=BC_type_edge,
    )

    rho2 = np.linspace(0, 1, 200)
    P = pres_scale * cubic_poly(rho2)
    iota = cubic_poly(rho2)

    state = State.new(parameters, tmp_path)
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

    # compare fortran and python implementation of the interpolation
    iota_params = parameters["iota"]
    iota_profile = gvec.util.profile_from_parameters(iota_params)
    np.testing.assert_allclose(iota_profile(rho2), iota_interpol)


@pytest.mark.parametrize("testcase", ["w7x_from_vmec_initLA_F"])
def test_vmec_profile_init(woutfile, parameters, tmp_path):
    """test whether the (fortran) profile implementation works when reading the profiles from a VMEC wout file."""
    wout = xr.load_dataset(woutfile)
    phi = wout.phi
    P_vmec = wout.presf
    iota_vmec = -1 * wout.iotaf
    rho_wout = np.sqrt((phi - phi[0]) / (phi[-1] - phi[0]))  # Normalized toroidal flux

    state = State.new(parameters, tmp_path)

    P_gvec = state.evaluate_profile("p", rho_wout, deriv=0)
    dP_dr_gvec = state.evaluate_profile("p", 0, deriv=1)
    iota_gvec = state.evaluate_profile("iota", rho_wout, deriv=0)
    diota_dr_gvec = state.evaluate_profile("iota", 0, deriv=1)

    np.testing.assert_allclose(P_gvec, P_vmec)
    np.testing.assert_allclose(iota_gvec, iota_vmec)
    np.testing.assert_allclose(diota_dr_gvec, 0, atol=tol)
    np.testing.assert_allclose(dP_dr_gvec, 0, atol=tol)


@pytest.mark.parametrize("testcase", ["w7x_from_vmec_initLA_F"])
@pytest.mark.parametrize("profile_type", ["pres", "iota"])
def test_vmec_with_custom_profile(woutfile, parameters, c_poly, profile_type, tmp_path):
    wout = xr.load_dataset(woutfile)
    phi = wout.phi
    rho_wout = np.sqrt((phi - phi[0]) / (phi[-1] - phi[0]))  # Normalized toroidal flux

    parameters[profile_type] = {}
    if profile_type == "pres":
        var = "p"
        scale = 1
        parameters["init_with_profile_pressure"] = "T"
    else:
        var = "iota"
        scale = -1
        parameters["init_with_profile_iota"] = "T"
    parameters[profile_type]["scale"] = scale
    parameters[profile_type]["coefs"] = c_poly
    parameters[profile_type]["type"] = "polynomial"

    state = State.new(parameters, tmp_path)

    rho = np.linspace(0, 1, 100)
    ref_poly = Polynomial(c_poly)
    profile_true = scale * ref_poly(rho**2)

    profile_gvec = state.evaluate_profile(var, rho, deriv=0)
    profile_gvec_at_vmec = state.evaluate_profile(var, rho_wout, deriv=0)

    with np.testing.assert_raises(AssertionError):
        np.testing.assert_allclose(profile_gvec_at_vmec, wout[f"{profile_type}f"])
    np.testing.assert_allclose(profile_gvec, profile_true)
