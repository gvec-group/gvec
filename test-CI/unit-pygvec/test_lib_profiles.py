import pytest
from abc import ABC, abstractmethod

try:
    import numpy as np

    import gvec.lib
except ImportError:
    pytest.skip("Import Error", allow_module_level=True)


tol = 1e-12
rho = np.linspace(0, 1, 21)
poly_coefs = [1, 2, 3, -1, -2, -3]  # reference polynomial coefficients to test against


# === Utility functions === #
# for transforming polynomials into B-Splines.


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


@pytest.fixture(scope="module", autouse=True)
def redirect_abort():
    """Redirect abort to raise an exception instead."""
    gvec.lib.modgvec_py_binding.redirect_abort()


# === Tests === #


class BaseTestProfile(ABC):
    @abstractmethod
    @pytest.fixture()
    def profile(self):
        """the profile object (c_rProfile) to test"""
        pass

    @pytest.fixture()
    def reference_rho2(self, deriv):
        polynomial = np.polynomial.Polynomial(poly_coefs).deriv(deriv)
        return polynomial(rho**2)

    @pytest.fixture()
    def reference_rho(self, deriv):
        coefs = [poly_coefs[i // 2] if i % 2 == 0 else 0 for i in range(len(poly_coefs) * 2)]
        polynomial = np.polynomial.Polynomial(coefs).deriv(deriv)
        return polynomial(rho)

    def test_init(self, profile):
        assert isinstance(profile, gvec.lib.modgvec_rprofile_base.c_rProfile)
        assert profile._alloc

    @pytest.mark.parametrize("deriv", [0])
    def test_eval_rho(self, profile, reference_rho):
        values = [profile.eval_at_rho(r) for r in rho]
        np.testing.assert_allclose(values, reference_rho, atol=tol)

    @pytest.mark.parametrize("deriv", [0])
    def test_eval_rho2(self, profile, reference_rho2):
        values = [profile.eval_at_rho2(r**2) for r in rho]
        np.testing.assert_allclose(values, reference_rho2, atol=tol)

    @pytest.mark.parametrize("deriv", range(5))
    def test_eval_rho_deriv(self, profile, deriv, reference_rho):
        values = [profile.eval_at_rho(r, deriv) for r in rho]
        np.testing.assert_allclose(values, reference_rho, atol=tol)

    @pytest.mark.parametrize("deriv", [5, 10])
    def test_eval_rho_deriv_error(self, profile, deriv, reference_rho):
        with pytest.raises(RuntimeError):
            profile.eval_at_rho(rho[0], deriv)

    @pytest.mark.parametrize("deriv", range(len(poly_coefs)))
    def test_eval_rho2_deriv(self, profile, deriv, reference_rho2):
        values = [profile.eval_at_rho2(r**2, deriv) for r in rho]
        np.testing.assert_allclose(values, reference_rho2, atol=tol)


class TestProfilePoly(BaseTestProfile):
    """Test the polynomial profile."""

    @pytest.fixture()
    def profile(self):
        return gvec.lib.modgvec_rprofile_poly.t_rProfile_poly(poly_coefs)

    def test_init(self, profile):
        super().test_init(profile)
        assert profile.deg == len(poly_coefs) - 1


class TestProfileBSpline(BaseTestProfile):
    """Test the B-spline profile."""

    @pytest.fixture()
    def profile(self):
        """B-spline constructed from polynomial coefficients"""
        deg = len(poly_coefs) - 1
        knots = np.concatenate([np.zeros(deg + 1), np.ones(deg + 1)])
        bspl_coefs = np.array([poly2bspl_coeff(poly_coefs, j, knots) for j in range(deg + 1)])

        return gvec.lib.modgvec_rprofile_bspl.t_rProfile_bspl(knots, bspl_coefs)

    def test_init(self, profile):
        super().test_init(profile)
        assert profile.deg == len(poly_coefs) - 1
