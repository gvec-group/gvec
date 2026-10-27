"""Tests related to MinimizationState"""

try:
    import numpy as np
    from gvec.core.minimize import MinimizationState
except ImportError:
    pass  # tests will be skipped via the `check_import` fixture


def test_circular_torus_current_on_axis(tmp_path):
    R0, minor, flux, current = 4.0, 0.5, 1.0, 1e5
    parameters = {
        "ProjectName": "test_circular_torus_current_on_axis",
        "which_hmap": 1,
        "nfp": 1,
        "PhiEdge": flux,
        "iota": {"type": "polynomial", "coefs": [0.25]},
        "pres": {"type": "polynomial", "coefs": [0.0]},
        "I_tor": {"type": "polynomial", "coefs": [0.0, current]},
        "sgrid": {"nElems": 4},
        "X1X2_deg": 5,
        "LA_deg": 5,
        "X1_mn_max": [1, 0],
        "X2_mn_max": [1, 0],
        "LA_mn_max": [1, 0],
        "X1_sin_cos": "_cos_",
        "X2_sin_cos": "_sin_",
        "LA_sin_cos": "_sin_",
        "X1_b_cos": {(0, 0): R0, (1, 0): minor},
        "X2_b_sin": {(1, 0): minor},
        "X1_a_cos": {(0, 0): R0},
        "init_LA": False,
        "totalIter": 1,
        "minimize_tol": 1e-6,
        "picard_current": {"iota_tol": 1e-7},
    }
    state = MinimizationState.new(parameters, rundir=tmp_path)
    state.knots_per_sgrid_element = 2
    fitted, *_ = state._get_new_iota(iota_tol=1e-7)
    # Ampere's law on a circular torus, lambda=0 and I(s)=current*s:
    # iota(s)=mu0*current*sqrt(R0**2-minor**2*s)/(2*flux).
    s = np.array([0.0, 1e-4, 0.1, 0.5, 1.0])
    expected = 4 * np.pi * 1e-7 * current * np.sqrt(R0**2 - minor**2 * s) / (2 * flux)
    np.testing.assert_allclose(fitted(s), expected, rtol=0, atol=1e-6)
