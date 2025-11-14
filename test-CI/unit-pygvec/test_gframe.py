import pytest

try:
    import numpy as np
    from gvec import util, gframe, fourier
except ImportError:
    pass  # tests will be skipped via the `check_import` fixture


# === Fixtures === #


@pytest.fixture(autouse=True)
def mark_xfail(request):
    if "test_gframe" in request.node.name and "_helix" in request.node.name:
        request.node.add_marker(
            pytest.mark.xfail(reason="_helix tests cases do not yet work with gframe script")
        )


# === Test GFRAME === #


@pytest.mark.parametrize(
    "case",
    util.boundary_generator_cases().keys(),
)
@pytest.mark.parametrize("nfp", [1, 2, 3])
def test_gframe_boundary(case, nfp):
    tol = 1e-8
    bnd_params = util.boundary_generator(case)
    bnd_params["nfp"] = nfp
    ntheta = 11
    nzeta = 21
    theta = np.linspace(0, 2 * np.pi, ntheta, endpoint=False)
    zeta = np.linspace(0, 2 * np.pi / nfp, nzeta, endpoint=False)
    M_in, N_in = bnd_params["X1_mn_max"]
    X1_fp, X2_fp = util.evaluate_boundary(theta, zeta, bnd_params)
    M, N = gframe.minimal_modes(X1_fp, X2_fp, tolerance=tol)
    assert M <= M_in and N <= N_in, (
        f"boundary surface modes (M={M},N={N}) are larger than they should be (M={M_in},N={N_in})"
    )
    zetafull = np.linspace(0, 2 * np.pi, nzeta * nfp, endpoint=False)
    X1, X2 = util.evaluate_boundary(theta, zetafull, bnd_params)
    # simply use RZ hmap here:
    xyz = np.zeros((nzeta * nfp, ntheta, 3))
    xyz[:, :, 0] = X1.T * np.cos(zetafull[:, None])
    xyz[:, :, 1] = -X1.T * np.sin(zetafull[:, None])
    xyz[:, :, 2] = X2.T

    params, dict_out = gframe.construct_gframe_from_surface(
        xyz, nfp, name=f"test_{case}", writeFiles=False, tolerance_clean_surface=tol
    )
    dict_surf = gframe.to_surface(dict_out, ntheta=ntheta, nzeta=nzeta, tolerance=tol)
    M, N = params["X1_mn_max"]
    dist = np.max(np.sqrt(np.sum((xyz - dict_surf["xyz"]) ** 2, axis=-1)))
    assert dist < tol, (
        f"xyz from gframe boundary surface does not match the input surface,max|sqrt((xvec-xvec_in)^2)||={dist}"
    )

    assert M <= M_in and N <= N_in, (
        f"gframe boundary surface '{case}' X1,X2 (M={M},N={N}) is larger than the input surface (M={M_in},N={N_in})"
    )
