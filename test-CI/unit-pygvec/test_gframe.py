import pytest

try:
    import numpy as np
    import gvec
    from gvec import gframe
except ImportError:
    # in this case, we already need the gvec package during test collection
    # therefore we skip the whole module
    pytest.skip("Skipping test_gframe.py: gvec package not available", allow_module_level=True)


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
    gvec.util.boundary_generator_cases().keys(),
)
@pytest.mark.parametrize("nfp", [1, 2, 3])
@pytest.mark.parametrize("nz_nt", [[20, 10], [21, 11]], ids=["even", "odd"])
@pytest.mark.parametrize("shft", [0.0, 0.5], ids=["noshift", "halfshift"])
def test_gframe_boundary(case, nz_nt, shft, nfp):
    def eval_bnd_xyz(theta_in, zeta_in, bnd_params):
        X1, X2 = gvec.util.evaluate_boundary(theta_in, zeta_in, bnd_params)
        # simply use RZ hmap here:
        xyz = np.zeros((len(zeta_in), len(theta_in), 3))
        xyz[:, :, 0] = X1.T * np.cos(zeta_in[:, None])
        xyz[:, :, 1] = -X1.T * np.sin(zeta_in[:, None])
        xyz[:, :, 2] = X2.T
        return xyz

    tol = 1e-8
    bnd_params = gvec.util.boundary_generator(case)
    M_in, N_in = bnd_params["X1_mn_max"]
    bnd_params["nfp"] = nfp
    ntheta = 11
    nzeta = 21
    theta = np.linspace(0, 2 * np.pi, ntheta, endpoint=False)
    zeta = np.linspace(0, 2 * np.pi / nfp, nzeta, endpoint=False)

    X1_fp, X2_fp = gvec.util.evaluate_boundary(theta, zeta, bnd_params)
    M, N = gframe.minimal_modes(X1_fp, X2_fp, tolerance=tol)
    assert M <= M_in and N <= N_in, (
        f"boundary surface modes (M={M},N={N}) are larger than they should be (M={M_in},N={N_in})"
    )
    zetafull = np.linspace(0, 2 * np.pi, nzeta * nfp, endpoint=False)
    # generate data to compare with
    xyz_target = eval_bnd_xyz(theta, zetafull, bnd_params)
    # test odd/even input points and shifts
    t = np.linspace(0, 2 * np.pi, nz_nt[1], endpoint=False)
    z = np.linspace(0, 2 * np.pi, nz_nt[0] * nfp, endpoint=False)
    t = t + shft * t[1]
    z = z + shft * z[1]

    xyz_in = eval_bnd_xyz(t, z, bnd_params)

    params, dict_out = gframe.construct_gframe_from_surface(
        xyz_in,
        nfp,
        theta0=t[0],
        zeta0=z[0],
        name=f"test_{case}",
        writeFiles=False,
        tolerance_clean_surface=tol,
    )
    dict_surf = gframe.to_surface(dict_out, ntheta=ntheta, nzeta=nzeta, tolerance=tol)

    M, N = params["X1_mn_max"]
    dist = np.max(np.sqrt(np.sum((xyz_target - dict_surf["xyz"]) ** 2, axis=-1)))
    assert dist < tol, (
        f"xyz from gframe boundary surface does not match the input surface,max|sqrt((xvec-xvec_in)^2)||={dist}"
    )

    assert M <= M_in and N <= N_in, (
        f"gframe boundary surface '{case}' X1,X2 (M={M},N={N}) is larger than the input surface (M={M_in},N={N_in})"
    )


@pytest.mark.parametrize(
    "case",
    [
        "ellip_cyl",
        "ellip_cyl_rot",
        "ellip_cyl_rot2",
        "ellip_cyl_helix",
        "ellip_cyl_helix_rot2",
        "ellip_cyl_helix_rot",
    ],
)
@pytest.mark.parametrize("nfp", [1, 2, 3])
def test_writhe_boundary(case, nfp):
    params = gvec.util.boundary_generator(case, X1_00=1.0 + 0.25 * nfp)
    params["nfp"] = nfp
    theta = np.linspace(0, 2 * np.pi, 10, endpoint=False)

    zeta = np.linspace(0, 2 * np.pi, nfp * 200, endpoint=False)
    X1, X2 = gvec.util.evaluate_boundary(theta, zeta, params)
    xyz_surf = np.zeros((len(zeta), len(theta), 3))
    xyz_surf[:, :, 0] = X1.T * np.cos(zeta[:, None])
    xyz_surf[:, :, 1] = X1.T * np.sin(zeta[:, None])
    xyz_surf[:, :, 2] = X2.T

    # X0 = 0.5*(xyz_surf[:, len(theta) // 2, :]+ xyz_surf[:, 0, :]) # mid-curve
    C_a = xyz_surf[:, 0, :]
    C_b = xyz_surf[:, len(theta) // 2, :]

    Lk = gvec.util.linking_number(C_a, C_b, endpoint=False)
    # Nvec = Nvec/np.sqrt(np.sum(Nvec**2,axis=-1))[:,None]
    # "exact" writhe from linking - twist
    Tw = gvec.gframe.twist_of_ribbon(C_a, C_b - C_a)
    # approximate writhe
    Wr = gvec.util.writhe(C_a, endpoint=False)
    assert np.abs(Wr - (Lk - Tw)) < 1e-3, (
        f"Compare approx. writhe {Wr:5.2e} with (Lk - Tw). Lk = {Lk:.0f} Tw={Tw:5.2e} ... Difference: {Wr - (Lk - Tw):5.2e}"
    )
