import pytest
from pathlib import Path

try:
    import numpy as np
    import gvec
    from gvec import gframe
except ImportError:
    # in this case, we already need the gvec package during test collection
    # therefore we skip the whole module
    pytest.skip("Skipping test_gframe.py: gvec package not available", allow_module_level=True)


DATA = Path(__file__).parent / "../data"

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
    # "exact" writhe from linking - twist
    Tw = gvec.gframe.twist_of_ribbon(C_a, C_b - C_a)
    # approximate writhe
    Wr_polygon = gvec.util.writhe_from_polygon(C_a, endpoint=False)
    assert np.abs(Wr_polygon - (Lk - Tw)) < 1e-3, (
        f"Compare writhe from polygon {Wr_polygon:5.2e} with writhe from riboon, (Lk - Tw). Lk = {Lk:.0f} Tw={Tw:5.2e} ... Difference: {Wr_polygon - (Lk - Tw):5.2e}"
    )


@pytest.mark.parametrize(
    "gframe_ncfile",
    [
        "2022_QH_nfp7.nc",
        "brezellator-Gframe.nc",
        "brezellator-stellsym-scaled-Gframe.nc",
        "N2.12-v3.1-hi.nc",
        "N2.12-v3.1-lo.nc",
        "N2-63.nc",
        "N3.471-v3.1-lo.nc",
        "N3.Knot-v3.1-lo.nc",
    ],
)
def test_gframe_files_writhe(gframe_ncfile):
    gframe_ncfile = DATA / gframe_ncfile
    dict_gframe = gvec.gframe.read_Gframe_ncfile(gframe_ncfile)
    assert "axis" in dict_gframe and "boundary" in dict_gframe, (
        f"Expected 'axis' and 'boundary' in gframe dict, but not found. Keys: {dict_gframe.keys()}"
    )
    dict_highres = gvec.gframe.to_axis(dict_gframe, nzeta=201)
    xyz = dict_highres["axis"]["xyz"].T
    Nxyz = dict_highres["axis"]["Nxyz"].T
    Wr_polygon = gvec.util.writhe_from_polygon(xyz)
    Wr_gframe = gvec.gframe.writhe(xyz, N=Nxyz)[0]
    assert np.isclose(Wr_polygon, Wr_gframe, atol=1e-3), (
        f"Writhe from polygon {Wr_polygon} and from gframe {Wr_gframe} differ by more than 1e-3 for file {gframe_ncfile}"
    )


def trefoil(t, scale=3.0, zscale=1.0):
    x = np.sin(t) + 2 * np.sin(2 * t)
    y = np.cos(t) - 2 * np.cos(2 * t)
    z = -np.sin(3 * t) * zscale
    X0 = np.stack((x, y, z), axis=1)

    xp = np.cos(t) + 4 * np.cos(2 * t)
    yp = -np.sin(t) + 4 * np.sin(2 * t)
    zp = -3 * np.cos(3 * t) * zscale
    X0p = np.stack((xp, yp, zp), axis=1)

    xpp = -np.sin(t) - 8 * np.sin(2 * t)
    ypp = -np.cos(t) + 8 * np.cos(2 * t)
    zpp = 9 * np.sin(3 * t) * zscale
    X0pp = np.stack((xpp, ypp, zpp), axis=1)

    xppp = -np.cos(t) - 16 * np.cos(2 * t)
    yppp = np.sin(t) - 16 * np.sin(2 * t)
    zppp = 27 * np.cos(3 * t) * zscale
    X0ppp = np.stack((xppp, yppp, zppp), axis=1)

    return scale * X0, scale * X0p, scale * X0pp, scale * X0ppp


def get_frenet_trefoil(nz, scale=3.0, zscale=1.0):
    t = np.linspace(0, 2 * np.pi, nz, endpoint=False)
    X0, X0p, X0pp, X0ppp = trefoil(t, scale=scale, zscale=zscale)
    return gframe.frenet_frame_evaluate(X0, X0p, X0pp, X0ppp)


def test_frenet_trefoil():
    dict_frenet_exact = get_frenet_trefoil(3 * 51)
    dict_frenet = gframe.frenet_frame(dict_frenet_exact["X0"])
    for key in dict_frenet_exact.keys():
        assert key in dict_frenet, (
            f"Expected key '{key}' in frenet dict, but not found. Keys: {dict_frenet.keys()}"
        )
        if dict_frenet[key] is None:
            assert dict_frenet_exact[key] is None, (
                f"Expected {key} to be None, but got {dict_frenet_exact[key]}"
            )
        else:
            np.testing.assert_allclose(
                dict_frenet[key],
                dict_frenet_exact[key],
                rtol=1e-5,
                atol=1e-8,
                err_msg=f"Mismatch in {key}",
            )


def test_frenet_zero_curvature():
    # test that for a closed curve with zero curvature, the frenet frame detects it.
    zeta = np.linspace(0, 2 * np.pi, 2 * 11, endpoint=False)
    x = np.sin(zeta)
    y = np.sin(zeta) * np.cos(zeta)
    z = np.sin(zeta)
    X0 = np.stack([x, y, z], axis=-1)
    dict_frenet = gvec.gframe.frenet_frame(X0)
    assert (
        dict_frenet["N"] is None and dict_frenet["B"] is None and dict_frenet["tau"] is None
    ), "Expected N,B,tau to be None for zero curvature"


@pytest.mark.parametrize("zscale", [-1.0, -0.45, -0.1, 0.05, 0.1, 0.5, 1.0])
def test_writhe_trefoil(zscale):
    dict_frenet_high = get_frenet_trefoil(3 * 81, zscale=zscale)
    Wr_high, Lk_high, Tw_high = gvec.gframe.writhe(
        dict_frenet_high["X0"], N=dict_frenet_high["N"], nint=481
    )
    nzetas = [3 * 10, 3 * 21]
    Wr_from_polygon = np.zeros_like(nzetas, dtype=float)
    Wr_frenet = np.zeros_like(nzetas, dtype=float)
    Wr_centroid = np.zeros_like(nzetas, dtype=float)
    for i, nz in enumerate(nzetas):
        dict_frenet_x10 = get_frenet_trefoil(nz * 10, zscale=zscale)
        Wr_from_polygon[i] = gvec.util.writhe_from_polygon(dict_frenet_x10["X0"])
        dict_frenet = get_frenet_trefoil(nz, zscale=zscale)
        Wr_centroid[i], _, _ = gvec.gframe.writhe(dict_frenet["X0"])
        Wr_frenet[i], _, _ = gvec.gframe.writhe(dict_frenet["X0"], N=dict_frenet["N"])

    for method, Wr in zip(
        ["centroid", "frenet", "polygon"], [Wr_centroid, Wr_frenet, Wr_from_polygon]
    ):
        Wr_diff = np.abs(Wr - Wr_high)
        assert Wr_diff[0] > Wr_diff[1], (
            f"Writhe from {method} for trefoil not converging: |diff(Wr) = {Wr_diff[0]} -> {Wr_diff[1]} with reference high-res Wr={Wr_high}"
        )

    np.testing.assert_allclose(
        Wr_frenet,
        Wr_high,
        rtol=1e-5,
        atol=1e-8,
        err_msg="writhe from frenet frame of trefoil not accurate <1e-8",
    )
    np.testing.assert_allclose(
        Wr_centroid,
        Wr_high,
        rtol=1e-5,
        atol=1e-6,
        err_msg="writhe from centroid frame of trefoil not accurate <1e-6",
    )

    np.testing.assert_allclose(
        Wr_high,
        Wr_from_polygon,
        rtol=1e-5,
        atol=1e-3,
        err_msg="writhe from polygon for trefoil not within 1e-3",
    )


def tilted_elliptic_torus(
    ntheta, nzeta, R0=2.5, a=0.5, b=0.25, shift=[0.1, -0.5, 0.3], angle=np.pi / 3
):
    """
    Generate a torus with elliptic cross-section, shifted and tilted.
    """
    theta = np.linspace(0, 2 * np.pi, ntheta, endpoint=False)
    zeta = np.linspace(0, 2 * np.pi, nzeta, endpoint=False)
    x = (R0 + a * np.cos(theta[None, :])) * np.cos(zeta[:, None])
    y = -(R0 + a * np.cos(theta[None, :])) * np.sin(zeta[:, None])
    z = b * np.sin(theta[None, :]) + 0 * zeta[:, None]
    xyz = np.stack((x, y, z), axis=-1)
    # shift xyz and tilt the torus around x axis by pi/3
    R = np.array(
        [[1, 0, 0], [0, np.cos(angle), -np.sin(angle)], [0, np.sin(angle), np.cos(angle)]]
    )
    xyz_tilt = (xyz - np.array(shift)) @ R.T
    return xyz_tilt


@pytest.mark.parametrize("ntheta", [41, 44])
@pytest.mark.parametrize("nzeta", [81, 86])
def test_gframe_surface_volume(ntheta, nzeta):
    """Test that the volume computed from the surface matches the expected volume for a known torus with elliptic cross-section."""
    # Parameters for the torus
    R0 = 2.5  # major radius
    a = 0.5  # minor radius in R-direction
    b = 0.25  # minor radius in Z-direction
    xyz_tilt = tilted_elliptic_torus(ntheta, nzeta, R0=R0, a=a, b=b)
    # compute the volume from the surface
    volume_computed = gvec.gframe.surface_volume(xyz_tilt)
    # compute the expected volume analytically:
    volume_expected = 2 * np.pi**2 * R0 * a * b
    assert np.isclose(volume_computed, volume_expected, rtol=1e-8, atol=1e-12), (
        f"Computed volume {volume_computed} does not match expected volume {volume_expected} for the torus with elliptic cross-section"
    )
    xyz_mirror = xyz_tilt.copy()
    xyz_mirror[:, :, 2] = -xyz_mirror[:, :, 2]
    volume_computed_mirror = gvec.gframe.surface_volume(xyz_mirror)
    assert np.isclose(volume_computed_mirror, -volume_computed, rtol=1e-8, atol=1e-12), (
        f"mirrored surface volume {volume_computed_mirror} does not match unmirrored surface volume {-volume_computed} for a torus with elliptic cross-section"
    )


def test_gframe_surface_orientation():
    """Test the construction from surface for a tilted elliptic torus, and check that it raises a ValueError for the mirrored surface."""
    xyz_tilt = tilted_elliptic_torus(ntheta=41, nzeta=81)
    # this should work without error:
    params, dict_out = gframe.construct_gframe_from_surface(
        xyz_tilt, nfp=1, name="test_tilted_elliptic_torus"
    )
    # now mirror the surface and check that it raises a ValueError due to negative volume:
    xyz_mirror = xyz_tilt.copy()
    xyz_mirror[:, :, 2] = -xyz_mirror[:, :, 2]
    with pytest.raises(ValueError, match="negative volume"):
        params, dict_out = gframe.construct_gframe_from_surface(
            xyz_mirror, nfp=1, name="test_tilted_elliptic_torus_mirror"
        )


@pytest.mark.parametrize(
    "errormsg",
    [
        ("pm_jac_h", "positive and negative"),
        ("neg_jac_h", "negative everywhere"),
        ("neg_jac", "ALL JACOBIANS NEGATIVE"),
        ("pm_jac", "JACOBIAN WITH SIGN CHANGE"),
    ],
    ids=["pm_jac_h", "neg_jac_h", "neg_jac", "pm_jac"],
)
def test_gframe_gvec_run_detects_negative_jacobian(errormsg, tmp_path):
    """
    Test that GVEC run fails in the jacobian check either when the gframe is wrong, or the boundary is oriented wrong, or the axis is initially set wrong, and check that the expected error message is raised in each case.
    """
    nzeta = 5
    zetafull = np.linspace(0, 2 * np.pi, nzeta, endpoint=False)

    xyz = np.zeros((3, zetafull.shape[0]))
    Nxyz = np.zeros((3, zetafull.shape[0]))
    Bxyz = np.zeros((3, zetafull.shape[0]))

    # circular curve
    xyz[0, :] = 5.0 * np.cos(zetafull)
    xyz[1, :] = -5.0 * np.sin(zetafull)

    # First vector N is outward pointing
    Nxyz[0, :] = np.cos(zetafull)
    Nxyz[1, :] = -np.sin(zetafull)

    # Second vector B chosen to produce different Jh behavior:
    match errormsg[0]:
        case "pm_jac_h":
            Bxyz[2, :] = 1.0 * np.sin(zetafull)  # +/- jac_h
            X1_a_cos_00 = 0.0
            X2_b_sin_10 = 0.25
        case "neg_jac_h":
            Bxyz[2, :] = -1.0  # jac_h<0 everywhere
            X1_a_cos_00 = 0.0
            X2_b_sin_10 = (
                -0.25
            )  # that would restore Jac>0, but hmap is checked first, so it should still fail.
        case "neg_jac":
            Bxyz[2, :] = 1.0  # jac_h>0
            X1_a_cos_00 = 0.0
            X2_b_sin_10 = -0.25  # jac<0 wrong orientation
        case "pm_jac":
            Bxyz[2, :] = 1.0  # jac_h>0
            X1_a_cos_00 = -0.4  # +/- jac wrong axis.
            X2_b_sin_10 = 0.25

    # Create dict_gframe
    dict_gframe = {
        "nfp": 1,
        "axis": {
            "nzeta": nzeta,
            "zetafull": zetafull,
            "xyz": xyz,
            "Nxyz": Nxyz,
            "Bxyz": Bxyz,
        },
    }

    params = dict(
        ProjectName="test_" + errormsg[0],
        which_hmap=21,
        hmap_ncfile=f"test-Gframe{errormsg[0]}.nc",
        X1X2_deg=5,
        LA_deg=5,
        sgrid=dict(
            grid_type=0,
            nElems=5,
        ),
        X1_mn_max=(1, 0),
        X2_mn_max=(1, 0),
        LA_mn_max=(1, 0),
        minimize_tol=1e-3,
        totalIter=10,
        pres=dict(type="polynomial", coefs=[0.0]),
        iota=dict(type="polynomial", coefs=[0.7]),
        X1_a_cos={(0, 0): X1_a_cos_00},
        X1_b_cos={(1, 0): 0.5},
        X2_b_sin={(1, 0): X2_b_sin_10},
    )
    with gvec.util.chdir(tmp_path):
        gvec.gframe.write_Gframe_ncfile(f"test-Gframe{errormsg[0]}.nc", dict_gframe)
        # now try to run GVEC with this gframe, and check that it fails in the hmap jacobian check:
        with pytest.raises(gvec.errors.InitializationError, match=errormsg[1]):
            gvec.run(params)
