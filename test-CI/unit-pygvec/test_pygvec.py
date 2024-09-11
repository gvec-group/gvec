import pytest

try:
    import numpy as np
    import xarray as xr
    import os
    from pathlib import Path
    import subprocess

    import gvec
    import helpers
except ImportError:
    pytest.skip("Import Error", allow_module_level=True)


# === Fixtures === #


@pytest.fixture(scope="session")
def testgroup():
    return "unit-pygvec"


@pytest.fixture(scope="session")
def testcase():
    return "ellipstell_lowres"


@pytest.fixture(scope="session")
def testcase_run(testgroup, testcaserundir, testcase, annotations, artifact_pages_path):
    """run the default testcase and store the output"""
    # assume pip-build: the gvec executable is in the PATH
    args = ["gvec", "parameter.ini"]
    # run gvec - adapted from test_all.test_run
    with helpers.chdir(testcaserundir):
        # skip if stdout & stderr are present and ok
        try:
            helpers.assert_empty_stderr()
            helpers.assert_stdout_finished(message="GVEC SUCESSFULLY FINISHED!")
        except (FileNotFoundError, AssertionError):
            pass
        else:
            return
        # run GVEC
        with open("stdout.txt", "w") as stdout:
            stdout.write(f"RUNNING: \n {args} \n")
        with open("stdout.txt", "a") as stdout, open("stderr.txt", "w") as stderr:
            subprocess.run(args, text=True, stdout=stdout, stderr=stderr)
        for filename in ["stdout", "stderr"]:
            if pages_rundir := os.environ.get("CASENAME"):
                pages_rundir = f"CIrun_{pages_rundir}"
            else:
                pages_rundir = "."
            annotations["gvec-output"].append(
                dict(
                    external_link=dict(
                        label=f"{testgroup}/{testcase}/{filename}",
                        url=f"{artifact_pages_path}/{pages_rundir}/{testgroup}/{testcase}/{filename}.txt",
                    )
                )
            )
        # check if GVEC was successful
        helpers.assert_empty_stderr()
        helpers.assert_stdout_finished(message="GVEC SUCESSFULLY FINISHED!")
    return


@pytest.fixture()
def testfiles(tmpdir, testcaserundir, testcase_run):
    """prepare the ellipstell parameters"""
    paramfile = "parameter.ini"
    statefile = "ELLIPSTELL_LOWRES_State_0000_00000000.dat"
    with helpers.chdir(tmpdir):
        os.symlink(testcaserundir / paramfile, paramfile)
        os.symlink(testcaserundir / statefile, statefile)
        yield (paramfile, statefile)


@pytest.fixture()
def teststate(testfiles):
    with gvec.post.State(*testfiles) as state:
        yield state


@pytest.fixture()
def evals_r(teststate):
    rho = np.linspace(0, 1, 6)
    ds = gvec.post.Evaluations(state=teststate, coords={"rho": rho})
    return ds


@pytest.fixture()
def evals_rtz(teststate):
    rho = np.linspace(0, 1, 6)
    theta = np.linspace(0, 2 * np.pi, 32, endpoint=False)
    zeta = np.linspace(0, 2 * np.pi, 10, endpoint=False)
    ds = gvec.post.Evaluations(
        state=teststate, coords={"rho": rho, "theta": theta, "zeta": zeta}
    )
    return ds


@pytest.fixture()
def evals_rtz_int(teststate):
    ds = gvec.post.Evaluations(
        state=teststate,
        rho="int",
        theta="int",
        zeta="int",
    )
    return ds


@pytest.fixture()
def evals_r_int_tz(teststate):
    ds = gvec.post.Evaluations(
        state=teststate,
        rho="int",
        theta=100,
        zeta=100,
    )
    return ds


# === Tests === #


def test_version():
    import pkg_resources

    assert isinstance(gvec.__version__, str)
    assert gvec.__version_tuple__ >= (0, 2, 1)
    assert pkg_resources.get_distribution("gvec").version == gvec.__version__


def test_state(testfiles):
    with gvec.post.State(*testfiles) as state:
        assert isinstance(state, gvec.post.State)


def test_state_args(testfiles):
    paramfile, statefile = testfiles

    with pytest.raises(FileNotFoundError):
        state = gvec.post.State(paramfile, "nonexistent.dat")

    with pytest.raises(FileNotFoundError):
        state = gvec.post.State("nonexistent.ini", statefile)


def test_state_explicit(testfiles):
    state = gvec.post.State(*testfiles)
    assert isinstance(state, gvec.post.State)
    assert state.initialized
    state.finalize()
    assert not state.initialized

    with pytest.raises(RuntimeError):
        state.evaluate_base_tens("X1", None, [0.5], [0.5], [0.5])


def test_state_twice(testfiles):
    # double context
    with gvec.post.State(*testfiles) as state:
        pass
    with gvec.post.State(*testfiles) as state:
        pass

    # double explicit
    state = gvec.post.State(*testfiles)
    state.finalize()
    state = gvec.post.State(*testfiles)
    state.finalize()

    # simultaneous context (not implemented yet)
    with pytest.raises(NotImplementedError):
        with gvec.post.State(*testfiles) as state1:
            with gvec.post.State(*testfiles) as state2:
                pass


@pytest.mark.parametrize("quantity", ["X1", "X2", "LA"])
def test_integration_points(teststate, quantity):
    r_GP, r_w, t_n, t_w, z_n, z_w = teststate.get_integration_points(quantity)
    assert np.isclose(t_w, 2 * np.pi / t_n)
    assert np.isclose(z_w, 2 * np.pi / z_n)
    assert np.isclose(np.sum(r_w), 1.0)
    assert np.isclose(np.sum(r_w * r_GP), 0.5)
    assert np.all(r_GP >= 0.0) and np.all(r_GP <= 1.0)
    assert np.all(np.diff(r_GP) > 0.0)


def test_evaluations_init(teststate):
    ds = gvec.post.Evaluations(state=teststate, coords={"rho": [0.5, 0.6]})
    assert np.allclose(ds.rho, [0.5, 0.6])
    assert {"rho"} == set(ds.coords)

    ds = gvec.post.Evaluations(teststate, rho=[0.5, 0.6])
    assert np.allclose(ds.rho, [0.5, 0.6])
    assert {"rho"} == set(ds.coords)

    ds = gvec.post.Evaluations(
        teststate, rho=np.array([0.5, 0.6]), theta="int", zeta="int"
    )
    assert np.allclose(ds.rho, [0.5, 0.6])
    assert {"rho", "theta", "zeta"} == set(ds.coords)
    assert ds.rho.attrs["integration_points"] == False
    ds.compute("iota_mean")

    ds = gvec.post.Evaluations(teststate, rho="int", theta=None, zeta=None)
    assert {"rho"} == set(ds.coords)
    assert ds.rho.attrs["integration_points"] == True
    ds.compute("iota_mean")
    ds.compute("V")

    ds = gvec.post.Evaluations(teststate, "int", "int", "int")
    for c in ["rho", "theta", "zeta"]:
        assert c in ds.coords
        assert ds[c].attrs["integration_points"] == True
    ds.compute("V")


def test_evaluate_base_tens(teststate):
    rho = np.linspace(0, 1, 6)
    theta = np.linspace(0, 2 * np.pi, 32, endpoint=False)
    zeta = np.linspace(0, 2 * np.pi, 10, endpoint=False)

    # base evaluation
    X1 = teststate.evaluate_base_tens("X1", None, rho, theta, zeta)
    X2 = teststate.evaluate_base_tens("X2", None, rho, theta, zeta)
    LA = teststate.evaluate_base_tens("LA", None, rho, theta, zeta)
    dX2_drtz = teststate.evaluate_base_tens("X2", "rtz", rho, theta, zeta)

    assert X1.shape == (6, 32, 10)
    assert not np.any(np.isnan(X1))
    assert not np.any(np.isnan(X2))
    assert not np.any(np.isnan(LA))
    # magnetic axis collapses to a line
    assert np.allclose(np.std(X1[0, ...], axis=0), 0.0)
    assert np.allclose(np.std(X2[0, ...], axis=0), 0.0)


def test_evaluate_base_tens_vectorize(teststate):
    """Test if the vectorization works properly (correct order of indices)"""
    X1_222 = teststate.evaluate_base_tens("X1", None, [0.5, 0.6], [0, 0.1], [0, 0.1])
    X1_122 = teststate.evaluate_base_tens("X1", None, [0.5], [0, 0.1], [0, 0.1])
    X1_212 = teststate.evaluate_base_tens("X1", None, [0.5, 0.6], [0.1], [0, 0.1])
    X1_221 = teststate.evaluate_base_tens("X1", None, [0.5, 0.6], [0, 0.1], [0.1])
    X1_112 = teststate.evaluate_base_tens("X1", None, [0.5], [0.1], [0, 0.1])
    assert np.allclose(X1_222[0, :, :], X1_122[0, :, :])
    assert np.allclose(X1_222[:, 1, :], X1_212[:, 0, :])
    assert np.allclose(X1_222[:, :, 1], X1_221[:, :, 0])
    assert np.allclose(X1_222[0, 1, :], X1_112[0, 0, :])


def test_evaluate_base_tens_bounds(teststate):
    rho = np.linspace(0, 1, 6)
    theta = np.linspace(0, 4 * np.pi, 8, endpoint=False)
    zeta = np.linspace(-2 * np.pi, 0, 10, endpoint=False)

    teststate.evaluate_base_tens("X1", None, rho, theta, zeta)
    with pytest.raises(ValueError):
        teststate.evaluate_base_tens("X3", None, rho, theta, zeta)
    rho = np.linspace(-1, 1, 6)
    with pytest.raises(ValueError):
        teststate.evaluate_base_tens("X1", None, rho, theta, zeta)
    rho = np.linspace(0, 2, 6)
    with pytest.raises(ValueError):
        teststate.evaluate_base_tens("X1", None, rho, theta, zeta)


def test_evaluate_base_tens_all(teststate):
    rho = np.linspace(0, 1, 6)
    theta = np.linspace(0, 2 * np.pi, 32, endpoint=False)
    zeta = np.linspace(0, 2 * np.pi, 10, endpoint=False)

    results = teststate.evaluate_base_tens_all("X1", rho, theta, zeta)
    assert len(results) == 10
    assert all(result.shape == (6, 32, 10) for result in results)
    with pytest.raises(ValueError):
        teststate.evaluate_base_tens_all("X3", rho, theta, zeta)


def test_evaluate_base_list_tz_all(teststate):
    rho = np.linspace(0, 1, 6)
    theta = np.linspace(0, 2 * np.pi, 32, endpoint=False)
    zeta = np.linspace(0, 2 * np.pi, 10, endpoint=False)
    T, Z = np.meshgrid(theta, zeta)
    thetazeta = np.stack([T.flatten(), Z.flatten()])

    results = teststate.evaluate_base_list_tz_all("X1", rho, thetazeta)
    assert len(results) == 10
    assert all(result.shape == (6, 32 * 10) for result in results)
    with pytest.raises(ValueError):
        results = teststate.evaluate_base_list_tz_all("X1", rho, thetazeta.T)


def test_evaluate_base_all_compare(teststate):
    rho = np.linspace(0, 1, 6)
    theta = np.linspace(0, 2 * np.pi, 32, endpoint=False)
    zeta = np.linspace(0, 2 * np.pi, 10, endpoint=False)
    T, Z = np.meshgrid(theta, zeta, indexing="ij")
    thetazeta = np.stack([T.flatten(), Z.flatten()])

    tens = teststate.evaluate_base_tens_all("X1", rho, theta, zeta)
    list_tz = teststate.evaluate_base_list_tz_all("X1", rho, thetazeta)
    for qt, ql in zip(tens, list_tz, strict=True):
        assert np.allclose(qt, ql.reshape(6, 32, 10))


def test_init_base_Boozer(teststate):
    rho = np.linspace(0, 1, 6)
    theta = np.linspace(0, 2 * np.pi, 32, endpoint=False)
    zeta = np.linspace(0, 2 * np.pi, 10, endpoint=False)

    assert not teststate.GB_available
    with pytest.raises(AttributeError):
        teststate.evaluate_base_tens_all("GB", rho, theta, zeta)

    teststate.construct_GB(degree=4, continuity=3, m=3, n=3, fourier="sin")
    assert teststate.GB_available
    teststate.evaluate_base_tens_all("GB", rho, theta, zeta)

    # teststate.project_GB(degree=5, continuity=4, m=3, n=3, fourier="both")

    with pytest.raises(ValueError):
        teststate.construct_GB(degree=5, continuity=3, m=3, n=3, fourier="both")


@pytest.mark.xfail(reason="GB is not accurate enough yet")
@pytest.mark.parametrize("method", ["projection", "interpolation"])
def test_compute_Boozer(teststate, evals_rtz_int, method):
    with pytest.raises(AttributeError):
        evals_rtz_int.compute("GB")
    teststate.construct_GB(degree=4, m=5, n=2, method=method)
    evals_rtz_int.compute("GB", "dGB_dt", "dGB_dz", "dGB_dt_def", "dGB_dz_def")
    assert "GB" in evals_rtz_int
    assert evals_rtz_int.GB.isel(rho=slice(1, None)).isnull().sum() == 0
    assert np.allclose(evals_rtz_int.dGB_dt, evals_rtz_int.dGB_dt_def)
    assert np.allclose(evals_rtz_int.dGB_dz, evals_rtz_int.dGB_dz_def)


def test_compute_base(evals_rtz):
    ds = evals_rtz

    ds.compute("X1")
    ds.compute("X2")
    ds.compute("LA")

    assert ds.X1.shape == (6, 32, 10)
    assert "dX1_dr" in ds
    assert "dX2_dt" in ds
    assert "dLA_dz" in ds
    assert "dX1_drt" in ds
    assert not np.any(np.isnan(ds.X1))


def test_compute_implicit(evals_rtz):
    ds = evals_rtz
    # multiple computes
    assert "X1" not in ds and "X2" not in ds
    ds.compute("X1", "X2")
    assert "X1" in ds and "X2" in ds
    # implict compute with __getitem__
    assert "LA" not in ds
    assert isinstance(ds["LA"], xr.DataArray)
    assert "LA" in ds
    # implicit compute with __getattr__
    assert "p" not in ds
    assert isinstance(ds.p, xr.DataArray)
    assert "p" in ds


def test_evaluate_hmap(teststate):
    rho = np.linspace(0, 1, 6)
    theta = np.linspace(0, 2 * np.pi, 32, endpoint=False)
    zeta = np.linspace(0, 2 * np.pi, 10, endpoint=False)
    R, T, Z = np.meshgrid(rho, theta, zeta, indexing="ij")
    # X1, X2, dX1_dr, dX2_dr, dX1_dt, dX2_dt, dX1_dz, dX2_dz
    X1 = teststate.evaluate_base_tens_all("X1", rho, theta, zeta)
    X2 = teststate.evaluate_base_tens_all("X2", rho, theta, zeta)
    inputs = sum([[x1, x2] for x1, x2 in zip(X1[:4], X2[:4])], [])
    inputs = inputs[:2] + [Z] + inputs[2:]
    inputs = [i.flatten() for i in inputs]

    outputs = teststate.evaluate_hmap(*inputs)
    assert len(outputs) == 4
    assert all(output.shape == (3, 6 * 32 * 10) for output in outputs)


def test_compute_hmap(evals_rtz):
    ds = evals_rtz
    ds.compute("pos")

    assert ds.pos.shape == (3, 6, 32, 10)
    assert "vector" in ds.dims
    assert "vector" in ds.coords
    assert "e_X1" in ds
    assert "e_X2" in ds
    assert "e_zeta3" in ds
    assert not np.any(np.isnan(ds.pos))

    ds.compute("e_rho")
    assert ds.e_rho.shape == (3, 6, 32, 10)


def test_compute_metric(evals_rtz):
    ds = evals_rtz

    ds.compute("g_tt")
    for ij in ("rr", "rt", "rz", "tt", "tz", "zz"):
        key = f"g_{ij}"
        assert key in ds
        assert set(ds[key].coords) == {"rho", "theta", "zeta"}
        for k in "rtz":
            assert f"dg_{ij}_d{k}" in ds

    ds.compute("e_rho")
    ds.compute("e_theta")
    ds.compute("e_zeta")
    idxs = {"r": "rho", "t": "theta", "z": "zeta"}
    for ij in "rr rt rz tt tz zz".split():
        key = f"g_{ij}"
        assert np.allclose(
            ds[f"g_{ij}"],
            xr.dot(ds[f"e_{idxs[ij[0]]}"], ds[f"e_{idxs[ij[1]]}"], dim="vector"),
        )


def test_evaluate_profile(teststate):
    rho = np.linspace(0, 1, 6)

    iota = teststate.evaluate_profile("iota", rho)
    assert iota.shape == rho.shape
    dp_dr = teststate.evaluate_profile("p_prime", rho)


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
def test_compute_profile(evals_r, quantity):
    ds = evals_r
    ds.compute(quantity)
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
def test_compute_derived(evals_rtz, quantity, on_axis):
    ds = evals_rtz
    ds.compute(quantity)
    assert quantity in ds
    if on_axis:
        assert not np.any(np.isnan(ds[quantity]))
    else:
        assert np.all(np.isnan(ds[quantity].isel(rho=0)))
        assert not np.any(np.isnan(ds[quantity].isel(rho=slice(1, None))))


def test_compute_basis(evals_rtz):
    ds = evals_rtz
    for coord in ["rho", "theta", "zeta"]:
        ds.compute(f"e_{coord}", f"grad_{coord}")
    ds = ds.isel(rho=slice(1, None))
    for coord in ["rho", "theta", "zeta"]:
        assert np.allclose(
            xr.dot(ds[f"e_{coord}"], ds[f"grad_{coord}"], dim="vector"), 1.0
        )
        for coord2 in ["rho", "theta", "zeta"]:
            if coord2 == coord:
                continue
            ds.compute(f"e_{coord2}", f"grad_{coord2}")
            assert np.allclose(
                xr.dot(ds[f"e_{coord}"], ds[f"grad_{coord2}"], dim="vector"), 0.0
            )
            assert np.allclose(
                xr.dot(ds[f"grad_{coord}"], ds[f"e_{coord2}"], dim="vector"), 0.0
            )


def test_volume_integral(evals_rtz_int, evals_rtz):
    ds = evals_rtz_int
    one = xr.ones_like(ds.Jac)
    V_default = ds.volume_integral(one)
    V_true = ds.volume_integral(one, Jac=True)
    V_false = ds.volume_integral(ds.Jac, Jac=False)
    V_str = ds.volume_integral("Jac", Jac=False)
    V_explicit = ds.volume_integral(one, Jac=ds.Jac)
    for V in [V_default, V_false, V_str, V_explicit]:
        assert np.allclose(V, V_true)
    assert np.allclose(4 * np.pi**2, ds.volume_integral(one, Jac=False))
    for Jac in [False, True, ds.Jac, ds.Jac_l, "Jac"]:
        assert np.allclose(1, ds.volume_average(one, Jac=Jac))

    with pytest.raises(ValueError):
        evals_rtz.volume_integral(one, Jac=ds.Jac_l)


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
    ],
)
def test_integral_quantities(evals_rtz_int, evals_r_int_tz, quantity):
    ds = evals_rtz_int
    ds.compute(quantity)
    # --- check that metadata is preserved --- #
    assert ds.rho.attrs["integration_points"] == True

    # automatic change to integration points if necessary
    ds = evals_r_int_tz
    ds.compute(quantity)
    assert ds.rho.attrs["integration_points"] == True
    assert ds.theta.attrs["integration_points"] == False
    assert ds.zeta.attrs["integration_points"] == False
    assert "theta_weight" not in ds and "zeta_weight" not in ds
    assert "rho_weights" in ds
