import pytest

try:
    from gvec import util
    from gvec.vtk import ev2vtk
    import gvec
    import numpy as np
    from copy import deepcopy
except ImportError:
    pass  # tests will be skipped via the `check_import` fixture

# === TESTS === #


def test_CaseInsensitiveDict():
    cid = util.CaseInsensitiveDict({"a": 1, "B": 2, 3: "C"})
    assert cid["A"] == 1
    assert cid["b"] == 2
    assert cid[3] == "C"
    assert list(cid.keys()) == ["a", "B", 3]
    assert list(cid.lower_keys()) == ["a", "b", 3]
    assert list(cid.items()) == [("a", 1), ("B", 2), (3, "C")]
    assert list(cid.lower_items()) == [("a", 1), ("b", 2), (3, "C")]
    assert cid == util.CaseInsensitiveDict(cid.lower_items())
    with pytest.raises(KeyError):
        cid["3"]
    with pytest.raises(ValueError):
        cid.update({"d": 1, "D": 2})
    with pytest.raises(ValueError):
        _ = util.CaseInsensitiveDict({"a": 1, "A": 2})
    with pytest.raises(ValueError):
        cid == {"a": 1, "A": 2}
    cid.update(e=5)
    assert "f" in cid | {"f": 6}
    cid |= {"g": 7}
    assert "a" in cid
    assert "A" in cid
    assert "b" in cid
    assert "e" in cid
    assert "f" not in cid
    assert "g" in cid


def test_copy_test_CaseInsensitiveDict():
    from collections.abc import MutableSequence, MutableMapping

    def recursive_is_not(data_left, data_right):
        if not isinstance(data_left, int) and not isinstance(data_left, str):
            assert data_left is not data_right
        assert data_left == data_right
        match data_left:
            case MutableSequence():
                for values in zip(data_left, data_right):
                    recursive_is_not(*values)
            case MutableMapping():
                for key in data_left:
                    recursive_is_not(data_left[key], data_right[key])
            case _:
                pass

    cid = util.CaseInsensitiveDict(
        cid_in=util.CaseInsensitiveDict({"a": 1, "B": 2}),
        list_cid_in=[
            util.CaseInsensitiveDict({"a": 1, "B": 2}),
            util.CaseInsensitiveDict({"c": 3, "D": 4}),
        ],
        NAME="test",
    )
    cid_copy = deepcopy(cid)
    recursive_is_not(cid_copy, cid)


def test_ev2vtk(testcaserundir, testfiles):
    rho = np.linspace(0, 1, 4)
    rho[0] = 1e-4  # avoid evaluation at rho=0
    theta = np.linspace(0, 2 * np.pi, 5)  # including endpoints
    zeta = np.linspace(0, 2 * np.pi, 6)  # full torus, including endpoints

    vars_out = ["X1", "X2", "LA", "iota", "p", "pos", "grad_zeta", "e_rho", "B"]
    state = gvec.State(*testfiles)
    ev = gvec.Evaluations(rho=rho, theta=theta, zeta=zeta, state=state)
    state.compute(ev, *vars_out)

    ev2vtk(testcaserundir / "test_ev2.vtk", ev[vars_out])
    ev2vtk(testcaserundir / "test_ev2_only_profile.vtk", ev[["pos", "iota"]])
    ev2vtk(testcaserundir / "test_ev2_only_scalar.vtk", ev[["pos", "X1"]])
    ev2vtk(testcaserundir / "test_ev2_only_vector.vtk", ev[["pos", "B"]])


def test_read_write_parameters_ini(testcaserundir, testfiles, tmp_path):
    """
    Test reading and writing parameters with an INI file
    """
    # Read parameters from a file (ini)
    params = util.read_parameter_file_ini(testfiles[0])
    assert isinstance(params, util.CaseInsensitiveDict)
    assert "projectname" in params

    # Write parameters to a new file (ini)
    new_params_file = tmp_path / "parameters-copy.ini"
    util.write_parameter_file_ini(params, new_params_file)
    assert new_params_file.exists()

    # Read the new parameters file
    new_params = util.read_parameter_file_ini(new_params_file)
    assert new_params == params


@pytest.mark.parametrize(
    "func, transform",
    [
        ("theta", lambda t, z: (-t, z)),
        ("zeta", lambda t, z: (t, -z)),
    ],
    ids=["theta", "zeta"],
)
def test_flip_parameters(func, transform):
    func = getattr(util, f"flip_boundary_{func}")
    # Test flipping parameters in theta direction
    parameters = {
        "X1_b_cos": {
            (0, 0): 1.0,
            (0, 1): 1.1,
            (1, 0): 1.2,
            (1, 1): 1.3,
        },
        "X2_b_cos": {
            (0, 0): 2.0,
            (0, 1): 2.1,
            (1, 0): 2.2,
            (1, 1): 2.3,
        },
        "X1_b_sin": {
            (0, 1): 3.0,
            (1, 0): 3.1,
            (1, 1): 3.2,
        },
        "X2_b_sin": {
            (0, 1): 4.0,
            (1, 0): 4.1,
            (1, 1): 4.2,
        },
    }
    flipped_parameters = func(parameters)
    theta = np.linspace(0, 2 * np.pi, 4, endpoint=False)
    zeta = np.linspace(0, 2 * np.pi, 3, endpoint=False)
    t, z = np.meshgrid(theta, zeta, indexing="ij")
    tref, zref = transform(t, z)
    for i in [1, 2]:
        result, reference = np.zeros((2, *t.shape))
        for key, func in zip([f"X{i}_b_sin", f"X{i}_b_cos"], [np.sin, np.cos]):
            for m, n in flipped_parameters[key]:
                result += flipped_parameters[key][m, n] * func(m * t - n * z)
            for m, n in parameters[key]:
                reference += parameters[key][m, n] * func(m * tref - n * zref)
        np.testing.assert_allclose(result, reference, err_msg=f"with variable X{i}")


def test_boundary_direction():
    parameters = {
        "X1_b_cos": {
            (0, 0): 1.0,
            (0, 1): 1.1,
            (1, 0): 1.2,
            (1, 1): 1.3,
        },
        "X2_b_cos": {
            (0, 0): 2.0,
            (0, 1): 2.1,
            (1, 0): 2.2,
            (1, 1): 2.3,
        },
        "X1_b_sin": {
            (0, 1): 3.0,
            (1, 0): 3.1,
            (1, 1): 3.2,
        },
        "X2_b_sin": {
            (0, 1): 4.0,
            (1, 0): 4.1,
            (1, 1): 4.2,
        },
    }
    flipped_parameters = util.flip_boundary_theta(parameters)
    A1 = util.signed_cross_sectional_area(parameters, 0.0)
    s1 = util.check_boundary_direction(parameters)
    A2 = util.signed_cross_sectional_area(flipped_parameters, 0.0)
    s2 = util.check_boundary_direction(flipped_parameters)
    assert s2 and not s1
    assert A1 < 0 and A2 > 0


@pytest.mark.parametrize("npoints", [(21, 31), (31, 30), (201, 101)])
@pytest.mark.parametrize("endpoint", [False, True])
def test_linking_number(npoints, endpoint):
    npoints_a = npoints[0]
    curve_a = np.zeros((npoints_a, 3))
    theta = np.linspace(0, 2 * np.pi, npoints_a, endpoint=endpoint)
    curve_a[:, 0] = 1.01 * np.cos(theta)
    curve_a[:, 1] = 0.99 * np.sin(theta)
    curve_a[:, 2] = 0.01
    npoints_b = npoints[1]
    curve_b = np.zeros((npoints_b, 3))
    theta = np.linspace(0, 2 * np.pi, npoints_b, endpoint=endpoint)
    curve_b[:, 1] = -0.01
    curve_b[:, 2] = 0.49 * np.sin(theta)

    curve_b[:, 0] = 1.4 - 0.51 * np.cos(theta)
    assert np.isclose(util.linking_number(curve_a, curve_b, endpoint=endpoint), 1.0), (
        "linking number of two linked circles should be 1.0"
    )

    curve_b[:, 0] = 1.4 + 0.51 * np.cos(theta)
    assert np.isclose(util.linking_number(curve_a, curve_b, endpoint=endpoint), -1.0), (
        "linking number of two linked circles with opposite orientation should be -1.0"
    )

    curve_b[:, 0] = 1.6 - 0.51 * np.cos(theta)
    assert np.isclose(util.linking_number(curve_a, curve_b, endpoint=endpoint), 0.0), (
        "linking number of two unlinked circles should be 0.0"
    )


@pytest.mark.parametrize(
    "case, Lk_expected",
    [
        ("ellip_cyl_helix", 0),
        ("ellip_cyl_helix_rot", 1),
        ("ellip_cyl_helix_rot2", 0),
    ],
    ids=["no_link", "link", "no_link2"],
)
@pytest.mark.parametrize("endpoint", [False, True])
def test_linking_number_boundary(case, Lk_expected, endpoint):
    params = gvec.util.boundary_generator(case)
    nfp = 3
    params["nfp"] = nfp
    theta = np.linspace(0, 2 * np.pi, 21, endpoint=False)
    for phi_dir in [-1, 1]:
        zeta = np.linspace(0, 2 * np.pi, params["nfp"] * 81, endpoint=endpoint)
        X1, X2 = gvec.util.evaluate_boundary(theta, zeta, params)
        xyz_surf = np.zeros((len(zeta), len(theta), 3))
        xyz_surf[:, :, 0] = X1.T * np.cos(phi_dir * zeta[:, None])
        xyz_surf[:, :, 1] = X1.T * np.sin(phi_dir * zeta[:, None])
        xyz_surf[:, :, 2] = X2.T
        Lk = gvec.util.linking_number(
            xyz_surf[:, 0, :], xyz_surf[:, len(theta) // 2, :], endpoint=endpoint
        )
        assert np.abs(Lk - phi_dir * Lk_expected * nfp) < 1e-8, (
            f"Linking number of the surface is not as expected! Got {Lk}, expected {phi_dir * Lk_expected * nfp}, for phi direction {-phi_dir}"
        )
