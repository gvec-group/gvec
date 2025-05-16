import pytest

try:
    from gvec import util
    from gvec.vtk import ev2vtk
    import gvec
    import numpy as np
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
    assert cid == cid.copy()
    assert cid == util.CaseInsensitiveDict(cid.lower_items())
    with pytest.raises(KeyError):
        cid["3"]
    with pytest.raises(ValueError):
        cid.update({"d": 1, "D": 2})
    with pytest.raises(ValueError):
        _ = util.CaseInsensitiveDict({"a": 1, "A": 2})
    with pytest.raises(ValueError):
        cid == {"a": 1, "A": 2}
    assert "a" in cid
    assert "A" in cid
    assert "b" in cid


def test_ev2vtk(testcaserundir, testfiles):
    rho = np.linspace(0, 1, 5)
    rho[0] = 1e-4  # avoid evaluation at rho=0
    theta = np.linspace(0, 2 * np.pi, 5)  # including endpoints
    zeta = np.linspace(0, 2 * np.pi, 5)  # full torus, including endpoints
    vars_out = ["X1", "X2", "LA", "iota", "p", "pos", "grad_zeta", "e_rho", "B"]

    with gvec.State(*testfiles) as state:
        ev = gvec.Evaluations(rho=rho, theta=theta, zeta=zeta, state=state)
        state.compute(ev, *vars_out)

    ev = ev[vars_out]
    ev2vtk(testcaserundir / "test_ev2.vtk", ev)


@pytest.mark.parametrize(
    "func, transform",
    [
        (util.flip_parameters_theta, lambda t, z: (-t, z)),
        (util.flip_parameters_zeta, lambda t, z: (t, -z)),
    ],
    ids=["theta", "zeta"],
)
def test_flip_parameters(func, transform):
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
    print(flipped_parameters)
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
