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
