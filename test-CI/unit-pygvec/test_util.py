import pytest

try:
    from gvec import util
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

def test_xr2vtk(tmpdir, testcaserundir,testfiles):
    with gvec.State(*testfiles) as state_new:
        ev = gvec.Evaluations(rho=4, theta=5, zeta=5, state=state_new)
        state_new.compute(ev, "X1", "X2", "LA", "iota", "p","pos","grad_zeta","e_rho","B")
    ev_out=ev[["X1", "X2", "LA", "iota", "p","pos","grad_zeta","e_rho","B"]]
    
    # check if any of the variables contain NaNs
    list_vars = []
    for var in ev_out:
        if ev_out[var].isnull().any():
            continue
        else:
            list_vars.append(var)
    ev_out = ev_out[list_vars]
    text = gvec.util.xr2vtk(str(testcaserundir)+"/test_xr2vtk", ev_out)