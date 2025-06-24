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
    rho = np.linspace(0, 1, 5)
    rho[0] = 1e-4  # avoid evaluation at rho=0
    theta = np.linspace(0, 2 * np.pi, 5)  # including endpoints
    zeta = np.linspace(0, 2 * np.pi, 5)  # full torus, including endpoints
    vars_out = ["X1", "X2", "LA", "iota", "p", "pos", "grad_zeta", "e_rho", "B"]

    state = gvec.State(*testfiles)
    ev = gvec.Evaluations(rho=rho, theta=theta, zeta=zeta, state=state)
    state.compute(ev, *vars_out)

    ev = ev[vars_out]
    ev2vtk(testcaserundir / "test_ev2.vtk", ev)


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
