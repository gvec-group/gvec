import pytest

try:
    import numpy as np
    import xarray as xr
    from gvec.coils import Coil, CoilSet
    import tempfile
except ImportError:
    pytest.skip("Import Error", allow_module_level=True)

mu_0 = 4 * np.pi * 0.99999999987 * 1e-7


def circular_coil_axis_mod_B(z, R, coil_current):
    mod_B = mu_0 * coil_current * R**2 / (2 * (z**2 + R**2) ** (3 / 2))
    return mod_B


def create_circular_coil(radius, current, n_coilpoints, shift=0):
    phi = np.linspace(0, 2 * np.pi, n_coilpoints, endpoint=True)
    coil_points = np.zeros((3, n_coilpoints))
    coil_points[0, :] = radius * np.cos(phi)
    coil_points[1, :] = radius * np.sin(phi)
    coil_points[2, :] = shift
    circ_coil = Coil(coil_points, current)
    return circ_coil


def test_BiotSavart_circular_coil():
    n_coilpoints = 10000
    n_evals = 1000
    coil_radius = 1
    coil_current = 1
    circ_coil = create_circular_coil(coil_radius, coil_current, n_coilpoints)
    z = np.linspace(-1, 1, n_evals)

    pos = np.zeros((3, len(z)))
    pos[2, :] = z

    mod_B_analytic = circular_coil_axis_mod_B(z, coil_radius, coil_current)
    ds = circ_coil.eval_mod_B(pos)

    assert np.all(abs(ds.B.sel(xyz=["x", "y"])) <= 1e-12)
    np.testing.assert_allclose(ds.mod_B, mod_B_analytic)


def test_BiotSavart_circular_coil_set():
    n_coilpoints = 10000
    n_evals = 1000
    coil_radius = 1
    coil_current = 1
    n_coils = 10
    z = np.linspace(-1, 1, n_evals)
    circ_coils = []
    mod_B_analytic = 0.0
    for shift in range(n_coils):
        mod_B_analytic += circular_coil_axis_mod_B(z - shift, coil_radius, coil_current)
        circ_coils.append(
            create_circular_coil(coil_radius, coil_current, n_coilpoints, shift=shift)
        )

    circ_coil_set = CoilSet(circ_coils)

    pos = np.zeros((3, len(z)))
    pos[2, :] = z

    ds = circ_coil_set.eval_mod_B(pos)

    assert np.all(abs(ds.B.sel(xyz=["x", "y"])) <= 1e-12)
    np.testing.assert_allclose(ds.mod_B, mod_B_analytic)


def test_coil_save_and_load():
    n_coilpoints = 100
    coil_radius = 1
    coil_current = 1
    circ_coil = create_circular_coil(coil_radius, coil_current, n_coilpoints)
    with tempfile.TemporaryDirectory() as path_tmp:
        circ_coil.save(path_tmp + "/coil_test.nc")
        circ_coil_2 = Coil.load(path_tmp + "/coil_test.nc")
        np.testing.assert_array_equal(circ_coil_2.coil_points, circ_coil.coil_points)
        assert circ_coil.coil_current == circ_coil_2.coil_current


def test_coilset_save_and_load():
    n_coilpoints = 100
    coil_radius = 1
    coil_current = 1
    n_coils = 10
    circ_coils = []
    for shift in range(n_coils):
        circ_coils.append(
            create_circular_coil(coil_radius, coil_current, n_coilpoints, shift=shift)
        )
    circ_coil_set = CoilSet(circ_coils)
    with tempfile.TemporaryDirectory() as path_tmp:
        try:
            circ_coil_set.save(path_tmp + "/coil_test.nc")
        except ModuleNotFoundError:
            circ_coil_set.save(path_tmp + "/coil_test.nc", engine="netcdf4")
        circ_coil_set_2 = CoilSet.load(path_tmp + "/coil_test.nc")
        for coil in circ_coil_set.coils:
            np.testing.assert_array_equal(
                circ_coil_set[coil].coil_points, circ_coil_set_2[coil].coil_points
            )
            assert circ_coil_set[coil].coil_current == circ_coil_set_2[coil].coil_current
