import pytest

try:
    import numpy as np
    from typing import Literal
    from gvec.coils import Coil, CoilSet, trace_fieldlines
    import tempfile
    from scipy.constants import mu_0
    from scipy.special import ellipk, ellipe

    rng = np.random.default_rng(0)
except ImportError:
    pass


@pytest.fixture
def coil_parameters():
    params = dict(n_coilpoints=50000, n_evals=1000, coil_radius=1, coil_current=1, n_coils=10)
    return params


@pytest.fixture
def off_axis_eval_pos(coil_parameters):
    pos = rng.uniform(-0.8, 0.8, [3, coil_parameters["n_evals"]])
    return pos


@pytest.fixture
def on_axis_eval_pos(coil_parameters):
    z = np.linspace(-1, 1, coil_parameters["n_evals"])
    pos = np.zeros((3, len(z)))
    pos[2, :] = z
    return (pos, z)


@pytest.fixture
def circular_coil(coil_parameters):
    params = coil_parameters
    circ_coil = create_circular_coil(
        params["coil_radius"], params["coil_current"], params["n_coilpoints"]
    )
    return circ_coil


@pytest.fixture
def circular_coil_set(coil_parameters):
    params = coil_parameters
    circ_coils = []
    shifts = np.arange(params["n_coils"])
    for shift in shifts:
        circ_coils.append(
            create_circular_coil(
                params["coil_radius"],
                params["coil_current"],
                params["n_coilpoints"],
                shift=shift,
            )
        )
    return CoilSet(circ_coils), shifts


def circular_coil_axis_mod_B(z, R, coil_current):
    mod_B = mu_0 * coil_current * R**2 / (2 * (z**2 + R**2) ** (3 / 2))
    return mod_B


def circular_coil_A(r, theta, phi, R, coil_current):
    """Vector potential of a circular coil in the xy plane

    Parameters
    ----------
    r : ArrayLike
        Spherical radius.
    theta : ArrayLike
        Polar angle
    phi : ArrayLike
        Azimuthal angle
    R : float
        Coil radius
    coil_current : float
        Coil current

    Returns
    -------
    np.ndarray
        A in Cartesian coordiantes.
    """
    k_2 = 4 * R * r * np.sin(theta) / (R**2 + r**2 + 2 * R * r * np.sin(theta))

    prefactor = mu_0 * coil_current * R / np.pi

    A_phi = (
        prefactor
        / np.sqrt(R**2 + r**2 + 2 * R * r * np.sin(theta))
        * (((2 - k_2) * ellipk(k_2) - 2 * ellipe(k_2)) / k_2)
    )
    e_phi = np.zeros([3, len(r)])
    e_phi[0, :] = -np.sin(phi)
    e_phi[1, :] = np.cos(phi)
    A = A_phi * e_phi
    return A


def circular_coil_B(r, theta, phi, R, coil_current):
    """Magnetic field of a circular coil in the xy plane

    Parameters
    ----------
    r : ArrayLike
        Spherical radius.
    theta : ArrayLike
        Polar angle
    phi : ArrayLike
        Azimuthal angle
    R : float
        Coil radius
    coil_current : float
        Coil current

    Returns
    -------
    np.ndarray
        B in Cartesian coordiantes.
    """

    prefactor = mu_0 * coil_current / np.pi
    alpha_2 = R**2 + r**2 - 2 * R * r * np.sin(theta)
    beta_2 = R**2 + r**2 + 2 * R * r * np.sin(theta)
    beta = np.sqrt(beta_2)

    k_2 = 1 - alpha_2 / beta_2

    B_r = prefactor * R**2 * np.cos(theta) / (alpha_2 * beta) * ellipe(k_2)
    B_theta = (
        prefactor
        / (2 * alpha_2 * beta * np.sin(theta))
        * ((r**2 + R**2 * np.cos(2 * theta)) * ellipe(k_2) - alpha_2 * ellipk(k_2))
    )

    e_r = np.zeros([3, len(r)])
    e_r[0, :] = np.sin(theta) * np.cos(phi)
    e_r[1, :] = np.sin(theta) * np.sin(phi)
    e_r[2, :] = np.cos(theta)

    e_theta = np.zeros([3, len(r)])
    e_theta[0, :] = np.cos(theta) * np.cos(phi)
    e_theta[1, :] = np.cos(theta) * np.sin(phi)
    e_theta[2, :] = -np.sin(theta)

    B = B_r * e_r + B_theta * e_theta

    return B


def create_circular_coil(radius, current, n_coilpoints, shift=0):
    phi = np.linspace(0, 2 * np.pi, n_coilpoints, endpoint=True)
    coil_points = np.zeros((3, n_coilpoints))
    coil_points[0, :] = radius * np.cos(phi)
    coil_points[1, :] = radius * np.sin(phi)
    coil_points[2, :] = shift
    circ_coil = Coil(coil_points, current)
    return circ_coil


def test_BiotSavart_B_axis_circular_coil(coil_parameters, on_axis_eval_pos, circular_coil):
    params = coil_parameters
    pos, z = on_axis_eval_pos

    mod_B_analytic = circular_coil_axis_mod_B(z, params["coil_radius"], params["coil_current"])
    ds = circular_coil.eval_mod_B(pos)

    assert np.all(abs(ds.B.sel(xyz=["x", "y"])) <= 1e-12)
    np.testing.assert_allclose(ds.mod_B, mod_B_analytic)


def test_BiotSavart_B_axis_circular_coil_set(
    coil_parameters, on_axis_eval_pos, circular_coil_set
):
    params = coil_parameters
    pos, z = on_axis_eval_pos
    circ_coil_set, shifts = circular_coil_set
    mod_B_analytic = 0.0
    for shift in shifts:
        mod_B_analytic += circular_coil_axis_mod_B(
            z - shift, params["coil_radius"], params["coil_current"]
        )
    pos = np.zeros((3, len(z)))
    pos[2, :] = z

    ds = circ_coil_set.eval_mod_B(pos)

    assert np.all(abs(ds.B.sel(xyz=["x", "y"])) <= 1e-12)
    np.testing.assert_allclose(ds.mod_B, mod_B_analytic)


def test_BiotSavart_A_circular_coil(coil_parameters, off_axis_eval_pos, circular_coil):
    params = coil_parameters
    pos = off_axis_eval_pos
    r = np.sqrt(np.sum(pos**2, axis=0))
    theta = np.arccos(pos[2, :] / r)
    phi = np.sign(pos[1, :] * np.arccos(pos[0, :] / np.sqrt(pos[0, :] ** 2, pos[1, :] ** 2)))

    A = circular_coil_A(
        r=r, theta=theta, phi=phi, R=params["coil_radius"], coil_current=params["coil_current"]
    )
    mod_A_analytic = np.sqrt(np.sum(A**2, axis=0))
    ds = circular_coil.eval_mod_A(pos)

    assert np.all(abs(ds.A.sel(xyz=["z"])) <= 1e-12)
    np.testing.assert_allclose(ds.mod_A, mod_A_analytic)


def test_BiotSavart_B_circular_coil(coil_parameters, off_axis_eval_pos, circular_coil):
    params = coil_parameters
    pos = off_axis_eval_pos

    r = np.sqrt(np.sum(pos**2, axis=0))
    theta = np.arccos(pos[2, :] / r)
    phi = np.sign(pos[1, :] * np.arccos(pos[0, :] / np.sqrt(pos[0, :] ** 2, pos[1, :] ** 2)))

    B = circular_coil_B(
        r=r, theta=theta, phi=phi, R=params["coil_radius"], coil_current=params["coil_current"]
    )
    mod_B_analytic = np.sqrt(np.sum(B**2, axis=0))
    ds = circular_coil.eval_mod_B(pos)

    np.testing.assert_allclose(ds.mod_B, mod_B_analytic)


def test_BiotSavart_A_circular_coil_set(coil_parameters, off_axis_eval_pos, circular_coil_set):
    params = coil_parameters
    circ_coil_set, shifts = circular_coil_set
    A_analytic = 0.0
    pos = off_axis_eval_pos
    pos_aux = pos.copy()
    for shift in shifts:
        pos_aux[2, :] = pos[2, :] - shift
        r = np.sqrt(np.sum(pos_aux**2, axis=0))
        theta = np.arccos(pos_aux[2, :] / r)
        phi = np.sign(
            pos_aux[1, :]
            * np.arccos(pos_aux[0, :] / np.sqrt(pos_aux[0, :] ** 2, pos_aux[1, :] ** 2))
        )
        A_analytic += circular_coil_A(
            r=r,
            theta=theta,
            phi=phi,
            R=params["coil_radius"],
            coil_current=params["coil_current"],
        )

    mod_A_analytic = np.sqrt(np.sum(A_analytic**2, axis=0))
    ds = circ_coil_set.eval_mod_A(pos)

    assert np.all(abs(ds.A.sel(xyz=["z"])) <= 1e-12)
    np.testing.assert_allclose(ds.mod_A, mod_A_analytic)


def test_BiotSavart_B_circular_coil_set(coil_parameters, off_axis_eval_pos, circular_coil_set):
    params = coil_parameters
    pos = off_axis_eval_pos
    circ_coil_set, shifts = circular_coil_set
    B_analytic = 0.0
    pos_aux = pos.copy()
    for shift in shifts:
        pos_aux[2, :] = pos[2, :] - shift
        r = np.sqrt(np.sum(pos_aux**2, axis=0))
        theta = np.arccos(pos_aux[2, :] / r)
        phi = np.sign(
            pos_aux[1, :]
            * np.arccos(pos_aux[0, :] / np.sqrt(pos_aux[0, :] ** 2, pos_aux[1, :] ** 2))
        )
        B_analytic += circular_coil_B(
            r=r,
            theta=theta,
            phi=phi,
            R=params["coil_radius"],
            coil_current=params["coil_current"],
        )

    mod_B_analytic = np.sqrt(np.sum(B_analytic**2, axis=0))
    ds = circ_coil_set.eval_mod_B(pos)

    np.testing.assert_allclose(ds.mod_B, mod_B_analytic)


def test_coil_save_and_load(circular_coil):
    with tempfile.TemporaryDirectory() as path_tmp:
        circular_coil.save(path_tmp + "/coil_test.nc")
        circ_coil_2 = Coil.load(path_tmp + "/coil_test.nc")
        np.testing.assert_array_equal(circ_coil_2.coil_points, circular_coil.coil_points)
        assert circular_coil.coil_current == circ_coil_2.coil_current


def test_coilset_save_and_load(circular_coil_set):
    circ_coil_set = circular_coil_set[0]
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


def test_fieldlines(circular_coil_set):
    coil_set = circular_coil_set[0]
    starts = np.zeros([3, 2])
    starts[0, :] = [0.0, 0.5]
    surface_normals = [np.array([0, 0, 1])]
    surface_points = [np.array([0, 0, 0.5])]

    dt = trace_fieldlines(
        coils=coil_set,
        starts=starts,
        time=10,
        surf_normals=surface_normals,
        surf_points=surface_points,
    )

    for fieldline in dt:
        ds = dt[fieldline]
        assert "pos" in ds
        assert "t" in ds
        assert "xyz" in ds
        for i in range(len(surface_normals)):
            assert f"event_{i}" in ds
