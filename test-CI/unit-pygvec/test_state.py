import pytest

try:
    import numpy as np

    import gvec
    from gvec.state import State
except ImportError:
    pytest.skip("Import Error", allow_module_level=True)


# === Fixtures === #


@pytest.fixture()
def teststate(testfiles):
    with State(*testfiles) as state:
        yield state


# === Tests === #


def test_version():
    import pkg_resources

    assert isinstance(gvec.__version__, str)
    assert gvec.__version_tuple__ >= (0, 2, 1)
    assert pkg_resources.get_distribution("gvec").version == gvec.__version__


def test_state(testfiles):
    with State(*testfiles) as state:
        assert isinstance(state, State)


def test_state_args(testfiles):
    paramfile, statefile = testfiles

    with pytest.raises(FileNotFoundError):
        state = State(paramfile, "nonexistent.dat")

    with pytest.raises(FileNotFoundError):
        state = State("nonexistent.ini", statefile)


def test_state_explicit(testfiles):
    state = State(*testfiles)
    assert isinstance(state, State)
    assert state.initialized
    state.finalize()
    assert not state.initialized

    with pytest.raises(RuntimeError):
        state.evaluate_base_tens("X1", None, [0.5], [0.5], [0.5])


def test_state_twice(testfiles):
    # double context
    with State(*testfiles) as state:
        pass
    with State(*testfiles) as state:
        pass

    # double explicit
    state = State(*testfiles)
    state.finalize()
    state = State(*testfiles)
    state.finalize()

    # simultaneous context (not implemented yet)
    with pytest.raises(NotImplementedError):
        with State(*testfiles) as state1:
            with State(*testfiles) as state2:
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


def test_evaluate_profile(teststate):
    rho = np.linspace(0, 1, 6)

    iota = teststate.evaluate_profile("iota", rho)
    assert iota.shape == rho.shape
    dp_dr = teststate.evaluate_profile("p_prime", rho)
