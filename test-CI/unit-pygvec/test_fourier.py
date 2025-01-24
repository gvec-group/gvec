import pytest

try:
    import numpy as np
    from numpy.random import random

    from gvec import fourier, State, Evaluations, compute
except ImportError:
    pytest.skip("Import Error", allow_module_level=True)


# === Fixtures === #


@pytest.fixture(
    params=[(8, 10), (9, 10), (8, 11), (9, 11)], ids=["8x10", "9x10", "8x11", "9x11"]
)
def shape2d(request):
    return request.param


@pytest.fixture(params=[(5, 3), (3, 5), (5, 5)], ids=["5x3", "3x5", "5x5"])
def MN(request):
    return request.param


@pytest.fixture(
    params=[(12, 14), (13, 14), (12, 15), (13, 15)],
    ids=["12x14", "13x14", "12x15", "13x15"],
)
def points2d(request):
    return request.param


@pytest.fixture()
def state(testfiles):
    with State(*testfiles) as state:
        yield state


@pytest.fixture()
def ev(state):
    rho = [0.1, 0.5, 0.9]
    ds = Evaluations(rho=rho, theta=20, zeta=50, state=state)
    compute(ds, "mod_B", state=state)
    return ds


# === Test FFT === #


@pytest.mark.parametrize(
    "c, s",
    [
        ([0.0, 1.0, 0.1], []),
        ([0.0, 1.0, -0.1], [0.0, 0.5, -0.2, 0.1]),
        ([2.0, 1.0, -0.1], [0.0, 1.0, -0.2, 0.1]),
        (random(5), random(5)),
        (random(5), random(5)),
    ],
    ids=["cos", "sincos", "offset", "random-1", "random-2"],
)
@pytest.mark.parametrize("points", [10, 11], ids=["even", "odd"])
def test_fft1d(c, s, points: tuple[int, int]):
    t = np.linspace(0, 2 * np.pi, points, endpoint=False)
    x = sum([ci * np.cos(i * t) for i, ci in enumerate(c)])
    x += sum([si * np.sin(i * t) for i, si in enumerate(s)])

    xc, xs = fourier.fft1d(x)
    assert np.allclose(xc[: len(c)], c)
    assert np.allclose(xs[1 : len(s)], s[1:])
    assert xs[0] == 0

    y = sum([ci * np.cos(i * t) for i, ci in enumerate(xc)])
    y += sum([si * np.sin(i * t) for i, si in enumerate(xs)])
    assert np.allclose(x, y)


def test_fft2d_modes(MN):
    m, n = fourier.fft2d_modes(*MN)
    assert m.shape == (MN[0] + 1,)
    assert n.shape == (2 * MN[1] + 1,)
    assert m.min() == 0 and m.max() == MN[0]
    assert abs(n.min()) - n.max() in [0, 1]


def test_fft2d_modes_grid(MN):
    M, N = MN
    m, n = fourier.fft2d_modes(M, N, grid=True)
    assert m.shape == n.shape == (M + 1, 2 * N + 1)
    assert m.min() == 0 and m.max() == M
    assert abs(n.min()) - n.max() in [0, 1]
    assert abs(n.max()) == N


@pytest.mark.parametrize(
    "dM2, dN2",
    [(2, 0), (0, 2), (2, 2), (-2, 0), (0, -2), (-2, -2)],
    ids=["M+", "N+", "MN+", "M-", "N-", "MN-"],
)
def test_scale_modes2d(MN, dM2: int, dN2: int):
    M1, N1 = MN
    M2, N2 = M1 + dM2, N1 + dN2

    m1, n1 = fourier.fft2d_modes(M1, N1, grid=True)
    m2, n2 = fourier.fft2d_modes(M2, N2, grid=True)
    Mmin, Nmin = min(M1, M2), min(N1, N2)

    c1 = random((M1 + 1, 2 * N1 + 1))
    c2 = fourier.scale_modes2d(c1, M2, N2)
    assert c2.shape == (M2 + 1, 2 * N2 + 1)
    assert np.all(
        c1[(m1 <= Mmin) & (np.abs(n1) <= Nmin)]
        == c2[(m2 <= Mmin) & (np.abs(n2) <= Nmin)]
    )
    assert np.all(c2[(m2 > M1) & (np.abs(n2) > N1)] == 0)


def test_fft2d(MN, points2d):
    t = np.linspace(0, 2 * np.pi, points2d[0], endpoint=False)
    z = np.linspace(0, 2 * np.pi, points2d[1], endpoint=False)
    T, Z = np.meshgrid(t, z, indexing="ij")
    M, N = MN
    ms, ns = fourier.fft2d_modes(M, N)
    c = random((M + 1, 2 * N + 1))
    c[0, ns < 0] = 0
    s = random((M + 1, 2 * N + 1))
    s[0, ns <= 0] = 0
    x = sum(
        [
            sum(
                [
                    c[m, n] * np.cos(m * T - n * Z) + s[m, n] * np.sin(m * T - n * Z)
                    for m in ms
                ]
            )
            for n in ns
        ]
    )

    xc, xs = fourier.fft2d(x)
    xM, xN = xc.shape[0] - 1, xc.shape[1] // 2
    assert xc.shape == xs.shape == (xM + 1, 2 * xN + 1)

    c = fourier.scale_modes2d(c, xM, xN)
    s = fourier.scale_modes2d(s, xM, xN)
    assert np.allclose(c, xc)
    assert np.allclose(s, xs)


def test_ifft2d(MN):
    M, N = MN
    ms, ns = fourier.fft2d_modes(M, N)
    c = random((M + 1, 2 * N + 1))
    c[0, ns < 0] = 0
    s = random((M + 1, 2 * N + 1))
    s[0, ns <= 0] = 0

    x = fourier.ifft2d(c, s)
    assert x.shape == (2 * N + 1, 2 * M + 1)

    t = np.linspace(0, 2 * np.pi, 2 * M + 1, endpoint=False)
    z = np.linspace(0, 2 * np.pi, 2 * N + 1, endpoint=False)
    T, Z = np.meshgrid(t, z, indexing="ij")
    ref = sum(
        [
            sum(
                [
                    c[m, n] * np.cos(m * T - n * Z) + s[m, n] * np.sin(m * T - n * Z)
                    for m in ms
                ]
            )
            for n in ns
        ]
    )
    assert np.allclose(x, ref.T)


def test_ifft2d_fft2d(MN):
    M, N = MN
    ms, ns = fourier.fft2d_modes(M, N)
    c = random((M + 1, 2 * N + 1))
    c[0, ns < 0] = 0
    s = random((M + 1, 2 * N + 1))
    s[0, ns <= 0] = 0

    x = fourier.ifft2d(c, s)
    xc, xs = fourier.fft2d(x)
    assert np.allclose(c, xc)
    assert np.allclose(s, xs)


@pytest.mark.parametrize(
    "loop",
    [True, False],
    ids=["loop", "mesh"],
)
def test_eval2d(MN, points2d, loop):
    t = np.linspace(0, 2 * np.pi, points2d[0], endpoint=False)
    z = np.linspace(0, 2 * np.pi, points2d[1], endpoint=False)
    T, Z = np.meshgrid(t, z, indexing="ij")
    M, N = MN
    ms, ns = fourier.fft2d_modes(M, N)
    c = random((M + 1, 2 * N + 1))
    c[0, ns < 0] = 0
    s = random((M + 1, 2 * N + 1))
    s[0, ns <= 0] = 0

    x = sum(
        [
            sum(
                [
                    c[m, n] * np.cos(m * T - n * Z) + s[m, n] * np.sin(m * T - n * Z)
                    for m in ms
                ]
            )
            for n in ns
        ]
    )
    xe = fourier.eval2d(c, s, T, Z, loop=loop)
    assert np.allclose(x, xe)


def test_ev2ft_2d(ev):
    ft = fourier.ev2ft(ev[["mod_B"]].isel(rad=0))
    assert set(ft.dims) == {"m", "n"}
    assert set(ft.data_vars) == {"mod_B_mnc", "mod_B_mns"}


def test_ev2ft_3d(ev):
    ft = fourier.ev2ft(ev[["mod_B"]])
    assert set(ft.dims) == {"rad", "m", "n"}
    assert set(ft.data_vars) == {"mod_B_mnc", "mod_B_mns"}
