# Copyright (c) 2025 GVEC Contributors, Max Planck Institute for Plasma Physics
# License: MIT
from collections.abc import Mapping, MutableMapping
from copy import deepcopy
from typing import Literal

import numpy as np

from .params import CaseInsensitiveDict


def signed_cross_sectional_area(
    parameters: Mapping, zeta: float, resolution: int = 1000
) -> float:
    t = np.linspace(0, 2 * np.pi, resolution, endpoint=False)
    x1 = np.zeros_like(t)
    dx1dt = np.zeros_like(t)
    x2 = np.zeros_like(t)
    dx2dt = np.zeros_like(t)
    nfp = parameters.get("nfp", 1)
    for (m, n), value in parameters.get("X1_b_cos", {}).items():
        x1 += value * np.cos(m * t - n * nfp * zeta)
        dx1dt -= value * m * np.sin(m * t - n * nfp * zeta)
    for (m, n), value in parameters.get("X1_b_sin", {}).items():
        x1 += value * np.sin(m * t - n * nfp * zeta)
        dx1dt += value * m * np.cos(m * t - n * nfp * zeta)
    for (m, n), value in parameters.get("X2_b_cos", {}).items():
        x2 += value * np.cos(m * t - n * nfp * zeta)
        dx2dt -= value * m * np.sin(m * t - n * nfp * zeta)
    for (m, n), value in parameters.get("X2_b_sin", {}).items():
        x2 += value * np.sin(m * t - n * nfp * zeta)
        dx2dt += value * m * np.cos(m * t - n * nfp * zeta)
    dA = x1 * dx2dt - x2 * dx1dt
    return np.sum(dA)


def check_boundary_direction(parameters: Mapping) -> bool:
    """Determine whether the boundary is described by right-handed logical coordinates (θ,ζ).

    GVEC requires a right-handed logical coordinate system (ρ,θ,ζ).
    The logical coordinate system of the poloidal plane, (ρ,θ) is also required to be right-handed,
    which requires the poloidal angle to increase in the counter-clockwise direction.
    As a consequence the toroidal angle has to increase in the clockwise direction when viewed from above.
    This is ensured in the definition of the h-maps.

    Returns:
        bool: True if (ρ,θ) is right-handed / θ increases counter-clockwise, False otherwise.
    """
    return signed_cross_sectional_area(parameters, 0.0) > 0


def effective_minor_radius(
    parameters: Mapping,
    resolution: tuple[int, int] = (1000, 100),
):
    nfp = parameters.get("nfp", 1)
    areas = np.zeros(resolution[1])
    for z, zeta in enumerate(np.linspace(0, 2 * np.pi / nfp, resolution[1], endpoint=False)):
        areas[z] = abs(signed_cross_sectional_area(parameters, zeta, resolution=resolution[0]))
    return np.sqrt(np.mean(areas) / np.pi)


def evaluate_boundary(
    theta: np.ndarray, zeta: np.ndarray, parameters: Mapping
) -> tuple[np.ndarray, np.ndarray]:
    """Evaluate the boundary at the given (theta, zeta) points.

    Args:
        theta (1D np.ndarray): The poloidal angles at which to evaluate the boundary.
        zeta (1D np.ndarray): The toroidal angles at which to evaluate the boundary.
        parameters (Mapping): The parameters defining the boundary.

    Returns:
        tuple[2D np.ndarray, 2D np.ndarray]: The (X^1, X^2) coordinates of the boundary at the given (theta, zeta) points.
    """
    theta = np.asarray(theta)
    zeta = np.asarray(zeta)
    if theta.ndim != 1 or zeta.ndim != 1:
        raise ValueError("theta and zeta must be 1D arrays")
    nfp = parameters.get("nfp", 1)
    theta, zeta = np.meshgrid(theta, zeta, indexing="ij")
    x1 = np.zeros_like(theta)
    x2 = np.zeros_like(theta)
    for (m, n), value in parameters.get("X1_b_cos", {}).items():
        x1 += value * np.cos(m * theta - n * nfp * zeta)
    for (m, n), value in parameters.get("X1_b_sin", {}).items():
        x1 += value * np.sin(m * theta - n * nfp * zeta)
    for (m, n), value in parameters.get("X2_b_cos", {}).items():
        x2 += value * np.cos(m * theta - n * nfp * zeta)
    for (m, n), value in parameters.get("X2_b_sin", {}).items():
        x2 += value * np.sin(m * theta - n * nfp * zeta)
    return x1, x2


def evaluate_axis(zeta: np.ndarray, parameters: Mapping) -> tuple[np.ndarray, np.ndarray]:
    """Evaluate the magnetic axis at the given zeta points.

    Args:
        zeta (1D np.ndarray): The toroidal angles at which to evaluate the axis.
        parameters (Mapping): The parameters defining the axis.

    Returns:
        tuple[1D np.ndarray, 1D np.ndarray]: The (X^1, X^2) coordinates of the axis at the given zeta points.
    """
    zeta = np.asarray(zeta)
    if zeta.ndim != 1:
        raise ValueError("zeta must be a 1D array")
    nfp = parameters.get("nfp", 1)
    x1 = np.zeros_like(zeta)
    x2 = np.zeros_like(zeta)
    for (m, n), value in parameters.get("X1_a_cos", {}).items():
        if m != 0:
            raise ValueError("Axis X1_a_cos should only have m=0 modes")
        x1 += value * np.cos(-n * nfp * zeta)
    for (m, n), value in parameters.get("X1_a_sin", {}).items():
        if m != 0:
            raise ValueError("Axis X1_a_sin should only have m=0 modes")
        x1 += value * np.sin(-n * nfp * zeta)
    for (m, n), value in parameters.get("X2_a_cos", {}).items():
        if m != 0:
            raise ValueError("Axis X2_a_cos should only have m=0 modes")
        x2 += value * np.cos(-n * nfp * zeta)
    for (m, n), value in parameters.get("X2_a_sin", {}).items():
        if m != 0:
            raise ValueError("Axis X2_a_sin should only have m=0 modes")
        x2 += value * np.sin(-n * nfp * zeta)
    return x1, x2


def compute_boundary_perturbation(
    base_parameters: Mapping, perturbed_parameters: Mapping
) -> tuple[CaseInsensitiveDict, CaseInsensitiveDict]:
    """Computes the difference between the perturbed and base boundary parameters as a perturbation."""
    new_base = CaseInsensitiveDict()
    new_perturbed = CaseInsensitiveDict()
    for i in [1, 2]:
        for sincos in ["sin", "cos"]:
            perturbed = {}
            base = {}
            # set boundary modes to values from restart
            for (m, n), v in base_parameters.get(f"{i}_b_{sincos}", {}).items():
                base[m, n] = v
            # set boundary perturbation to difference between current and restart
            for (m, n), v in perturbed_parameters.get(f"X{i}_b_{sincos}", {}).items():
                v = v - base_parameters.get(f"X{i}_b_{sincos}", {}).get((m, n), 0)
                if v != 0.0:
                    perturbed[m, n] = v
            if base or perturbed:
                new_base[f"X{i}_b_{sincos}"] = base
                new_perturbed[f"X{i}pert_b_{sincos}"] = perturbed
    return new_base, new_perturbed


def flip_boundary_theta(parameters: MutableMapping) -> MutableMapping:
    """Flip the boundary parameters in the poloidal direction. θ → -θ."""
    output_params = deepcopy(parameters)
    for var in ["X1_b", "X2_b"]:
        if f"{var}_cos" in parameters:
            output_params[f"{var}_cos"] = {}
            for (m, n), value in parameters[f"{var}_cos"].items():
                if m == 0:
                    output_params[f"{var}_cos"][m, n] = value
                else:
                    output_params[f"{var}_cos"][m, -n] = value
        if f"{var}_sin" in parameters:
            output_params[f"{var}_sin"] = {}
            for (m, n), value in parameters[f"{var}_sin"].items():
                if m == 0:
                    output_params[f"{var}_sin"][m, n] = value
                else:
                    output_params[f"{var}_sin"][m, -n] = -value
    return output_params


def flip_boundary_zeta(parameters: MutableMapping) -> MutableMapping:
    output_params = deepcopy(parameters)
    for var in ["X1_b", "X2_b", "X1_a", "X2_a"]:
        if f"{var}_cos" in parameters:
            output_params[f"{var}_cos"] = {}
            for (m, n), value in parameters[f"{var}_cos"].items():
                if m == 0:
                    output_params[f"{var}_cos"][m, n] = value
                else:
                    output_params[f"{var}_cos"][m, -n] = value
        if f"{var}_sin" in parameters:
            output_params[f"{var}_sin"] = {}
            for (m, n), value in parameters[f"{var}_sin"].items():
                if m == 0:
                    output_params[f"{var}_sin"][m, n] = -value
                else:
                    output_params[f"{var}_sin"][m, -n] = value
    return output_params


def shift_boundary_theta_pi(parameters: MutableMapping) -> MutableMapping:
    """
    Shift the theta origin of the boundary by pi.

    cos(m(θ+π)-nζ) = (-1)^m cos(mθ-nζ)
    sin(m(θ+π)-nζ) = (-1)^m sin(mθ-nζ)
    """
    parameters = deepcopy(parameters)
    for var in ["X1", "X2"]:
        for sc in ["cos", "sin"]:
            key = f"{var}_b_{sc}"
            if key in parameters:
                for (m, n), value in list(parameters[key].items()):
                    if m % 2 == 1:
                        parameters[key][m, n] = -value
    return parameters


def flip_parameters_theta(parameters: MutableMapping) -> MutableMapping:
    parameters = flip_boundary_theta(parameters)

    for profile in ["iota"]:
        if profile in parameters:
            parameters[profile]["scale"] = -parameters[profile].get("scale", 1.0)

    return parameters


def flip_parameters_zeta(parameters: MutableMapping) -> MutableMapping:
    parameters = flip_boundary_zeta(parameters)

    if "phiedge" in parameters:
        parameters["phiedge"] = -parameters["phiedge"]
    for profile in ["iota", "I_tor"]:
        if profile in parameters:
            parameters[profile]["scale"] = -parameters[profile].get("scale", 1.0)

    return parameters


def axis_from_boundary(parameters: MutableMapping) -> MutableMapping:
    parameters2 = deepcopy(parameters)
    N = parameters["X1_mn_max"][1]
    parameters2["X1_a_cos"] = {(0, n): parameters["X1_b_cos"][0, n] for n in range(N + 1)}
    parameters2["X2_a_sin"] = {(0, n): parameters["X2_b_sin"][0, n] for n in range(N + 1)}
    if "X1_b_sin" in parameters:
        parameters2["X1_a_sin"] = {(0, n): parameters["X1_b_sin"][0, n] for n in range(N + 1)}
    if "X2_b_cos" in parameters:
        parameters2["X2_a_cos"] = {(0, n): parameters["X2_b_cos"][0, n] for n in range(N + 1)}
    return parameters2


def ellipse_circumference_factor(epsilon: float) -> float:
    """
    Compute the circumference factor of an ellipse with elongation epsilon.
    This uses the approximation by Ramanujan, accurate up to h^5 (6.5% error at ε=10).

    A = a b π = aeff^2 π
    ε = a / b >= 1
    C = 2 π aeff Cf(ε)
    Cf ~ (1 + ε) / (2 √ε) [1 + 3 h / (10 + √(4 - 3 h))]
    h = (ε - 1)^2 / (ε + 1)^2
    """
    h = (epsilon - 1) ** 2 / (epsilon + 1) ** 2
    Cf = (1 + epsilon) / (2 * np.sqrt(epsilon)) * (1 + 3 * h / (10 + np.sqrt(4 - 3 * h)))
    return Cf


def boundary_generator_cases():
    return {
        "ellip_cyl": "elliptic cross-section with no change in zeta",
        "ellip_cyl_breathe": "elliptic cross-section with cross-section area changing with zeta",
        "ellip_cyl_rot": "elliptic cross section of constant ellipticity that only rotates with zeta",
        "ellip_cyl_rot2": "cross section of constant ellipticity that only rotates with zeta, but theta=0 origin remains on positive X1 direction",
        "ellip_cyl_stretch": "elliptic cross section of  changing ellipticity with zeta, not orientation of theta=0",
        "ellip_cyl_helix": "cross section of constant ellipticity, where the axis/boundary moves like a helix over zeta",
        "ellip_cyl_helix_rot": "cross section of constant ellipticity, where the axis/boundary moves like a helix, and rotates along it, over zeta",
        "ellip_cyl_helix_rot2": "cross section of constant ellipticity, where the axis/boundary moves like a helix, and rotates along it, over zeta, but theta=0 origin remains on positive X1 direction",
    }


def boundary_generator(case: str, X1_00=1.0, a0=0.5, ellipticity=0.4, helix_r=0.5):
    """
    Define parameters for some simple boundaries for testing.

    Parameters
    ----------
    case : str
        The name of the boundary:
            see `boundary_generator_cases` dictionary

    X1_00 : float, optional
        =major radius if $X^1=R$
    a0 : float, optional
        =minor radius scale
    ellipticity : float, optional
        =ellipticity of the cross section

    Returns
    -------
    parameter dictionary describing $X^1$ and $X^2$
    """
    params = {}
    match case:
        case "ellip_cyl":
            params["X1_mn_max"] = [1, 0]
            params["X2_mn_max"] = [1, 0]
            params["X1_a_cos"] = {(0, 0): X1_00}
            params["X1_b_cos"] = {
                (0, 0): X1_00,
                (1, 0): a0 * (1.0 - ellipticity),
            }
            params["X2_b_sin"] = {
                (1, 0): a0 * (1.0 + ellipticity),
            }
        case "ellip_cyl_breathe":
            breathe = 0.1
            params["X1_mn_max"] = [1, 1]
            params["X2_mn_max"] = [1, 1]
            params["X1_a_cos"] = {(0, 0): X1_00}
            params["X1_b_cos"] = {
                (0, 0): X1_00,
                (1, 0): a0 * (1.0 - ellipticity) * (1 + breathe),
                (1, -1): -a0 * (1.0 - ellipticity) * 0.5 * breathe,
                (1, 1): -a0 * (1.0 - ellipticity) * 0.5 * breathe,
            }
            params["X2_b_sin"] = {
                (1, 0): a0 * (1.0 + ellipticity) * (1 + breathe),
                (1, -1): -a0 * (1.0 + ellipticity) * 0.5 * breathe,
                (1, 1): -a0 * (1.0 + ellipticity) * 0.5 * breathe,
            }
        case "ellip_cyl_rot":
            params["X1_mn_max"] = [1, 1]
            params["X2_mn_max"] = [1, 1]
            params["X1_a_cos"] = {(0, 0): X1_00}
            params["X1_b_cos"] = {
                (0, 0): X1_00,
                (1, 1): a0,
                (1, -1): -a0 * ellipticity,
            }
            params["X2_b_sin"] = {
                (1, 1): a0,
                (1, -1): a0 * ellipticity,
            }
        case "ellip_cyl_rot2":
            params["X1_mn_max"] = [1, 2]
            params["X2_mn_max"] = [1, 2]
            params["X1_a_cos"] = {(0, 0): X1_00}
            params["X1_b_cos"] = {
                (0, 0): X1_00,
                (1, 0): a0,
                (1, -2): -a0 * ellipticity,
            }
            params["X2_b_sin"] = {
                (1, 0): a0,
                (1, -2): a0 * ellipticity,
            }
        case "ellip_cyl_stretch":
            params["X1_mn_max"] = [1, 1]
            params["X2_mn_max"] = [1, 1]
            params["X1_a_cos"] = {(0, 0): X1_00}
            params["X1_b_cos"] = {
                (0, 0): X1_00,
                (1, 0): a0,
                (1, 1): -0.5 * a0 * ellipticity,
                (1, -1): -0.5 * a0 * ellipticity,
            }
            params["X2_b_sin"] = {
                (1, 0): a0,
                (1, 1): 0.5 * a0 * ellipticity,
                (1, -1): 0.5 * a0 * ellipticity,
            }

        case "ellip_cyl_helix":
            params["X1_mn_max"] = [1, 1]
            params["X2_mn_max"] = [1, 1]
            params["X1_a_cos"] = {
                (0, 0): X1_00,
                (0, 1): helix_r,
                (0, -1): helix_r,
            }
            params["X2_a_sin"] = {
                (0, 1): helix_r,
                (0, -1): helix_r,
            }
            params["X1_b_cos"] = {
                (0, 0): X1_00,
                (0, 1): helix_r,
                (0, -1): helix_r,
                (1, 0): a0 * (1.0 - ellipticity),
            }
            params["X2_b_sin"] = {
                (1, 0): a0 * (1.0 + ellipticity),
                (0, 1): helix_r,
                (0, -1): helix_r,
            }
        case "ellip_cyl_helix_rot":
            params["X1_mn_max"] = [1, 1]
            params["X2_mn_max"] = [1, 1]
            params["X1_a_cos"] = {
                (0, 0): X1_00,
                (0, 1): helix_r,
                (0, -1): helix_r,
            }
            params["X2_a_sin"] = {
                (0, 1): helix_r,
                (0, -1): helix_r,
            }
            params["X1_b_cos"] = {
                (0, 0): X1_00,
                (0, 1): helix_r,
                (0, -1): helix_r,
                (1, 1): a0,
                (1, -1): -a0 * ellipticity,
            }
            params["X2_b_sin"] = {
                (1, 1): a0,
                (1, -1): a0 * ellipticity,
                (0, 1): helix_r,
                (0, -1): helix_r,
            }
        case "ellip_cyl_helix_rot2":
            params["X1_mn_max"] = [1, 2]
            params["X2_mn_max"] = [1, 2]
            params["X1_a_cos"] = {
                (0, 0): X1_00,
                (0, 1): helix_r,
                (0, -1): helix_r,
            }
            params["X2_a_sin"] = {
                (0, 1): helix_r,
                (0, -1): helix_r,
            }
            params["X1_b_cos"] = {
                (0, 0): X1_00,
                (0, 1): helix_r,
                (0, -1): helix_r,
                (1, 0): a0,
                (1, -2): -a0 * ellipticity,
            }
            params["X2_b_sin"] = {
                (1, 1): a0,
                (1, -1): a0 * ellipticity,
                (0, 0): helix_r,
                (0, -2): helix_r,
            }
        case _:
            raise ValueError(f"request boundary '{case}', does not exist!")

    return params


def linking_number(curve_a: np.ndarray, curve_b: np.ndarray, tol=1e-15, endpoint=False):
    r"""
    Compute the linking number of two non-intersecting curves $C_a(\zeta_a),C_b(\zeta_b)$, solving the (non-singular) double integral over two curves

    $$
    \text{Lk} = \frac{1}{4\pi}\int_0^{2\pi} \int_0^{2\pi} \frac{ r_b - r_a}{\abs{r_b - r_a}^3} \cdot \left(\frac{\partial r_a}{\partial\zeta_a}  \times \frac{\partial r_b}{\partial\zeta_b}\right) d{\zeta_a} d{\zeta_b}
    $$

    Here, the double integral is approximated by representing the curves as polygons, and computing the linking number as the sum of solid angles between all linear segments of polygon curves. The solid angle of two linear segments can be computed exactly. See paper by  [Klenin and Langowski](https://onlinelibrary.wiley.com/doi/10.1002/1097-0282(20001015)54:5%3C307::AID-BIP20%3E3.0.CO;2-Y)
    See GVEC documentation for more details.

    Parameters
    ----------
    curve_a : np.ndarray
        the x,y,z point positions of the first curve. First and last point must coincide, if `endpoint=True`. shape must be [npoints_a,3]
    curve_b : np.ndarray
        the x,y,z point positions of the first curve. First and last point must coincide, if `endpoint=True`. shape must be [npoints_b,3]
    tol : float
        tolerance to consider points coincident (default: 1e-15)
    endpoint : bool
        `True`:  the first and last point of each curve coincide, else the last point is connected to the first point to close the curve. Default is `False`.

    Returns
    -------
    Lk : float
        the linking number of the two curves
    """
    if not curve_a.ndim == 2:
        raise ValueError("curve_a must be a 2D array")
    if not curve_b.ndim == 2:
        raise ValueError("curve_b must be a 2D array")
    if not curve_a.shape[1] == 3:
        raise ValueError("second dimension of curve_a must be of size 3")
    if not curve_b.shape[1] == 3:
        raise ValueError("second dimension of curve_b must be of size 3")
    closed_a = np.allclose(curve_a[0, :], curve_a[-1, :])
    closed_b = np.allclose(curve_b[0, :], curve_b[-1, :])
    if endpoint:
        if not closed_a:
            raise ValueError("first and last point of curve_a must coincide (closed curve)")
        if not closed_b:
            raise ValueError("first and last point of curve_b must coincide (closed curve)")
        _curve_a = curve_a
        _curve_b = curve_b
    else:
        if closed_a:
            raise ValueError(
                "first and last point of curve_a coincide, but endpoint=False was chosen"
            )
        if closed_b:
            raise ValueError(
                "first and last point of curve_b coincide, but endpoint=False was chosen"
            )
        _curve_a = np.vstack([curve_a, curve_a[0, :]])
        _curve_b = np.vstack([curve_b, curve_b[0, :]])
    # nseg_a = curve_a.shape[0]-1
    # nseg_b = curve_b.shape[0]-1
    # Lk = 0.0
    # for i in range(nseg_a):
    #    for j in range(nseg_b):
    #        r1 = curve_a[i,:]
    #        r2 = curve_a[i+1,:]
    #        r3 = curve_b[j,:]
    #        r4 = curve_b[j+1,:]
    #        Lk += solid_angle_between_segments(r1,r2,r3,r4)
    # Lk /= (4.0*np.pi)
    return (
        2
        * np.sum(
            solid_angle_between_segments(
                _curve_a[0:-1, None, :],
                _curve_a[1:, None, :],
                _curve_b[None, 0:-1, :],
                _curve_b[None, 1:, :],
                tol=tol,
            )
        )
        / (4 * np.pi)
    )


def writhe_from_polygon(curve: np.ndarray, endpoint=False):
    r"""
    Compute the writhe of a closed curve $C(\zeta)$, solving the (singular!) double integral over the curve

    $$
    \text{Wr} = \frac{1}{4\pi}\int_0^{2\pi} \int_0^{2\pi} \frac{ r(\zeta^\prime) - r(\zeta)}{\abs{r(\zeta^\prime) - r(\zeta)}^3} \cdot \left(\frac{\partial r}{\partial\zeta}  \times \frac{\partial r}{\partial\zeta^\prime} \right) d{\zeta} d{\zeta^\prime}
    $$

    Here, the double integral is approximated by representing the curve as a polygon, and computing the writhe as the sum of solid angles between all linear segments of polygon curve. The solid angle of two linear segments can be computed exactly. See paper by  [Klenin and Langowski](https://onlinelibrary.wiley.com/doi/10.1002/1097-0282(20001015)54:5%3C307::AID-BIP20%3E3.0.CO;2-Y)
    See GVEC documentation for more details.

    Parameters
    ----------
    curve : np.ndarray
        the x,y,z point positions of the curve. First and last point must coincide, if `endpoint=True`. shape must be [npoints,3]
    endpoint : bool
        `True`: the first and last point of the curve coincide, else the last point is connected to the first point to close the curve. Default is `False`

    Returns
    -------
    Wr : float
        the approximate writhe of the curve

    Warnings
    --------
    The algorithm converges very slowly with the number of line segments.
    """
    if not curve.ndim == 2:
        raise ValueError("curve must be a 2D array")
    if not curve.shape[1] == 3:
        raise ValueError("second dimension of curve must be of size 3")
    if endpoint:
        if not np.allclose(curve[0, :], curve[-1, :]):
            raise ValueError("first and last point of curve must coincide (closed curve)")
        _curve = curve
    else:
        _curve = np.vstack([curve, curve[0, :]])
    nseg = _curve.shape[0] - 1
    # mask to select only pairs of segments (i,j) with j>i+1 to avoid coincident or adjacent segments
    mask = np.triu(np.ones((nseg, nseg), dtype=np.int8), k=2)
    # apply mask to segment(i) x segment(j), with j>i+1 remaining
    r1 = (_curve[0:-1, None, :] * mask[:, :, None])[(mask == 1), :]
    r2 = (_curve[1:, None, :] * mask[:, :, None])[(mask == 1), :]
    r3 = (_curve[None, 0:-1, :] * mask[:, :, None])[(mask == 1), :]
    r4 = (_curve[None, 1:, :] * mask[:, :, None])[(mask == 1), :]

    # factor 2 to account for symmetry of (i,j) and (j,i)
    return 2 * np.sum(solid_angle_between_segments(r1, r2, r3, r4)) / (2 * np.pi)


def solid_angle_between_segments(
    r1: np.ndarray, r2: np.ndarray, r3: np.ndarray, r4: np.ndarray, tol=1e-15
):
    """
    Compute the solid angle between two line segments r1->r2 and r3->r4.
    See paper by  [Klenin and Langowski](https://onlinelibrary.wiley.com/doi/10.1002/1097-0282(20001015)54:5%3C307::AID-BIP20%3E3.0.CO;2-Y) and direct summation from https://doi.org/10.1145/3450626.3459778
    r1,r2,r3,r4 must have the same shape, with last dimension of size 3 (x,y,z)
    Args:
        r1 (np.ndarray): last dimension is x,y,z position of the first  point of the first segment
        r2 (np.ndarray): last dimension is x,y,z position of the second point of the first segment
        r3 (np.ndarray): last dimension is x,y,z position of the first  point of the second segment
        r4 (np.ndarray): last dimension is x,y,z position of the second point of the second segment
        tol (float): tolerance to consider points coincident (default: 1e-15)
    Returns:
        omega (float): solid angle between the two segments, without 1/(4pi) factor
    """
    # # algorithm with 4 arctan
    # N1 = np.cross(r3 - r1, r4 - r1)
    # N2 = np.cross(r4 - r1, r4 - r2)
    # N3 = np.cross(r4 - r2, r3 - r2)
    # N4 = np.cross(r3 - r2, r3 - r1)
    # omega = (
    #     np.arctan2(np.vecdot(N1, N2), np.sqrt(np.sum(np.cross(N1, N2) ** 2, axis=-1)))
    #     + np.arctan2(np.vecdot(N2, N3), np.sqrt(np.sum(np.cross(N2, N3) ** 2, axis=-1)))
    #     + np.arctan2(np.vecdot(N3, N4), np.sqrt(np.sum(np.cross(N3, N4) ** 2, axis=-1)))
    #     + np.arctan2(np.vecdot(N4, N1), np.sqrt(np.sum(np.cross(N4, N1) ** 2, axis=-1)))
    # )
    # sigma = np.vecdot(np.cross(r4 - r3, r2 - r1), (r3 - r1))
    # #set very small sigma to zero to avoid  sign(+0)->1 and sign(-0)->-1 issues
    # sigma[np.abs(sigma) < tol] = 0
    # return omega * np.sign(sigma)

    # faster algorithm for direct summation, with two atan, from https://doi.org/10.1145/3450626.3459778
    # where they define a= (r3-r1) , b= (r3-r2) , c= (r4-r2), d= (r4-r1)
    # omega = 2( atan(b.(c x a) / (|a||c||b| +c.a|b| + b.a|c| + b.c|a| ) )
    #           -atan(d.(c x a) / (|a||c||d| +c.a|d| + d.a|c| + d.c|a| ) ) )

    cxa = np.cross((r4 - r2), (r3 - r1))
    abs_a = np.sqrt(np.vecdot((r3 - r1), (r3 - r1)))
    abs_c = np.sqrt(np.vecdot((r4 - r2), (r4 - r2)))
    ac2 = abs_a * abs_c + np.vecdot((r3 - r1), (r4 - r2))
    ac2_vec = (r3 - r1) * np.expand_dims(abs_c, axis=-1) + (r4 - r2) * np.expand_dims(
        abs_a, axis=-1
    )

    omega_half = np.arctan2(
        np.vecdot((r3 - r2), cxa),
        ac2 * np.sqrt(np.vecdot((r3 - r2), (r3 - r2))) + np.vecdot((r3 - r2), ac2_vec),
    ) - np.arctan2(
        np.vecdot((r4 - r1), cxa),
        ac2 * np.sqrt(np.vecdot((r4 - r1), (r4 - r1))) + np.vecdot((r4 - r1), ac2_vec),
    )
    return omega_half


def scale_boundary(
    parameters: MutableMapping, scale: float, radius: Literal["major", "minor", "both"] = "both"
) -> MutableMapping:
    """
    Scale the boundary parameters by a given factor, with the 'radius' argument specifying which modes to scale:
    - "major": scale only the m=0 modes, which approximately corresponds to scaling the effective major radius
    - "minor": scale only the m!=0 modes, which corresponds to scaling the effective minor radius
    - "both": scale all modes, which corresponds to scaling the overall size of the configuration

    Note
    ----
    This does not scale the reference frame (hmap), which may be required for more advanced reference frames like the GFrame.
    """
    parameters = deepcopy(parameters)
    for var in ["X1_b", "X2_b"]:
        for sc in ["cos", "sin"]:
            key = f"{var}_{sc}"
            if key in parameters:
                for (m, n), value in parameters[key].items():
                    if m == 0 and radius in ["major", "both"]:
                        parameters[key][m, n] = value * scale
                    elif m != 0 and radius in ["minor", "both"]:
                        parameters[key][m, n] = value * scale
    return parameters
