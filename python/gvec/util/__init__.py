"""GVEC utility module

This module is part of the gvec python package, but also used directly in the tests.
"""

from .util import (
    chdir,
    get_compile_options,
    logging_setup,
    version_info,
    compute_FD,
)
from .geometry import (
    axis_from_boundary,  # Needs docstrings
    boundary_generator,
    boundary_generator_cases,
    check_boundary_direction,
    compute_boundary_perturbation,  # Needs docstrings
    effective_minor_radius,  # Needs docstrings
    ellipse_circumference_factor,
    evaluate_axis,
    evaluate_boundary,
    flip_boundary_theta,  # Needs docstrings
    flip_boundary_zeta,  # Needs docstrings
    flip_parameters_theta,  # Needs docstrings
    flip_parameters_zeta,  # Needs docstrings
    linking_number,
    shift_boundary_theta_pi,
    signed_cross_sectional_area,  # Needs docstrings
    solid_angle_between_segments,
    writhe_from_polygon,
)
from .params import (
    CaseInsensitiveDict,
    adapt_parameter_file,
    bspl2gvec,  # should probably be moved
    flatten_parameters,
    get_boundary_from_statefile,
    parameters_from_vmec,
    read_parameter_file_ini,
    read_parameters,
    stack_parameters,
    stringify_mn_parameters,
    unstringify_mn_parameters,
    write_parameter_file_ini,
    write_parameters,
)
