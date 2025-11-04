# Copyright (c) 2025 GVEC Contributors, Max Planck Institute for Plasma Physics
# License: MIT
"""GVEC exception classes

These exceptions are used throughout GVEC, in particular the Fortran part, to
signal errors back to calling library or user. With f90wrap, we always initially
get a RuntimeError, which is then converted to a more specific exception type here.
"""

import sys
import contextlib
import re
import logging
import functools

# === Exception Classes === #


class GVECError(Exception):
    """Base class for all GVEC exceptions."""

    pass


class FortranError(GVECError):
    """General exception raised from Fortran code."""

    pass


class MissingParameterError(GVECError, KeyError):
    """Exception raised for missing input parameters which are required."""

    pass


class InvalidParameterError(GVECError, ValueError):
    """Exception raised for invalid input parameters."""

    pass


class InitializationError(GVECError):
    """Exception raised during initialization of GVEC."""

    pass


class UnknownError(GVECError):
    """Exception raised for unknown errors."""

    pass


ERROR_MAP = dict(
    RuntimeError=FortranError,
    MissingParameterError=MissingParameterError,
    InvalidParameterError=InvalidParameterError,
    FileNotFoundError=FileNotFoundError,
    InitializationError=InitializationError,
)
ERROR_PATTERN = r"(.+) \| (" + "|".join(ERROR_MAP.keys()) + r") \| (.+)"

# === Exception Handler === #


@contextlib.contextmanager
def catch_gvec_errors():
    try:
        yield
    except RuntimeError as e:
        if match := re.match(ERROR_PATTERN, str(e)):
            active_region = match.group(1)
            error_type = match.group(2)
            error_msg = match.group(3)

            logger = logging.getLogger(f"gvec.lib.{active_region}")
            logger.error(f"{error_type}: {error_msg.strip()}")

            if error_type in ERROR_MAP:
                error = ERROR_MAP[error_type](error_msg)
            else:
                error = UnknownError(f"{error_type}: {error_msg}")
            error.logged = True
            error.active_region = active_region
            raise error from None
        else:
            raise


@contextlib.contextmanager
def log_errors(logger=None, exit=False):
    try:
        yield
    except Exception as error:
        if not getattr(error, "logged", False):
            if logger is None:
                logger = logging.getLogger("gvec")
            logger.error(f"{type(error).__name__}: {str(error).strip()}")
        if exit:
            sys.exit(1)
        raise


def without_traceback(func):
    """Decorator to suppress tracebacks and log errors and exit instead."""

    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        with log_errors(exit=True):
            return func(*args, **kwargs)

    return wrapper


def with_log_errors(func):
    """Decorator to log errors."""

    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        with log_errors():
            return func(*args, **kwargs)

    return wrapper
