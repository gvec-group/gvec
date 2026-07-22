# Copyright (c) 2026 GVEC Contributors, Max Planck Institute for Plasma Physics
# License: MIT
"""gvec.core.state._binding - utilities for handling the python-fortran interface for gvec.State"""

# === Imports === #

import functools
import logging
import tempfile

import gvec.util
from gvec.errors import catch_gvec_errors
from gvec.lib import _binding, _run, _state, _globals, t_sfl_boozer

# === Globals === #

logger = logging.getLogger("gvec.state")

# this variable tracks the global state of the fortran library
# it tracks which State object is currently bound/loaded/initialized
# _post.initialized tracks additionally if the State was properly initialized
bound_state = None

# === Decorators & Utility Functions === #


flush_stdout = _binding.flush_stdout


def with_binding(method):
    """This decorator wraps a method to ensure the State object is bound
    and the fortran library is initialized.

    If another state is bound, it is first unbound and the current state is bound.
    This ensures that unnecessary re-initialization of the fortran library is avoided,
    but switching between different State objects is still possible with minimal (syntactic) overhead.
    """

    @functools.wraps(method)
    def wrapped(self, *args, **kwargs):
        global bound_state

        if bound_state is not self:
            self.bind(force=True)

        try:
            with catch_gvec_errors():
                return method(self, *args, **kwargs)
        except:
            self.unbind(cleanup=True)
            raise

    return wrapped


# === Binding to the Fortran library === #


def bind(self, force: bool = False):
    """Bind this State object to the Fortran library. Allocate & initialize everything."""
    global bound_state
    if bound_state is self:
        return
    if bound_state is not None:
        if not force:
            raise RuntimeError(
                f"Another state {bound_state!r} is already bound to the fortran library. Unbind that first."
            )
        else:
            bound_state.unbind()

    if _run.initialized:
        logger.debug("Fortran library still initialized by run, attempting cleanup.")
        _run.cleanup()
    if _state.initialized:
        raise RuntimeError("Fortran library is initialized, but no state is bound!?")

    bound_state = self
    logger.debug(f"Binding state {self!r} to the fortran library.")

    _binding.redirect_abort()  # redirect 'abort' to raise a 'RuntimeError'
    if self._stdout is None:
        self._stdout = tempfile.NamedTemporaryFile(mode="r", prefix="gvec-stdout-")
    _binding.redirect_stdout(self._stdout.name)
    logger.debug(f"Redirected stdout to {self._stdout.name}")

    MAXLEN = _globals.maxlen
    with open(self.parameterfile, "r") as f:
        content = f.readlines()
        if any([len(line) > MAXLEN for line in content]):
            raise ValueError(
                f"Parameter file {self.parameterfile} contains lines longer than {MAXLEN} characters."
            )

    with gvec.util.chdir(self.rundir):
        try:
            with catch_gvec_errors():
                logger.debug("Initialize from parameterfile")
                _state.init(self.parameterfile.name)
                if self.statefile is not None:
                    logger.debug("Read state from statefile")
                    if not self.statefile.is_absolute():
                        statefile = self.statefile.relative_to(self.rundir)
                    else:
                        statefile = self.statefile
                    _state.readstate(statefile)
                else:
                    logger.debug("Initialize solution without statefile")
                    _state.initsolution()
        except:
            logger.debug("Error during binding, attempting cleanup.")
            self.unbind(cleanup=True)
            raise

    self._children = []

    if not _state.initialized:
        raise RuntimeError("Failed to initialize fortran library.")
    logger.debug(f"Bound state {self!r} to the fortran library.")


def unbind(self, cleanup: bool = False):
    """Unbind this State object from the Fortran library. Finalize & deallocate everything."""
    global bound_state

    if bound_state is not self:
        raise RuntimeError(
            f"State {self!r} is not bound to the fortran library, but {bound_state!r} is."
        )
    if not _state.initialized and not cleanup:
        raise RuntimeError("Fortran library is not initialized, but state is bound!?")

    bound_state = None
    logger.debug(f"Unbinding state {self!r} from the fortran library.")

    with catch_gvec_errors():
        for child in self._children:
            if isinstance(child, t_sfl_boozer):
                logger.debug(f"Unbinding child {child!r} from the fortran library.")
                del child
            else:
                logger.error(f"Unknown child: {child!r}")
        self._children = []

        _state.finalize()

    if _state.initialized:
        raise RuntimeError("Failed to finalize fortran library.")
    logger.debug(f"Unbound state {self!r} from the fortran library.")
