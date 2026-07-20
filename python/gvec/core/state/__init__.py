# Copyright (c) 2026 GVEC Contributors, Max Planck Institute for Plasma Physics
# License: MIT
"""gvec.state - low-level python API for postprocessing

The gvec.state module provides the `State` class which wraps the src/pygvec/state.f90 fortran module.
It provides an encapsulation for the state of the GVEC library, allowing safe postprocessing routines.
This module checks the input arguments and handles the initialization and finalization of the Fortran state.

The `State` class is distributed across multiple files, grouping the methods and properties by their functionality.
"""

from ._class import State
from .util import load_state, find_state, find_states
