# Copyright (c) 2026 GVEC Contributors, Max Planck Institute for Plasma Physics
# License: MIT
"""gvec.lib - intermediate python API for the python-fortran interface

The python-fortran interface uses a hierarchy of layers to connect the fortran library to the python package with varying degrees of abstraction and safety.
1) fortran library (src/*)
2) fortran wrapping layer (gen. by f90wrap)
3) C wrapping layer (gen. by f2py)
4) python wrapping layer (gen. by f90wrap: _lib.py)
5) final wrapping layer (lib.py)
6) low-level python API (core.state, core.run)
7) high-level python API (core.compute, quantities)
"""

from . import _lib
from gvec._lib import modgvec_py_binding as _binding
from gvec._lib import modgvec_py_run as _run
from gvec._lib import modgvec_py_state as _state
from gvec._lib import modgvec_globals as _globals
from gvec._lib import modgvec_sfl_boozer as _sfl_boozer
from gvec._lib import modgvec_biotsavart as _biotsavart

t_sfl_boozer = _lib.Modgvec_Sfl_Boozer.t_sfl_boozer

c_rProfile = _lib.Modgvec_Rprofile_Base.c_rProfile
t_rProfile_poly = _lib.Modgvec_Rprofile_Poly.t_rProfile_poly
t_rProfile_bspl = _lib.Modgvec_Rprofile_Bspl.t_rProfile_bspl

t_fBase = _lib.Modgvec_Fbase.t_fBase
