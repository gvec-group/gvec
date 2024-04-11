#================================================================================================================================#
# Copyright (C) 2024 Robert Babin <robert.babin@ipp.mpg.de>
#
# This file is part of GVEC. GVEC is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 
# of the License, or (at your option) any later version.
#
# GVEC is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
# of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License v3.0 for more details.
#
# You should have received a copy of the GNU General Public License along with GVEC. If not, see <http://www.gnu.org/licenses/>.
#================================================================================================================================#
"""pygvec postprocessing"""

from ._post import modpygvec_post as _post

from pathlib import Path

_status = "pre"

class State:
    def __init__(self, parameterfile, statefile):
        global _status
        if _status != "pre":
            raise NotImplementedError("Only one instance of State is allowed.")
        _status = "init"
        if not Path(parameterfile).exists():
            raise FileNotFoundError(f"Parameter file {parameterfile} does not exist.")
        if not Path(statefile).exists():
            raise FileNotFoundError(f"State file {statefile} does not exist.")
        
        _post.init(parameterfile)
        _post.readstate(statefile)
    
    def __enter__(self):
        global _status
        if _status != "init":
            raise NotImplementedError("State is not initialized.")
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        global _status
        if _status == "init":
            _post.finalize()
            _status = "final"
    
    def __del__(self):
        global _status
        if _status == "init":
            _post.Finalize()
            _status = "final"

