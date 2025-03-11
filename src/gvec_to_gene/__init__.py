# Copyright (c) 2025 GVEC Contributors, Max Planck Institute for Plasma Physics
# License: MIT

from subprocess import check_output
from cffi import FFI
import os
import sys
from configparser import ConfigParser
from pathlib import Path


def get_lib_handle(definitions, header_file, library_file):
    ffi = FFI()
    command = ['cc', '-E'] + definitions + [header_file]
    interface = check_output(command).decode('utf-8')

    # remove possible \r characters on windows which
    # would confuse cdef
    _interface = [l.strip('\r') for l in interface.split('\n')]

    ffi.cdef('\n'.join(_interface))
    lib = ffi.dlopen(library_file, ffi.RTLD_DEEPBIND | ffi.RTLD_GLOBAL | ffi.RTLD_LAZY)
    return lib


# this interface requires the header file and library file
# and these can be either provided by interface_file_names.cfg
# in the same path as this file
# or if this is not found then using environment variables
_this_path = Path(os.path.dirname(os.path.realpath(__file__)))
_cfg_file = _this_path / 'interface_file_names.cfg'
if _cfg_file.exists():
    config = ConfigParser()
    config.read(_cfg_file)
    header_file_name = config.get('configuration', 'header_file_name')
    _header_file = _this_path / header_file_name
    _header_file = str(_header_file)
    library_file_name = config.get('configuration', 'library_file_name')
    _library_file = _this_path / library_file_name
    _library_file = str(_library_file)
else:
    _header_file = os.getenv('GVEC_TO_GENE_HEADER_FILE')
    assert _header_file is not None
    _library_file = os.getenv('GVEC_TO_GENE_LIBRARY_FILE')
    assert _library_file is not None

_lib = get_lib_handle(definitions=['-DGVEC_TO_GENE_API=', '-DGVEC_TO_GENE_NOINCLUDE'],
                      header_file=_header_file,
                      library_file=_library_file)


# we change names to obtain a more pythonic API
init = _lib.init_gvec_to_gene
coords = _lib.gvec_to_gene_coords
finalize = _lib.finalize_gvec_to_gene
write_data_to_vtk = _lib.write_data_to_vtk
test_print_char_rank2_array = _lib.test_print_char_rank2_array
test_print_file_name = _lib.test_print_file_name
test_pass_arrays_shift = _lib.test_pass_arrays_shift
test_int_array = _lib.test_int_array

__all__ = [
    '__version__',
    'init',
    'coords',
    'finalize',
    'write_data_to_vtk',
    'test_print_char_rank2_array',
    'test_print_file_name',
    'test_pass_arrays_shift',
    'test_int_array',
]
