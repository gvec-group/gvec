# Copyright (c) 2025 GVEC Contributors, Max Planck Institute for Plasma Physics
# License: MIT

import os
import sys
sys.path.append(os.getenv('GVEC_TO_GENE_MODULE_PATH'))

import gvec_to_gene
from cffi import FFI
import numpy as np

ffi = FFI()

#---------------------------------------------------------------------
# test 1

name_val = 'esta_caso.vtu'
p = ffi.new("char[]", name_val.encode())


#p = ffi.new("char[]", b"here should be a file name :)")
gvec_to_gene.test_print_file_name(p)

#---------------------------------------------------------------------
# test 2

strl = 16
arrc_type = "S{}".format(strl)
arrc_np = np.array(['arg-1\x00blaa', 'arg--2\x00blaa', 'arg---3\x00blaa'],
                   dtype=np.dtype(arrc_type))
arrc = ffi.cast("char(*)[]", arrc_np.ctypes.data)
gvec_to_gene.test_print_char_rank2_array(strl,len(arrc_np),arrc)

#---------------------------------------------------------------------
# test 3

nthet = 2
nzeta = 5

arr_in = np.zeros((nthet,nzeta),dtype=np.double,order='F')
arr_out = np.zeros((nthet,nzeta),dtype=np.double,order='F')
arr_out_e = np.array([i for i in range(1,nthet*nzeta+1)],
                     dtype=np.double).reshape((nthet,nzeta),order='F')

# get c pointers to arrays
arr_in_c = ffi.cast("double(*)[]", arr_in.ctypes.data)
arr_out_c = ffi.cast("double(*)[]", arr_out.ctypes.data)

gvec_to_gene.test_pass_arrays_shift(nthet,nzeta,arr_in_c,arr_out_c)

if(np.array_equal(arr_out_e, arr_out)):
    sys.exit(0)
else:
    print("arr_in:\n" + str(arr_in))
    print("arr_out:\n" + str(arr_out))
    print("arr_out_e:\n" + str(arr_out_e))
    sys.exit(1)
