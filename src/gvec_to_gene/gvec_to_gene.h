// Copyright (c) 2025 GVEC Contributors, Max Planck Institute for Plasma Physics
// License: MIT

/* CFFI would issue warning with pragma once */
#ifndef GVEC_TO_GENE_H_INCLUDED
#define GVEC_TO_GENE_H_INCLUDED

#ifndef GVEC_TO_GENE_API
#include "gvec_to_gene_export.h"
#define GVEC_TO_GENE_API GVEC_TO_GENE_EXPORT
#endif

// cffi doesn't like 3rd c-arrays :(
GVEC_TO_GENE_API
void write_data_to_vtk(const int dim1,
		       const int vecdim,
		       const int nVal,
		       const int NPlot[], // [dim1]
		       const int nElems,
		       const int strlen,
		       const char VarNames_c[][], // [strlen][nVal]
		       const double Coord[][], // [vecdim][1:PRODUCT(Nplot+1)][nElems]
		       const double Values[][], // [nVal][1:PRODUCT(Nplot+1)][nElems]
		       char FileString_c[]);


GVEC_TO_GENE_API
void init_gvec_to_gene(char fileName[]);

// cffi doesn't like 3rd c-arrays :(
GVEC_TO_GENE_API
void gvec_to_gene_coords(const int nthet,
			 const int nzeta,
			 const double spos_in,
			 const double theta_star_in[][], // [nthet][nzeta]
			 const double zeta_in[][], // [nthet][nzeta]
			 double theta_out[][], // [nthet][nzeta]
			 double cart_coords[][]); // [3][nthet][nzeta]

GVEC_TO_GENE_API
void finalize_gvec_to_gene();

GVEC_TO_GENE_API
void test_print_char_rank2_array(const int strlen, const int nval,
				 const char varnames_c[][]);

GVEC_TO_GENE_API
void test_print_file_name(char fileName[]);

GVEC_TO_GENE_API
void test_pass_arrays_shift(const int nthet, const int nzeta,
			    const double arr_in[][],
			    double arr_out[][]);

GVEC_TO_GENE_API
void test_int_array(const int dim1,
			    const int Nplot[]); //[dim1]);

#endif /* GVEC_TO_GENE_H_INCLUDED */
