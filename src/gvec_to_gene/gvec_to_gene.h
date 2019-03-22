/* CFFI would issue warning with pragma once */
#ifndef GVEC_TO_GENE_H_INCLUDED
#define GVEC_TO_GENE_H_INCLUDED

#ifndef GVEC_TO_GENE_API
#include "gvec_to_gene_export.h"
#define GVEC_TO_GENE_API GVEC_TO_GENE_EXPORT
#endif

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
void test_print_file_name(char fileName[]);

GVEC_TO_GENE_API
void test_pass_arrays_shift(const int nthet, const int nzeta,
			   const double arr_in[][],
			   double arr_out[][]);

#endif /* GVEC_TO_GENE_H_INCLUDED */
