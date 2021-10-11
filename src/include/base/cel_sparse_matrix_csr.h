/*     
 * File:   cel_sparse_matrix_coo_3d.h
 * Project CRESTA (see details on https://cresta-project.eu) 
 * CRESTA Exascale library
 * Send you email about the bugs to
 * The High Performance Computing Center Stuttgart (HLRS) of the University of Stuttgart
 * Dmitry Khabi,  (email: khabi@hlrs.de)
 * Created on March 25, 2013
*/
#ifndef CEL_SPARSE_MATRIX_CSR_H
#define CEL_SPARSE_MATRIX_CSR_H
#include "cel_types.h"
#include "cel_domain_3d.h"
#include "cel_sparse_matrix_coo.h"
//! \struct Sparse matrix in CSR format
typedef struct 
{
	INDEX xadj_length; //!length of row pointers
	INDEX* xadj;//!array with row pointers
	INDEX* adjncy;//!array with column indices
	INDEX* vtxdist;//!array with distribution of rows between processes
	REAL* values;//!array with values (if is_diagonal_sperate==1 only off diagonal values)
	REAL* diagonal;//!array with diagonal values or nil if is_diagonal_sperate==0
	INDEX* rows_diag_ptrs;//!array with row pointers of diagonal part of the matrix
	INDEX* cols_diag_indices;//!array with column indices of diagonal part of the matrix
	INDEX block_indx;//! block_indx block index of this part of matrix
	BOOLEAN is_diagonal_separat;//! indicate whether the diagonal part of matrix is stored in separe arrays
} cel_sparse_matrix_csr;

//! Initialize new sparse matricies in compresed row storage format(CSR)
//! @param (in) matrix_num number of matrices to initialiez
//! @param (in) output_on 1 - print error to sdterr; 0 - no output
//! @return  pointer to 1 dimensional array with \struct cel_sparse_matrix_coo with nill initialized fields OR nil if no memory
cel_sparse_matrix_csr* cel_sparse_matrix_csr_new(INDEX matrix_num, int output_on);

//! Allocate memory for internal fields in \struct sparse matrix in compresed row storage format(CSR)
//! @param (in) elements_num minimal number of non zero elements in matrix
//! @param (in) rows_num minimal number of rows in matrix
//! @param (in) is_diagonal_separat indicate whether the diagonal part of matrix is stored in separe arrays
//! @param (out)   pointer to \struct cel_sparse_matrix_csr
//! @param (in) output_on 1 - print error to sdterr; 0 - no output
//! @return  0 - no error ; -1 - not enought memory
int cel_sparse_matrix_csr_allocate(INDEX elements_num, INDEX rows_num, BOOLEAN is_diagonal_separat, cel_sparse_matrix_csr* sparse_matrix_csr, int output_on);

//! Free sparse matricies in compresed row storage format (CSR) - only internal structure are realised
//! @param (in) pointer to 1 dimensional array with \struct cel_sparse_matrix_csr with nill initialized fields OR nil if no memory
//! @param (in) matrix_num number of matrices to initialiez
void cel_sparse_matrix_csr_free(cel_sparse_matrix_csr* sparse_matrix_csr, INDEX matrix_num);

//! Output of the sparse matrix in compresed row storage format(CSR) to stderr
//! @param (in) matrix_coo pointer to \struct sparse matrix in compresed row storage format(CSR)
//! @param (in) output_on 1 - print error to sdterr; 0 - no output
void cel_sparse_matrix_csr_print(cel_sparse_matrix_csr* sparse_matrix_csr, int output_on);

//! convert sparse matrix in coo format to sparse matrix in csr format
//! @param (in) sparse_matrix_coo pointer to sparse matrix in coordinate format(COO)
//! @param (out) sparse_matrix_csr pointer to sparse matrix in compresed row storage format(CSR)
//! @param (in) output_on 1 - print error to sdterr; 0 - no output
//! @return  0 - no error ; -1 - not enought memory
int cel_sparse_matrix_scr_from_sparse_matrix_coo(cel_sparse_matrix_coo* sparse_matrix_coo, cel_sparse_matrix_csr* sparse_matrix_csr, int output_on);

//! initialize vtxdist in sparse matrix in compresed row storage format(CSR)
//! @param (in) domain_3D description of regular grid on domain [0;1]x[0;1]x[0;1]
//! @param (in) processes_num number of processes
//! @param (out) sparse_matrix_csr pointer to \struct sparse matrix in compresed row storage format(CSR)
//! @param (in) output_on 1 - print error to sdterr; 0 - no output
//! @return  0 - no error ; -1 - not enought memory
int cel_sparse_matrix_csr_init_vtxdist(cel_domain_3D* domain_3D, INDEX processes_num, cel_sparse_matrix_csr* sparse_matrix_csr, int output_on);

//! get info about the matrix
//! @param (out) sparse_matrix_csr pointer to \struct sparse matrix in compresed row storage format(CSR)
//! @param (out) rows_number number of rows
//! @param (out) elements_number number of non-zero elements
void cel_sparse_matrix_csr_info(cel_sparse_matrix_csr* sparse_matrix_csr, INDEX *rows_number, INDEX *elements_number);

#endif

