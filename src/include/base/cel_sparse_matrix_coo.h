/*     
 * File:   cel_sparse_matrix_coo_3d.h
 * Project CRESTA (see details on https://cresta-project.eu) 
 * CRESTA Exascale library
 * Send you email about the bugs to
 * The High Performance Computing Center Stuttgart (HLRS) of the University of Stuttgart
 * Dmitry Khabi,  (email: khabi@hlrs.de)
 * Created on March 25, 2013
*/
#ifndef CEL_SPARSE_MATRIX_COO_H
#define CEL_SPARSE_MATRIX_COO_H
#include "cel_types.h"

//! \struct Sparse matrix in coordinate format (COO)
typedef struct 
{
  INDEX block_indx;//! block_indx block index of this part of matrix
	INDEX elements_num;//! elements_num number of elements
	INDEX rows_num;//! rows_num number of rows
  INDEX cols_num;//! cols_num number of columns
  INDEX block_rows_num;//! rows_num number of rows in block 
  INDEX block_cols_num;//! cols_num number of columns in block 
	INDEX* rows;//! row index array of matrix
	INDEX* cols;//! col index array of matrix
	REAL* values;//! data array of matrix
} cel_sparse_matrix_coo;

//! Initialize new sparse matricies in coordinate format (COO)
//! @param (in) matrix_num number of matrices to initialiez
//! @param (in) output_on 1 - print error to sdterr; 0 - no output
//! @return  pointer to 1 dimensional array with \struct cel_sparse_matrix_coo with nill initialized fields OR nil if no memory
cel_sparse_matrix_coo* cel_sparse_matrix_coo_new(INDEX matrix_num, int output_on);

//! Allocate memory for internal fields in \struct sparse matrix in coordinates format (COO)
//! @param (in) elements_num minimal number of non zero elements in matrix
//! @param (out)   pointer to \struct cel_sparse_matrix_coo
//! @param (in) output_on 1 - print error to sdterr; 0 - no output
//! @return  0 - no error ; -1 - not enought memory
int cel_sparse_matrix_coo_allocate(INDEX elements_num, cel_sparse_matrix_coo* sparse_matrix_coo, int output_on);

//! Free sparse matricies in coordinate format (COO) - only internal structure are realised
//! @param (in) pointer to 1 dimensional array with \struct cel_sparse_matrix_coo with nill initialized fields OR nil if no memory
//! @param (in) matrix_num number of matrices to initialiez
void cel_sparse_matrix_coo_free(cel_sparse_matrix_coo* sparse_matrix_coo, INDEX matrix_num);

//! Output of the sparse matrix in in coordinate format (COO) to stderr
//! @param (in) matrix_coo pointer to sparse matrix in coo format
//! @param (in) output_on 1 - print error to sdterr; 0 - no output
void cel_sparse_matrix_coo_print(cel_sparse_matrix_coo* sparse_matrix_coo, int output_on);

#endif
