/*     
 * File:   cel_vector.h
 * Project CRESTA (see details on https://cresta-project.eu) 
 * CRESTA Exascale library
 * Send you email about the bugs to
 * The High Performance Computing Center Stuttgart (HLRS) of the University of Stuttgart
 * Dmitry Khabi,  (email: khabi@hlrs.de)
 * Created on March 25, 2013
*/
#ifndef CEL_VECTOR_H
#define CEL_VECTOR_H
#include "cel_types.h"
//! \struct vector
typedef struct
{
	INDEX length;//! elements_num number of elements
	REAL* values;//! values value array of the vector
} cel_vector;

//! Initialize new vectors
//! @param (in) vector_num number of vectors
//! @param (in) output_on 1 - print error to sdterr; 0 - no output
//! @return  pointer to 1 dimensional array with \struct cel_sparse_matrix_coo with nill initialized fields OR nil if no memory
cel_vector* cel_vector_new(INDEX vector_num, int output_on);

//! Allocate memory for internal fields in \struct vector
//! @param (in) elements_num number of elements
//! @param (out)   pointer to \struct cel_vector
//! @param (in) output_on 1 - print error to sdterr; 0 - no output
//! @return  0 - no error ; -1 - not enought memory
int cel_vector_allocate(INDEX elements_num, cel_vector* vector, int output_on);

//! Free \struct cel_vector - only internal structure are realised
//! @param (in) pointer to 1 dimensional array with \struct cel_vector with nill initialized fields OR nil if no memory
//! @param (in) matrix_num number of matrices to initialiez
void cel_vector_free(cel_vector* vector, INDEX vector_num);

//! Output of the vector to stderr
//! @param (in) vector pointer to \struct vector
//! @param (in) output_on 1 - print error to sdterr; 0 - no output
void cel_vector_print(cel_vector* vector, int output_on);

#endif

