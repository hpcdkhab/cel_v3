/*     
 * File:   cel_sort.h
 * Project CRESTA (see details on https://cresta-project.eu) Exascale library
 * Send you email about the bugs to
 * The High Performance Computing Center Stuttgart (HLRS) of the University of Stuttgart
 * Dmitry Khabi,  (email: khabi@hlrs.de)
 * Created on March 25, 2013
*/
#ifndef CEL_SORT_H
#define CEL_SORT_H
#include "cel_types.h"
//! Sort two arrays indices (type INDEX*) and values (type REAL*) accroding to indices in ascending order
//! @param (out) indices array of type INDEX*
//! @param (out) values  array of type INDEX*
//! @param (out) length length of arrays
void cell_sort_asc_indices_with_values_bubble(INDEX* indices, REAL* values, INDEX length);

#endif

