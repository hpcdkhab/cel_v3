/*     
 * File:   cel_sort.c -- see headers for details and function's comments
 * Project CRESTA (see details on https://cresta-project.eu) Exascale library
 * Send you email about the bugs to
 * The High Performance Computing Center Stuttgart (HLRS) of the University of Stuttgart
 * Dmitry Khabi,  (email: khabi@hlrs.de)
 * Created on March 25, 2013
*/
#include "cel_sort.h"
void cell_sort_asc_indices_with_values_bubble(INDEX* indices, REAL* values, INDEX length)
{
        int not_sorted = 1;
	INDEX jj=0;
	INDEX index_tmp;
	REAL real_tmp;
	INDEX ii;
			 
	while (not_sorted==1) 
	{
		not_sorted = 0;
		jj++;
		for (ii = 0; ii < length - jj; ii++) 
		{
			if (indices[ii] > indices[ii + 1]) 
			{
				index_tmp = indices[ii];
				indices[ii] = indices[ii + 1];
				indices[ii + 1] = index_tmp;
				real_tmp = values[ii];
				values[ii] = values[ii + 1];
				values[ii + 1] = real_tmp;
				not_sorted = 1;
			}
		}
	}
}

