/*     
 * File:   cel_sparse_matrix_coo.c -- see headers for details and function's comments
 * Project CRESTA (see details on https://cresta-project.eu)
 * CRESTA Exascale library
 * Send you email about the bugs to
 * The High Performance Computing Center Stuttgart (HLRS) of the University of Stuttgart
 * Dmitry Khabi,  (email: khabi@hlrs.de)
 * Created on March 25, 2013
*/

#include "cel_sparse_matrix_coo.h"
#include <stdlib.h>
#include <stdio.h>
cel_sparse_matrix_coo* cel_sparse_matrix_coo_new(INDEX matrix_num, int output_on)
{
        INDEX ii;
        void* temp_alloc;
	int err_code;
	cel_sparse_matrix_coo* sparse_matrix_coo;
        err_code=0;
        sparse_matrix_coo = Nullptr(cel_sparse_matrix_coo);     
        err_code += posix_memalign(&temp_alloc, (size_t)POSIX_MEM_ALIGN, sizeof(cel_sparse_matrix_coo)*(matrix_num));
        if(err_code > 0)
	{
		err_code = 1;
		if(output_on>1)
			fprintf(stderr,"ERROR in cel_new_sparse_matrix_coo -> not enought memory\n");
		return Nullptr(cel_sparse_matrix_coo);
	}
        sparse_matrix_coo =  (cel_sparse_matrix_coo*)temp_alloc;
        for(ii=0; ii<matrix_num;ii++)
        {
            sparse_matrix_coo[ii].block_indx = INDEX_MAX;
	        sparse_matrix_coo[ii].elements_num = 0;
	        sparse_matrix_coo[ii].rows_num = 0;
	        sparse_matrix_coo[ii].rows = Nullptr(INDEX);
	        sparse_matrix_coo[ii].cols = Nullptr(INDEX);
	        sparse_matrix_coo[ii].values = Nullptr(REAL);
        }
        return sparse_matrix_coo;
}

void cel_sparse_matrix_coo_free(cel_sparse_matrix_coo* sparse_matrix_coo, INDEX matrix_num)
{
        INDEX ii;
        for(ii=0; ii<matrix_num;ii++)
        {
	        free(sparse_matrix_coo[ii].rows);
	        free(sparse_matrix_coo[ii].cols);
	        free(sparse_matrix_coo[ii].values);
	        sparse_matrix_coo[ii].rows = Nullptr(INDEX);
	        sparse_matrix_coo[ii].cols = Nullptr(INDEX);
	        sparse_matrix_coo[ii].values = Nullptr(REAL);        
	    }
}

int cel_sparse_matrix_coo_allocate(INDEX elements_num, cel_sparse_matrix_coo* sparse_matrix_coo, int output_on)
{
        int err_code;
        void* temp_alloc;
        err_code=0;

	err_code += posix_memalign(&temp_alloc, (size_t)32, sizeof(INDEX)*elements_num);
	sparse_matrix_coo->rows =(INDEX*) temp_alloc;

	err_code += posix_memalign(&temp_alloc, (size_t)32, sizeof(INDEX)*elements_num);
	sparse_matrix_coo->cols =(INDEX*) temp_alloc;
	
	err_code += posix_memalign(&temp_alloc, (size_t)32, sizeof(REAL)*elements_num);
	sparse_matrix_coo->values =(REAL*) temp_alloc;
	if(err_code > 0)
	{
		err_code = 1;
		if(output_on>0)
			fprintf(stderr,"ERROR in cel_sparse_matrix_coo_allocate -> not enought memory\n");
		return err_code;
	}
        return err_code;
}

void cel_sparse_matrix_coo_print(cel_sparse_matrix_coo* sparse_matrix_coo, int output_on)
{
	INDEX ii;
	INDEX curr_row;
	if(output_on>0)
	{
		fprintf(stderr,"SPARSE MATRIX IN COO FORMAT (ELEMENTS: 0 - "INDEX_FORMAT" ) BEGIN:", sparse_matrix_coo->elements_num-1);
		curr_row = INDEX_MAX;
		for(ii=0; ii<sparse_matrix_coo->elements_num; ii++)				
		{
			if(curr_row != sparse_matrix_coo->rows[ii])
			{
				curr_row = sparse_matrix_coo->rows[ii];
				fprintf(stderr,"\nrow: "INDEX_FORMAT";\n",curr_row);
			}
			fprintf(stderr,INDEX_FORMAT":"REAL_FORMAT"; ",sparse_matrix_coo->cols[ii],sparse_matrix_coo->values[ii]);
		}
		fprintf(stderr,"\nEND SPARSE MATRIX IN COO FORMAT\n");
	}
}


