/*     
 * File:   cel_sparse_matrix_coo.c -- see headers for details and function's comments
 * Project CRESTA (see details on https://cresta-project.eu)
 * CRESTA Exascale library
 * Send you email about the bugs to
 * The High Performance Computing Center Stuttgart (HLRS) of the University of Stuttgart
 * Dmitry Khabi,  (email: khabi@hlrs.de)
 * Created on March 25, 2013
*/

#include "cel_sparse_matrix_csr.h"
#include <stdlib.h>
#include <stdio.h>
cel_sparse_matrix_csr* cel_sparse_matrix_csr_new(INDEX matrix_num, int output_on)
{
        INDEX ii;
        void* temp_alloc;
	int err_code;
	cel_sparse_matrix_csr* sparse_matrix_csr;
        err_code=0;
        sparse_matrix_csr = Nullptr(cel_sparse_matrix_csr);     
        err_code += posix_memalign(&temp_alloc, (size_t)POSIX_MEM_ALIGN, sizeof(cel_sparse_matrix_csr)*(matrix_num));
        if(err_code > 0)
		{
			err_code = 1;
			if(output_on>0)
				fprintf(stderr,"ERROR in cel_sparse_matrix_csr_new -> not enought memory\n");
			return Nullptr(cel_sparse_matrix_csr);
		}
        sparse_matrix_csr =  (cel_sparse_matrix_csr*)temp_alloc;
        for(ii=0; ii<matrix_num;ii++)
        {
	        sparse_matrix_csr[ii].xadj_length = 0;
	        sparse_matrix_csr[ii].vtxdist = Nullptr(INDEX);
	        sparse_matrix_csr[ii].values = Nullptr(REAL);
	        sparse_matrix_csr[ii].diagonal = Nullptr(REAL);
	        sparse_matrix_csr[ii].rows_diag_ptrs = Nullptr(INDEX);
	        sparse_matrix_csr[ii].cols_diag_indices = Nullptr(INDEX);
			sparse_matrix_csr[ii].block_indx = INDEX_MAX;
			sparse_matrix_csr[ii].is_diagonal_separat= 0;
        }
        return sparse_matrix_csr;
}

int cel_sparse_matrix_csr_allocate(INDEX elements_num, INDEX rows_num, BOOLEAN is_diagonal_separat, cel_sparse_matrix_csr* sparse_matrix_csr, int output_on)
{
        INDEX length_el;
        int err_code;
        void* temp_alloc;

        err_code=0;
        sparse_matrix_csr->is_diagonal_separat=is_diagonal_separat;
	if(sparse_matrix_csr->is_diagonal_separat==1)
		length_el = elements_num-rows_num;
	else
		length_el = elements_num;
	err_code += posix_memalign(&temp_alloc, (size_t)32, sizeof(INDEX)*(rows_num+1));
	sparse_matrix_csr->xadj = temp_alloc;
	err_code += posix_memalign(&temp_alloc, (size_t)32, sizeof(INDEX)*length_el);
	sparse_matrix_csr->adjncy = (INDEX*)temp_alloc;
	err_code += posix_memalign(&temp_alloc, (size_t)32, sizeof(REAL)*length_el);
	sparse_matrix_csr->values = (REAL*)temp_alloc;
	if(sparse_matrix_csr->is_diagonal_separat==1)
		err_code += posix_memalign(&temp_alloc, (size_t)32, sizeof(REAL)*rows_num);
	else
		temp_alloc = (void*)0;
	sparse_matrix_csr->diagonal = (REAL*)temp_alloc;
	if(sparse_matrix_csr->is_diagonal_separat==1)
		err_code += posix_memalign(&temp_alloc, (size_t)32, sizeof(INDEX)*rows_num);
	else
		temp_alloc = (void*)0;
	sparse_matrix_csr->cols_diag_indices = (INDEX*)temp_alloc;
	if(sparse_matrix_csr->is_diagonal_separat==1)
		err_code += posix_memalign(&temp_alloc, (size_t)32, sizeof(INDEX)*(rows_num+1));
	else
		temp_alloc = (void*)0;
	sparse_matrix_csr->rows_diag_ptrs = (INDEX*)temp_alloc;
	
	if(err_code > 0)
	{
		err_code = 1;
		if(output_on>0)
			fprintf(stderr,"ERROR in cel_sparse_matrix_csr_allocate -> not enought memory\n");
		return err_code;
	}
        return err_code;
}

void cel_sparse_matrix_csr_free(cel_sparse_matrix_csr* sparse_matrix_csr, INDEX matrix_num)
{
        INDEX ii;
        for(ii=0; ii<matrix_num;ii++)
        {
	        free(sparse_matrix_csr[ii].vtxdist);
	        free(sparse_matrix_csr[ii].values);
	        free(sparse_matrix_csr[ii].diagonal);
	        free(sparse_matrix_csr[ii].rows_diag_ptrs);
	        free(sparse_matrix_csr[ii].cols_diag_indices);
	        sparse_matrix_csr[ii].xadj_length = 0;
	        sparse_matrix_csr[ii].vtxdist = Nullptr(INDEX);
	        sparse_matrix_csr[ii].values = Nullptr(REAL);
	        sparse_matrix_csr[ii].diagonal = Nullptr(REAL);
	        sparse_matrix_csr[ii].rows_diag_ptrs = Nullptr(INDEX);
	        sparse_matrix_csr[ii].cols_diag_indices = Nullptr(INDEX);
			sparse_matrix_csr[ii].block_indx = INDEX_MAX;
			sparse_matrix_csr[ii].is_diagonal_separat= 0;
        }
}

void cel_sparse_matrix_csr_print(cel_sparse_matrix_csr* sparse_matrix_csr, int output_on)
{
	INDEX ii,jj;
	if(output_on>0)
	{
		fprintf(stderr,"SPARSE MATRIX IN CSR FORMAT: ( 0 - "INDEX_FORMAT" ) BEGIN:", sparse_matrix_csr->xadj_length-1);
		
		for(ii=0; ii<sparse_matrix_csr->xadj_length-1; ii++)				
		{
			fprintf(stderr,"\nrow: "INDEX_FORMAT"("INDEX_FORMAT"-"INDEX_FORMAT");\n",ii,sparse_matrix_csr->xadj[ii],sparse_matrix_csr->xadj[ii+1]);
			for(jj=sparse_matrix_csr->xadj[ii]; jj<sparse_matrix_csr->xadj[ii+1]; jj++)				
				fprintf(stderr,INDEX_FORMAT";",sparse_matrix_csr->adjncy[jj]);
			fprintf(stderr,"\n");
			for(jj=sparse_matrix_csr->xadj[ii]; jj<sparse_matrix_csr->xadj[ii+1]; jj++)				
				fprintf(stderr,REAL_FORMAT";",sparse_matrix_csr->values[jj]);	
		}
		if(sparse_matrix_csr->is_diagonal_separat==1)
		{
			fprintf(stderr,"DIAGONAL:\n");
			for(ii=0; ii<sparse_matrix_csr->xadj_length-1; ii++)
			{
				fprintf(stderr,REAL_FORMAT";",sparse_matrix_csr->diagonal[ii]);	
			}
		}
		fprintf(stderr,"\nEND PARSE MATRIX IN CSR FORMAT\n");
	}
}


int cel_sparse_matrix_scr_from_sparse_matrix_coo(cel_sparse_matrix_coo* sparse_matrix_coo, cel_sparse_matrix_csr* sparse_matrix_csr, int output_on)
{
	int err_code;
	INDEX ii;
	INDEX end, length_rows;
	INDEX row, col, current_row;
	INDEX xadj_index, adjncy_index, row_index;
	REAL value;
	err_code=0;
	end = sparse_matrix_coo->elements_num;
	length_rows = sparse_matrix_coo->rows_num;
        err_code = cel_sparse_matrix_csr_allocate(end, length_rows, sparse_matrix_csr->is_diagonal_separat,sparse_matrix_csr,output_on);
	
	if(err_code > 0)
	{
		err_code = 1;
		if(output_on>0)
			fprintf(stderr,"ERROR in cel_sparse_matrix_coo_to_sparse_matrix_csr -> not enought memory\n");
		return err_code;
	}
	
	sparse_matrix_csr->block_indx = sparse_matrix_coo->block_indx;
	current_row = INDEX_MAX;
	xadj_index = 0;
	adjncy_index = 0;
	row_index = 0;
	if(sparse_matrix_csr->is_diagonal_separat==1)
	{
		for(ii=0; ii < end; ii++)
		{
			value = sparse_matrix_coo->values[ii];
			if(value>EPSILON || value < -EPSILON)
			{
				row = sparse_matrix_coo->rows[ii];
				col = sparse_matrix_coo->cols[ii];
				if(row!=current_row)
				{
					sparse_matrix_csr->xadj[xadj_index++]=adjncy_index;
					current_row = row;
				}
				if(row!=col)
				{
					sparse_matrix_csr->adjncy[adjncy_index]=col;
					sparse_matrix_csr->values[adjncy_index++]=value;
				}else
				{
					sparse_matrix_csr->diagonal[row_index++]=value;
				}
			}
		}
		for(ii=0; ii < length_rows; ii++)
		{
			sparse_matrix_csr->cols_diag_indices[ii]=ii;
			sparse_matrix_csr->rows_diag_ptrs[ii]=ii;
		}
		sparse_matrix_csr->rows_diag_ptrs[length_rows]=length_rows;
	
	}else
	{
		for(ii=0; ii < end; ii++)
		{
			value = sparse_matrix_coo->values[ii];
			if(value>EPSILON || value < -EPSILON)
			{
				row = sparse_matrix_coo->rows[ii];
				col = sparse_matrix_coo->cols[ii];
				if(row!=current_row)
				{
					sparse_matrix_csr->xadj[xadj_index++]=adjncy_index;
					current_row = row;
				}
				sparse_matrix_csr->adjncy[adjncy_index]=col;
				sparse_matrix_csr->values[adjncy_index++]=value;
			}
		}
	}
	sparse_matrix_csr->xadj[xadj_index]=adjncy_index;	
	sparse_matrix_csr->xadj_length=xadj_index+1;

	return err_code;
}

int cel_sparse_matrix_csr_init_vtxdist(cel_domain_3D* domain_3D, INDEX processes_num, cel_sparse_matrix_csr* sparse_matrix_csr, int output_on)
{
	int err_code;
	INDEX ii, rows_num;
	void* temp_alloc;
	INDEX* vtxdist;
	
	err_code=0;
	
	err_code += posix_memalign(&temp_alloc, (size_t)32, sizeof(INDEX)*(processes_num+1));
	sparse_matrix_csr->vtxdist = (INDEX*)temp_alloc;
	if(err_code > 0)
	{
		err_code = 1;
		if(output_on>0)
			fprintf(stderr,"ERROR in cel_sparse_matrix_csr_init_vtxdist -> not enought memory\n");
		return err_code;
	}
	vtxdist = sparse_matrix_csr->vtxdist;
	rows_num = sparse_matrix_csr->xadj_length-1;
	vtxdist[0]=0;
	for(ii=1; ii<processes_num+1; ii++)
	{
		vtxdist[ii]=vtxdist[ii-1]+rows_num;
	}
	
	return err_code;
}

void cel_sparse_matrix_csr_info(cel_sparse_matrix_csr* sparse_matrix_csr, INDEX *rows_number, INDEX *elements_number)
{
	*rows_number = sparse_matrix_csr->xadj_length - 1;
	*elements_number = sparse_matrix_csr->xadj[*rows_number] - 1;
	if(sparse_matrix_csr->is_diagonal_separat==1)
		*elements_number += *rows_number;
	
}
