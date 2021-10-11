/*     
 * File:   cel_discretisation_poisson_3d.c -- see headers for details and function's comments
 * Project CRESTA (see details on https://www.cresta-project.eu)
 * CRESTA Exascale library
 * Send you email about the bugs to
 * The High Performance Computing Center Stuttgart (HLRS) of the University of Stuttgart
 * Dmitry Khabi,  (email: khabi@hlrs.de)
 * Created on March 25, 2013
*/

#include "cel_discretisation_poisson_3d.h"
#include "cel_sort.h"
#include "cel_domain_3d.h"
#include <stdio.h>
int cel_discretisation_poisson_3D(  cel_discretisation_task* desc_task,   cel_domain_3D* domain_3D, INDEX block_indx,   cel_sparse_matrix_coo* sparse_matrix_coo, int output_on)
{
	int err_code;
	INDEX* rows;
	INDEX* cols;
	REAL* values;
	
	err_code = 0;
	INDEX ii,jj,kk;
	INDEX elements_num, curr_indx, length, rows_num;
	INDEX xx_start,yy_start,zz_start,xx_end,yy_end,zz_end;
        elements_num = 0;
	cel_domain_3D_get_block_boundary(domain_3D,block_indx,&xx_start,&yy_start,&zz_start,&xx_end,&yy_end,&zz_end);
	if(output_on>10)
	{
		fprintf(stderr,"xx_start: "INDEX_FORMAT"; xx_end: "INDEX_FORMAT"\n",xx_start,xx_end);
		fprintf(stderr,"yy_start: "INDEX_FORMAT"; yy_end: "INDEX_FORMAT"\n",yy_start,yy_end);
		fprintf(stderr,"zz_start: "INDEX_FORMAT"; zz_end: "INDEX_FORMAT"\n",zz_start,zz_end);
	}

	length = (zz_end-zz_start)*(yy_end-yy_start)*(xx_end-xx_start);
	switch (desc_task->stencil)
	{
			case STENCIL_7:
				   elements_num = length*7;
				   break;
			case STENCIL_19:
				   elements_num = length*19;
				   break;
			case STENCIL_27:
				   elements_num = length*27;
				   break;
	}

	err_code = cel_sparse_matrix_coo_allocate(elements_num, sparse_matrix_coo,output_on);	
	if(err_code > 0 )
	{
		err_code = 1;
		if(output_on>0)
			fprintf(stderr,"ERROR in cel_discretisation_poisson_3D -> not enought memory\n");
		return err_code;
	}
        rows = sparse_matrix_coo->rows;
        cols = sparse_matrix_coo->cols;
        values = sparse_matrix_coo->values;
	curr_indx = 0;
	rows_num = 0;
	domain_3D->block_indx=block_indx;
	sparse_matrix_coo->block_indx=block_indx;
	switch (desc_task->stencil)
	{
			case STENCIL_7:
			{
				   for(kk=zz_start; kk<zz_end; kk++)
						for(jj=yy_start; jj<yy_end; jj++)
							for(ii=xx_start; ii<xx_end; ii++)
								_poisson_3D_node_discr_7_stencil(desc_task ,domain_3D, ii, jj, kk, rows, cols, values, &curr_indx, &rows_num);
				   break;
			}
			case STENCIL_19:
			{
				   for(kk=zz_start; kk<zz_end; kk++)
						for(jj=yy_start; jj<yy_end; jj++)
							for(ii=xx_start; ii<xx_end; ii++)
								_poisson_3D_node_discr_19_stencil(desc_task ,domain_3D, ii, jj, kk, rows, cols, values, &curr_indx, &rows_num);
				   break;
			}
			case STENCIL_27:
			{
				   for(kk=zz_start; kk<zz_end; kk++)
						for(jj=yy_start; jj<yy_end; jj++)
							for(ii=xx_start; ii<xx_end; ii++)
								_poisson_3D_node_discr_27_stencil(desc_task ,domain_3D, ii, jj, kk, rows, cols, values, &curr_indx, &rows_num);
				   break;
			}
	}

	sparse_matrix_coo->elements_num = curr_indx;
	sparse_matrix_coo->rows_num = rows_num;
	sparse_matrix_coo->cols_num = rows_num;

	return err_code;
}


void _poisson_3D_node_discr_7_stencil(  cel_discretisation_task* desc_task,   cel_domain_3D* domain_3D, INDEX ii, INDEX jj, INDEX kk, INDEX* rows, INDEX* cols, REAL* values, INDEX *current_indx, INDEX* rows_num)
{
	INDEX curr_row_index, cell_indx_global, curr_indx, block_idx_global, cell_indx_local;
	REAL diag_value;
	REAL hX = 1.0/((double)(domain_3D->nx-1));
	REAL hY = 1.0/((double)(domain_3D->ny-1));
	REAL hZ = 1.0/((double)(domain_3D->nz-1));
	REAL hx2 = 1.0/(hX*hX), hy2 = 1.0/(hY*hY), hz2 = 1.0/(hZ*hZ);
	
	INDEX length;
	
	curr_indx = *current_indx;
	curr_row_index = cel_domain_3D_get_index(domain_3D, ii, jj, kk, &block_idx_global, &cell_indx_local);
	if(curr_row_index<INDEX_MAX)
	{
		diag_value = 0.0;
		//X
		cell_indx_global = cel_domain_3D_get_index(domain_3D, ii-1, jj, kk, &block_idx_global, &cell_indx_local);
		if(cell_indx_global<INDEX_MAX)
		{
			rows[curr_indx]=curr_row_index;
			cols[curr_indx]=cell_indx_global;
			values[curr_indx++]= -hx2;
			diag_value += hx2;
		}
		cell_indx_global = cel_domain_3D_get_index(domain_3D, ii+1, jj, kk, &block_idx_global, &cell_indx_local);
		if(cell_indx_global<INDEX_MAX)
		{
			rows[curr_indx]=curr_row_index;
			cols[curr_indx]=cell_indx_global;
			values[curr_indx++]= -hx2;
			diag_value += hx2;
		}
		//Y
		cell_indx_global = cel_domain_3D_get_index(domain_3D, ii, jj-1, kk, &block_idx_global, &cell_indx_local);
		if(cell_indx_global<INDEX_MAX)
		{
			rows[curr_indx]=curr_row_index;
			cols[curr_indx]=cell_indx_global;
			values[curr_indx++]= -hy2;
			diag_value += hy2;
		}
		cell_indx_global = cel_domain_3D_get_index(domain_3D, ii, jj+1, kk, &block_idx_global, &cell_indx_local);
		if(cell_indx_global<INDEX_MAX)
		{
			rows[curr_indx]=curr_row_index;
			cols[curr_indx]=cell_indx_global;
			values[curr_indx++]= -hy2;
			diag_value += hy2;
		}
		//Z
		cell_indx_global = cel_domain_3D_get_index(domain_3D, ii, jj, kk-1, &block_idx_global, &cell_indx_local);
		if(cell_indx_global<INDEX_MAX)
		{
			rows[curr_indx]=curr_row_index;
			cols[curr_indx]=cell_indx_global;
			values[curr_indx++]= -hz2;
			diag_value += hz2;
		}
		cell_indx_global = cel_domain_3D_get_index(domain_3D, ii, jj, kk+1, &block_idx_global, &cell_indx_local);
		if(cell_indx_global<INDEX_MAX)
		{
			rows[curr_indx]=curr_row_index;
			cols[curr_indx]=cell_indx_global;
			values[curr_indx++]= -hz2;
			diag_value += hz2;
		}

		

		rows[curr_indx]=curr_row_index;
		cols[curr_indx]=curr_row_index;
		values[curr_indx++]=diag_value;
		
		length = curr_indx - *current_indx;
		cell_sort_asc_indices_with_values_bubble(cols+(*current_indx),values+(*current_indx),length);		
		
		*current_indx = curr_indx;
		*rows_num += 1;
	}

}

/* According to
 * A Family of Large-Stencil Discrete Laplacian Approximations inThree Dimensions
 * R. C. O’REILLY AND J. M. BECK
 * ToDo: change the shema: works only for hx==hy==hz
 */
void _poisson_3D_node_discr_19_stencil(  cel_discretisation_task* desc_task,   cel_domain_3D* domain_3D, INDEX ii, INDEX jj, INDEX kk, INDEX* rows, INDEX* cols, REAL* values, INDEX *current_indx, INDEX* rows_num)
{
	INDEX curr_row_index, cell_indx_global, curr_indx, block_idx_global, cell_indx_local;
	REAL diag_value;
	REAL hX = 1.0/((double)(domain_3D->nx-1));
	REAL hY = 1.0/((double)(domain_3D->ny-1));
	REAL hZ = 1.0/((double)(domain_3D->nz-1));
	REAL hx2 = 1.0/(6.0*hX*hX), hy2 = 1.0/(6.0*hY*hY), hz2 = 1.0/(6.0*hZ*hZ);
	REAL hxy2=hx2,hxz2=hx2,hyz2=hy2;
	
	INDEX length;
	
	curr_indx = *current_indx;
	curr_row_index = cel_domain_3D_get_index(domain_3D, ii, jj, kk, &block_idx_global, &cell_indx_local);
	if(curr_row_index<INDEX_MAX)
	{
		diag_value = 0.0;
		//X
		cell_indx_global = cel_domain_3D_get_index(domain_3D, ii-1, jj, kk, &block_idx_global, &cell_indx_local);
		if(cell_indx_global<INDEX_MAX)
		{
			rows[curr_indx]=curr_row_index;
			cols[curr_indx]=cell_indx_global;
			values[curr_indx++]= -2.0*hx2;
			diag_value += hx2;
		}
		cell_indx_global = cel_domain_3D_get_index(domain_3D, ii+1, jj, kk, &block_idx_global, &cell_indx_local);
		if(cell_indx_global<INDEX_MAX)
		{
			rows[curr_indx]=curr_row_index;
			cols[curr_indx]=cell_indx_global;
			values[curr_indx++]= -2.0*hx2;
			diag_value += hx2;
		}
		//Y
		cell_indx_global = cel_domain_3D_get_index(domain_3D, ii, jj-1, kk, &block_idx_global, &cell_indx_local);
		if(cell_indx_global<INDEX_MAX)
		{
			rows[curr_indx]=curr_row_index;
			cols[curr_indx]=cell_indx_global;
			values[curr_indx++]= -2.0*hy2;
			diag_value += hy2;
		}
		cell_indx_global = cel_domain_3D_get_index(domain_3D, ii, jj+1, kk, &block_idx_global, &cell_indx_local);
		if(cell_indx_global<INDEX_MAX)
		{
			rows[curr_indx]=curr_row_index;
			cols[curr_indx]=cell_indx_global;
			values[curr_indx++]= -2.0*hy2;
			diag_value += hy2;
		}
		//Z
		cell_indx_global = cel_domain_3D_get_index(domain_3D, ii, jj, kk-1, &block_idx_global, &cell_indx_local);
		if(cell_indx_global<INDEX_MAX)
		{
			rows[curr_indx]=curr_row_index;
			cols[curr_indx]=cell_indx_global;
			values[curr_indx++]= -2.0*hz2;
			diag_value += hz2;
		}
		cell_indx_global = cel_domain_3D_get_index(domain_3D, ii, jj, kk+1, &block_idx_global, &cell_indx_local);
		if(cell_indx_global<INDEX_MAX)
		{
			rows[curr_indx]=curr_row_index;
			cols[curr_indx]=cell_indx_global;
			values[curr_indx++]= -2.0*hz2;
			diag_value += hz2;
		}
		//XY
		cell_indx_global = cel_domain_3D_get_index(domain_3D, ii-1, jj+1, kk, &block_idx_global, &cell_indx_local);
		if(cell_indx_global<INDEX_MAX)
		{
			rows[curr_indx]=curr_row_index;
			cols[curr_indx]=cell_indx_global;
			values[curr_indx++]= -hxy2;
			diag_value += hxy2;
		}
		cell_indx_global = cel_domain_3D_get_index(domain_3D, ii+1, jj-1, kk, &block_idx_global, &cell_indx_local);
		if(cell_indx_global<INDEX_MAX)
		{
			rows[curr_indx]=curr_row_index;
			cols[curr_indx]=cell_indx_global;
			values[curr_indx++]= -hxy2;
			diag_value += hxy2;
		}
		cell_indx_global = cel_domain_3D_get_index(domain_3D, ii-1, jj-1, kk, &block_idx_global, &cell_indx_local);
		if(cell_indx_global<INDEX_MAX)
		{
			rows[curr_indx]=curr_row_index;
			cols[curr_indx]=cell_indx_global;
			values[curr_indx++]= -hxy2;
			diag_value += hxy2;
		}
		cell_indx_global = cel_domain_3D_get_index(domain_3D, ii+1, jj+1, kk, &block_idx_global, &cell_indx_local);
		if(cell_indx_global<INDEX_MAX)
		{
			rows[curr_indx]=curr_row_index;
			cols[curr_indx]=cell_indx_global;
			values[curr_indx++]= -hxy2;
			diag_value += hxy2;
		}
		//XZ
		cell_indx_global = cel_domain_3D_get_index(domain_3D, ii-1, jj, kk+1, &block_idx_global, &cell_indx_local);
		if(cell_indx_global<INDEX_MAX)
		{
			rows[curr_indx]=curr_row_index;
			cols[curr_indx]=cell_indx_global;
			values[curr_indx++]= -hxz2;
			diag_value += hxz2;
		}
		cell_indx_global = cel_domain_3D_get_index(domain_3D, ii+1, jj, kk-1, &block_idx_global, &cell_indx_local);
		if(cell_indx_global<INDEX_MAX)
		{
			rows[curr_indx]=curr_row_index;
			cols[curr_indx]=cell_indx_global;
			values[curr_indx++]= -hxz2;
			diag_value += hxz2;
		}
		cell_indx_global = cel_domain_3D_get_index(domain_3D, ii-1, jj, kk-1, &block_idx_global, &cell_indx_local);
		if(cell_indx_global<INDEX_MAX)
		{
			rows[curr_indx]=curr_row_index;
			cols[curr_indx]=cell_indx_global;
			values[curr_indx++]= -hxz2;
			diag_value += hxz2;
		}
		cell_indx_global = cel_domain_3D_get_index(domain_3D, ii+1, jj, kk+1, &block_idx_global, &cell_indx_local);
		if(cell_indx_global<INDEX_MAX)
		{
			rows[curr_indx]=curr_row_index;
			cols[curr_indx]=cell_indx_global;
			values[curr_indx++]= -hxz2;
			diag_value += hxz2;
		}
		//YZ
		cell_indx_global = cel_domain_3D_get_index(domain_3D, ii, jj-1, kk+1, &block_idx_global, &cell_indx_local);
		if(cell_indx_global<INDEX_MAX)
		{
			rows[curr_indx]=curr_row_index;
			cols[curr_indx]=cell_indx_global;
			values[curr_indx++]= -hyz2;
			diag_value += hyz2;
		}
		cell_indx_global = cel_domain_3D_get_index(domain_3D, ii, jj+1, kk-1, &block_idx_global, &cell_indx_local);
		if(cell_indx_global<INDEX_MAX)
		{
			rows[curr_indx]=curr_row_index;
			cols[curr_indx]=cell_indx_global;
			values[curr_indx++]= -hyz2;
			diag_value += hyz2;
		}
		cell_indx_global = cel_domain_3D_get_index(domain_3D, ii, jj-1, kk-1, &block_idx_global, &cell_indx_local);
		if(cell_indx_global<INDEX_MAX)
		{
			rows[curr_indx]=curr_row_index;
			cols[curr_indx]=cell_indx_global;
			values[curr_indx++]= -hyz2;
			diag_value += hyz2;
		}
		cell_indx_global = cel_domain_3D_get_index(domain_3D, ii, jj+1, kk+1, &block_idx_global, &cell_indx_local);
		if(cell_indx_global<INDEX_MAX)
		{
			rows[curr_indx]=curr_row_index;
			cols[curr_indx]=cell_indx_global;
			values[curr_indx++]= -hyz2;
			diag_value += hyz2;
		}
			
		

		rows[curr_indx]=curr_row_index;
		cols[curr_indx]=curr_row_index;
		values[curr_indx++]=24*hx2;//diag_value;
		
		length = curr_indx - *current_indx;
		cell_sort_asc_indices_with_values_bubble(cols+(*current_indx),values+(*current_indx),length);		
		
		*current_indx = curr_indx;
		*rows_num += 1;
	}

}
/* According to
 * A Family of Large-Stencil Discrete Laplacian Approximations inThree Dimensions
 * R. C. O’REILLY AND J. M. BECK
 * ToDo: change the shema: works only for hx==hy==hz
 */
void _poisson_3D_node_discr_27_stencil(  cel_discretisation_task* desc_task,   cel_domain_3D* domain_3D, INDEX ii, INDEX jj, INDEX kk, INDEX* rows, INDEX* cols, REAL* values, INDEX *current_indx, INDEX* rows_num)
{
	INDEX curr_row_index, cell_indx_global, curr_indx, block_idx_global, cell_indx_local;
	REAL diag_value;
	REAL hX = 1.0/((double)(domain_3D->nx-1));
	REAL hY = 1.0/((double)(domain_3D->ny-1));
	REAL hZ = 1.0/((double)(domain_3D->nz-1));
	REAL hx2 = 3.0/(13.0*hX*hX), hy2 = 3.0/(13.0*hY*hY), hz2 = 3.0/(13.0*hZ*hZ);
	REAL hxy2=(1.0/2.0)*hx2,hxz2=(1.0/2.0)*hy2,hyz2=(1.0/2.0)*hy2, hxyz2=(1.0/3.0)*hy2;
	
	INDEX length;
	
	curr_indx = *current_indx;
	curr_row_index = cel_domain_3D_get_index(domain_3D, ii, jj, kk, &block_idx_global, &cell_indx_local);
	if(curr_row_index<INDEX_MAX)
	{
		diag_value = 0.0;
		//X
		cell_indx_global = cel_domain_3D_get_index(domain_3D, ii-1, jj, kk, &block_idx_global, &cell_indx_local);
		if(cell_indx_global<INDEX_MAX)
		{
			rows[curr_indx]=curr_row_index;
			cols[curr_indx]=cell_indx_global;
			values[curr_indx++]= -hx2;
			diag_value += hx2;
		}
		cell_indx_global = cel_domain_3D_get_index(domain_3D, ii+1, jj, kk, &block_idx_global, &cell_indx_local);
		if(cell_indx_global<INDEX_MAX)
		{
			rows[curr_indx]=curr_row_index;
			cols[curr_indx]=cell_indx_global;
			values[curr_indx++]= -hx2;
			diag_value += hx2;
		}
		//Y
		cell_indx_global = cel_domain_3D_get_index(domain_3D, ii, jj-1, kk, &block_idx_global, &cell_indx_local);
		if(cell_indx_global<INDEX_MAX)
		{
			rows[curr_indx]=curr_row_index;
			cols[curr_indx]=cell_indx_global;
			values[curr_indx++]= -hy2;
			diag_value += hy2;
		}
		cell_indx_global = cel_domain_3D_get_index(domain_3D, ii, jj+1, kk, &block_idx_global, &cell_indx_local);
		if(cell_indx_global<INDEX_MAX)
		{
			rows[curr_indx]=curr_row_index;
			cols[curr_indx]=cell_indx_global;
			values[curr_indx++]= -hy2;
			diag_value += hy2;
		}
		//Z
		cell_indx_global = cel_domain_3D_get_index(domain_3D, ii, jj, kk-1, &block_idx_global, &cell_indx_local);
		if(cell_indx_global<INDEX_MAX)
		{
			rows[curr_indx]=curr_row_index;
			cols[curr_indx]=cell_indx_global;
			values[curr_indx++]= -hz2;
			diag_value += hz2;
		}
		cell_indx_global = cel_domain_3D_get_index(domain_3D, ii, jj, kk+1, &block_idx_global, &cell_indx_local);
		if(cell_indx_global<INDEX_MAX)
		{
			rows[curr_indx]=curr_row_index;
			cols[curr_indx]=cell_indx_global;
			values[curr_indx++]= -hz2;
			diag_value += hz2;
		}

		//XY
		cell_indx_global = cel_domain_3D_get_index(domain_3D, ii-1, jj+1, kk, &block_idx_global, &cell_indx_local);
		if(cell_indx_global<INDEX_MAX)
		{
			rows[curr_indx]=curr_row_index;
			cols[curr_indx]=cell_indx_global;
			values[curr_indx++]= -hxy2;
			diag_value += hxy2;
		}
		cell_indx_global = cel_domain_3D_get_index(domain_3D, ii+1, jj-1, kk, &block_idx_global, &cell_indx_local);
		if(cell_indx_global<INDEX_MAX)
		{
			rows[curr_indx]=curr_row_index;
			cols[curr_indx]=cell_indx_global;
			values[curr_indx++]= -hxy2;
			diag_value += hxy2;
		}
		cell_indx_global = cel_domain_3D_get_index(domain_3D, ii-1, jj-1, kk, &block_idx_global, &cell_indx_local);
		if(cell_indx_global<INDEX_MAX)
		{
			rows[curr_indx]=curr_row_index;
			cols[curr_indx]=cell_indx_global;
			values[curr_indx++]= -hxy2;
			diag_value += hxy2;
		}
		cell_indx_global = cel_domain_3D_get_index(domain_3D, ii+1, jj+1, kk, &block_idx_global, &cell_indx_local);
		if(cell_indx_global<INDEX_MAX)
		{
			rows[curr_indx]=curr_row_index;
			cols[curr_indx]=cell_indx_global;
			values[curr_indx++]= -hxy2;
			diag_value += hxy2;
		}
		//XZ
		cell_indx_global = cel_domain_3D_get_index(domain_3D, ii-1, jj, kk+1, &block_idx_global, &cell_indx_local);
		if(cell_indx_global<INDEX_MAX)
		{
			rows[curr_indx]=curr_row_index;
			cols[curr_indx]=cell_indx_global;
			values[curr_indx++]= -hxz2;
			diag_value += hxz2;
		}
		cell_indx_global = cel_domain_3D_get_index(domain_3D, ii+1, jj, kk-1, &block_idx_global, &cell_indx_local);
		if(cell_indx_global<INDEX_MAX)
		{
			rows[curr_indx]=curr_row_index;
			cols[curr_indx]=cell_indx_global;
			values[curr_indx++]= -hxz2;
			diag_value += hxz2;
		}
		cell_indx_global = cel_domain_3D_get_index(domain_3D, ii-1, jj, kk-1, &block_idx_global, &cell_indx_local);
		if(cell_indx_global<INDEX_MAX)
		{
			rows[curr_indx]=curr_row_index;
			cols[curr_indx]=cell_indx_global;
			values[curr_indx++]= -hxz2;
			diag_value += hxz2;
		}
		cell_indx_global = cel_domain_3D_get_index(domain_3D, ii+1, jj, kk+1, &block_idx_global, &cell_indx_local);
		if(cell_indx_global<INDEX_MAX)
		{
			rows[curr_indx]=curr_row_index;
			cols[curr_indx]=cell_indx_global;
			values[curr_indx++]= -hxz2;
			diag_value += hxz2;
		}
		//YZ
		cell_indx_global = cel_domain_3D_get_index(domain_3D, ii, jj-1, kk+1, &block_idx_global, &cell_indx_local);
		if(cell_indx_global<INDEX_MAX)
		{
			rows[curr_indx]=curr_row_index;
			cols[curr_indx]=cell_indx_global;
			values[curr_indx++]= -hyz2;
			diag_value += hyz2;
		}
		cell_indx_global = cel_domain_3D_get_index(domain_3D, ii, jj+1, kk-1, &block_idx_global, &cell_indx_local);
		if(cell_indx_global<INDEX_MAX)
		{
			rows[curr_indx]=curr_row_index;
			cols[curr_indx]=cell_indx_global;
			values[curr_indx++]= -hyz2;
			diag_value += hyz2;
		}
		cell_indx_global = cel_domain_3D_get_index(domain_3D, ii, jj-1, kk-1, &block_idx_global, &cell_indx_local);
		if(cell_indx_global<INDEX_MAX)
		{
			rows[curr_indx]=curr_row_index;
			cols[curr_indx]=cell_indx_global;
			values[curr_indx++]= -hyz2;
			diag_value += hyz2;
		}
		cell_indx_global = cel_domain_3D_get_index(domain_3D, ii, jj+1, kk+1, &block_idx_global, &cell_indx_local);
		if(cell_indx_global<INDEX_MAX)
		{
			rows[curr_indx]=curr_row_index;
			cols[curr_indx]=cell_indx_global;
			values[curr_indx++]= -hyz2;
			diag_value += hyz2;
		}
		//XYZ
		cell_indx_global = cel_domain_3D_get_index(domain_3D, ii-1, jj-1, kk-1, &block_idx_global, &cell_indx_local);
		if(cell_indx_global<INDEX_MAX)
		{
			rows[curr_indx]=curr_row_index;
			cols[curr_indx]=cell_indx_global;
			values[curr_indx++]= -hxyz2;
			diag_value += hxyz2;
		}
		cell_indx_global = cel_domain_3D_get_index(domain_3D, ii+1, jj+1, kk+1, &block_idx_global, &cell_indx_local);
		if(cell_indx_global<INDEX_MAX)
		{
			rows[curr_indx]=curr_row_index;
			cols[curr_indx]=cell_indx_global;
			values[curr_indx++]= -hxyz2;
			diag_value += hxyz2;
		}
		cell_indx_global = cel_domain_3D_get_index(domain_3D, ii-1, jj+1, kk-1, &block_idx_global, &cell_indx_local);
		if(cell_indx_global<INDEX_MAX)
		{
			rows[curr_indx]=curr_row_index;
			cols[curr_indx]=cell_indx_global;
			values[curr_indx++]= -hxyz2;
			diag_value += hxyz2;
		}
		cell_indx_global = cel_domain_3D_get_index(domain_3D, ii+1, jj-1, kk+1, &block_idx_global, &cell_indx_local);
		if(cell_indx_global<INDEX_MAX)
		{
			rows[curr_indx]=curr_row_index;
			cols[curr_indx]=cell_indx_global;
			values[curr_indx++]= -hxyz2;
			diag_value += hxyz2;
		}
		cell_indx_global = cel_domain_3D_get_index(domain_3D, ii-1, jj+1, kk+1, &block_idx_global, &cell_indx_local);
		if(cell_indx_global<INDEX_MAX)
		{
			rows[curr_indx]=curr_row_index;
			cols[curr_indx]=cell_indx_global;
			values[curr_indx++]= -hxyz2;
			diag_value += hxyz2;
		}
		cell_indx_global = cel_domain_3D_get_index(domain_3D, ii+1, jj-1, kk-1, &block_idx_global, &cell_indx_local);
		if(cell_indx_global<INDEX_MAX)
		{
			rows[curr_indx]=curr_row_index;
			cols[curr_indx]=cell_indx_global;
			values[curr_indx++]= -hxyz2;
			diag_value += hxyz2;
		}
		cell_indx_global = cel_domain_3D_get_index(domain_3D, ii-1, jj-1, kk+1, &block_idx_global, &cell_indx_local);
		if(cell_indx_global<INDEX_MAX)
		{
			rows[curr_indx]=curr_row_index;
			cols[curr_indx]=cell_indx_global;
			values[curr_indx++]= -hxyz2;
			diag_value += hxyz2;
		}
		cell_indx_global = cel_domain_3D_get_index(domain_3D, ii+1, jj+1, kk-1, &block_idx_global, &cell_indx_local);
		if(cell_indx_global<INDEX_MAX)
		{
			rows[curr_indx]=curr_row_index;
			cols[curr_indx]=cell_indx_global;
			values[curr_indx++]= -hxyz2;
			diag_value += hxyz2;
		}
			
		rows[curr_indx]=curr_row_index;
		cols[curr_indx]=curr_row_index;
		values[curr_indx++]=(44.0/3.0)*hx2;
		
		length = curr_indx - *current_indx;
		cell_sort_asc_indices_with_values_bubble(cols+(*current_indx),values+(*current_indx),length);		
		
		*current_indx = curr_indx;
		*rows_num += 1;
	}

}

int cel_discretisation_poisson_3D_boundary_conditions(  cel_discretisation_task* desc_task,   cel_domain_3D* domain_3D, INDEX block_indx,   cel_vector* boundary_vector,   cel_vector* x_vector,   cel_vector* solution_vector, int output_on)
{
	REAL hX = 1.0/((REAL)(domain_3D->nx-1));
	REAL hY = 1.0/((REAL)(domain_3D->ny-1));
	REAL hZ = 1.0/((REAL)(domain_3D->nz-1));
	INDEX xx_start, xx_end, yy_start, yy_end, zz_start, zz_end;
	INDEX elements_num;
	int err_code;
	REAL* values;
	REAL node_x, node_y, node_z;
	INDEX ii,jj,kk;
	INDEX cell_indx_local,block_idx_global;
	
	err_code = 0;
	cel_domain_3D_get_block_boundary(domain_3D,block_indx,&xx_start,&yy_start,&zz_start,&xx_end,&yy_end,&zz_end);
	
	elements_num = (zz_end-zz_start)*(yy_end-yy_start)*(xx_end-xx_start);
	err_code += cel_vector_allocate(elements_num,boundary_vector,output_on); 
	err_code += cel_vector_allocate(elements_num,x_vector,output_on); 
	err_code += cel_vector_allocate(elements_num,solution_vector,output_on); 
	if(err_code > 0)
	{
		err_code = 1;
		if(output_on>0)
			fprintf(stderr,"ERROR in cel_discretisation_poisson_3D_boundary_conditions\n");
		return err_code;
	}
	values = boundary_vector->values;
	for(kk=zz_start; kk<zz_end; kk++)				
	{
		for(jj=yy_start; jj<yy_end; jj++)
		{
			for(ii=xx_start; ii<xx_end; ii++)
			{
        cel_domain_3D_get_index(domain_3D, ii, jj, kk, &block_idx_global, &cell_indx_local);
				node_x=ii*hX;
				node_y=jj*hY;
				node_z=kk*hZ;
				solution_vector->values[cell_indx_local]=_poisson_3D_solution(desc_task,node_x,node_y,node_z);
				//x_vector->values[cell_indx_local]=_poisson_3D_solution(desc_task,node_x,node_y,node_z);
				if(ii == 0 || ii == domain_3D->nx-1 || jj == 0 || jj == domain_3D->ny-1 || kk == 0 || kk == domain_3D->nz-1)
				{
					values[cell_indx_local] = _poisson_3D_rhs_boundary(desc_task,node_x,node_y,node_z);
				}
				else
				{
					values[cell_indx_local] = _poisson_3D_rhs_inner(desc_task,node_x,node_y,node_z);
				}
			}
		}
	}
	
	values = x_vector->values;
	for(ii=0; ii<elements_num; ii++)
		values[ii]=(REAL)0.0;
	
	return err_code;
}

int cel_discretisation_poisson_3D_eliminate_dirichlet_values(  cel_discretisation_task* desc_task,   cel_domain_3D* domain_3D,   cel_vector* boundary_vector,  cel_sparse_matrix_coo* sparse_matrix_coo, int output_on)
{
	INDEX end;
	INDEX ii;
	INDEX block_idx_global,cell_indx_local;
	INDEX row, col;
	INDEX xx_row,yy_row,zz_row;
	INDEX xx_col,yy_col,zz_col;
	REAL node_x, node_y, node_z;
	REAL hX = 1.0/((REAL)(domain_3D->nx-1));
	REAL hY = 1.0/((REAL)(domain_3D->ny-1));
	REAL hZ = 1.0/((REAL)(domain_3D->nz-1));
	int err_code;
	
	err_code=0;
	end = sparse_matrix_coo->elements_num;
	
	for(ii=0; ii < end; ii++)
	{
		row = sparse_matrix_coo->rows[ii];
		col = sparse_matrix_coo->cols[ii];
		cel_domain_3D_get_node_coordinate_indicies(domain_3D,row,&xx_row,&yy_row,&zz_row);
		if(xx_row == 0 || xx_row == domain_3D->nx-1 || yy_row == 0 || yy_row == domain_3D->ny-1 || zz_row == 0 || zz_row == domain_3D->nz-1)
		{
			if(row!=col)
				sparse_matrix_coo->values[ii]=0.0;
			else
				sparse_matrix_coo->values[ii]=1.0;
		}else
		{
			if(row!=col)
			{
				cel_domain_3D_get_node_coordinate_indicies(domain_3D,col,&xx_col,&yy_col,&zz_col);
				if(xx_col == 0 || xx_col == domain_3D->nx-1 || yy_col == 0 || yy_col == domain_3D->ny-1 || zz_col == 0 || zz_col == domain_3D->nz-1)
				{
					node_x=((REAL)xx_col)*hX;
					node_y=((REAL)yy_col)*hY;
					node_z=((REAL)zz_col)*hZ;
					cel_domain_3D_get_index(domain_3D, xx_row, yy_row, zz_row, &block_idx_global, &cell_indx_local);
			/*		fprintf(stderr,"______________cell_indx_local:"INDEX_FORMAT"\n", cell_indx_local);
					fprintf(stderr,"col: "INDEX_FORMAT"; xx_col="INDEX_FORMAT"; yy_col="INDEX_FORMAT"; zz_col="INDEX_FORMAT";\n",col,xx_col,yy_col,zz_col);
					fprintf(stderr,"col: "INDEX_FORMAT"; node_x="REAL_FORMAT"; node_y="REAL_FORMAT"; node_z="REAL_FORMAT";\n",col,node_x,node_y,node_z);
					fprintf(stderr,"OLD B(ROW)="REAL_FORMAT"; A(COL)="REAL_FORMAT"; B(COL)="REAL_FORMAT";\n",boundary_vector->values[cell_indx_local],sparse_matrix_coo->values[ii],_poisson_3D_rhs_boundary(desc_task,node_x,node_y,node_z));
			*/		
					boundary_vector->values[cell_indx_local]-=sparse_matrix_coo->values[ii]*_poisson_3D_rhs_boundary(desc_task,node_x,node_y,node_z);
					sparse_matrix_coo->values[ii]=0.0;
			/*		fprintf(stderr,"NEW B(ROW)="REAL_FORMAT";\n",boundary_vector->values[cell_indx_local]);
					fprintf(stderr,"______________\n");*/
				}
			}
		}
	}
	return err_code;
}

REAL cel_discretisation_poisson_3D_error(  cel_discretisation_task* desc_task,   cel_domain_3D* domain_3D,   cel_vector* x_vector)
{
	REAL hX = 1.0/((REAL)(domain_3D->nx-1));
	REAL hY = 1.0/((REAL)(domain_3D->ny-1));
	REAL hZ = 1.0/((REAL)(domain_3D->nz-1));
	REAL h3 = hX*hY*hZ;
	REAL* values;
	REAL err;
	INDEX ii, first_b_index, end, xx_node, yy_node, zz_node;
	REAL xx, yy, zz;
	REAL result = 0;
	first_b_index = domain_3D->block_indx*domain_3D->bxyz;
	values = x_vector->values;
	end = x_vector->length;
	for(ii=0; ii < end; ii++)
	{
		cel_domain_3D_get_node_coordinate_indicies(domain_3D,first_b_index+ii,&xx_node,&yy_node,&zz_node);
		xx = (REAL) xx_node*hX; yy = (REAL) yy_node*hY; zz = (REAL) zz_node*hZ;
		err = _poisson_3D_solution(desc_task,xx,yy,zz) - values[ii];
		result += err*err *h3;
	}
	
	return result;
}

REAL _poisson_3D_rhs_inner(  cel_discretisation_task*  desc_task, REAL xx, REAL yy, REAL zz)
{
	switch (desc_task->boundary_cond)
	{
			case X_2_PLUS_Y_2_PLUS_Z_2:
				return -6.0;
				break;
                        case USER_DEFINED:
				return -6.0;
				break;
	}
	return REAL_UNDEF;                
}

REAL _poisson_3D_rhs_boundary(  cel_discretisation_task*  desc_task, REAL xx, REAL yy, REAL zz)
{
	switch (desc_task->boundary_cond)
	{
			case X_2_PLUS_Y_2_PLUS_Z_2:
				return (xx*xx+yy*yy+zz*zz);
                                break;
                        case USER_DEFINED:
				return (xx*xx+yy*yy+zz*zz);
				break;
	}
	return REAL_UNDEF;  
}

REAL _poisson_3D_solution(  cel_discretisation_task*  desc_task, REAL xx, REAL yy, REAL zz)
{
	switch (desc_task->boundary_cond)
	{
			case X_2_PLUS_Y_2_PLUS_Z_2:
				return  (xx*xx+yy*yy+zz*zz);
				break;
                        case USER_DEFINED:
				return (xx*xx+yy*yy+zz*zz);
				break;
	}
	return REAL_UNDEF; 
}

