/*     
 * File:   cel_domain_3d.c -- see headers for details and function's comments
 * Project CRESTA (see details on https://www.cresta-project.eu)
 * CRESTA Exascale library
 * Send you email about the bugs to
 * The High Performance Computing Center Stuttgart (HLRS) of the University of Stuttgart
 * Dmitry Khabi,  (email: khabi@hlrs.de)
 * Created on March 25, 2013
*/

#include "cel_domain_3d.h"
#include <stdio.h>
int cel_domain_3D_init(INDEX nx, INDEX ny, INDEX nz, INDEX dx, INDEX dy, INDEX dz,   cel_domain_3D* domain_3D, int output_on)
{
	int err_code;
	
	err_code = 0;
	domain_3D->nx = nx;
	domain_3D->ny = ny;
	domain_3D->nz = nz;
	domain_3D->dx = dx;
	domain_3D->dy = dy;
	domain_3D->dz = dz;
	domain_3D->bx = nx/dx;
	domain_3D->by = ny/dy;
	domain_3D->bz = nz/dz;
	domain_3D->bxyz = (nx/dx)*(ny/dy)*(nz/dz);
	domain_3D->blocks = dx*dy*dz;
	if(domain_3D->blocks*domain_3D->bxyz != domain_3D->nx*domain_3D->ny*domain_3D->nz)
	{
		err_code = 1;
		if (output_on>0)
		{
			fprintf(stderr,"ERROR in disct_get_domain_descr_3D -> wrong decomposition's parameters:\n");
		}
	}
	if (output_on>0)
	{
		cel_domain_3D_print(domain_3D,1);
	}
	return err_code;
}

void cel_domain_3D_print(  cel_domain_3D* domain_3D, int output_on)
{
	if(output_on>10)
	{
		fprintf(stderr,"nx number of nodes in X-direction: "INDEX_FORMAT";\n",domain_3D->nx);
		fprintf(stderr,"ny number of nodes in Y-direction: "INDEX_FORMAT";\n",domain_3D->ny);
		fprintf(stderr,"nz number of nodes in Z-direction: "INDEX_FORMAT";\n",domain_3D->nz);
		fprintf(stderr,"dx number of sub-domains in X-direction: "INDEX_FORMAT";\n",domain_3D->dx);
		fprintf(stderr,"dy number of sub-domains in Y-direction: "INDEX_FORMAT";\n",domain_3D->dy);
		fprintf(stderr,"dz number of sub-domains in Z-direction: "INDEX_FORMAT";\n",domain_3D->dz);
		fprintf(stderr,"blocks number blocks (sub-domains) in 3D domain: "INDEX_FORMAT";\n",domain_3D->blocks);
		fprintf(stderr,"bx number nodes in X-direction in sub-domain: "INDEX_FORMAT";\n",domain_3D->bx);
		fprintf(stderr,"by number nodes in Y-direction in sub-domain: "INDEX_FORMAT";\n",domain_3D->by);
		fprintf(stderr,"bz number nodes in Z-direction in sub-domain: "INDEX_FORMAT";\n",domain_3D->bz);
		fprintf(stderr,"bxyz number nodes in sub-domain:"INDEX_FORMAT";\n",domain_3D->bxyz);
	}
}

INDEX cel_domain_3D_get_index(  cel_domain_3D* domain_3D, INDEX node_x_indx, INDEX node_y_indx, INDEX node_z_indx,INDEX* block_idx, INDEX* cell_indx_local)
{
	INDEX cell_indx;
	INDEX plane_x, plane_y, plane_z;
	INDEX cell_indx_global;
	
	if(node_x_indx >= 0 && node_y_indx >= 0 && node_z_indx >= 0 && node_x_indx<domain_3D->nx && node_y_indx<domain_3D->ny && node_z_indx<domain_3D->nz)
	{
		cell_indx = node_x_indx + domain_3D->nx*node_y_indx + domain_3D->nx*domain_3D->ny*node_z_indx;
		
		*cell_indx_local = ( ( cell_indx % (domain_3D->nx*domain_3D->by) ) / domain_3D->nx)*(domain_3D->nx/domain_3D->dx);
		*cell_indx_local += cell_indx % (domain_3D->nx/domain_3D->dx);
		*cell_indx_local += ( (cell_indx/(domain_3D->nx*domain_3D->ny))*(domain_3D->nx/domain_3D->dx*domain_3D->by) )% domain_3D->bxyz;
		
		plane_z = (INDEX)cell_indx/((domain_3D->bx*domain_3D->by)*domain_3D->dx*domain_3D->dy);
		plane_x = (INDEX)cell_indx%domain_3D->nx;
		plane_y = (INDEX)(cell_indx%((domain_3D->bx*domain_3D->by)*domain_3D->dx*domain_3D->dy ))/domain_3D->nx;
		
		*block_idx = (INDEX)(plane_z/domain_3D->bz)*domain_3D->dx*domain_3D->dy;
		*block_idx += (INDEX)(plane_y/domain_3D->by)*domain_3D->dx;
		*block_idx += (INDEX)(plane_x/domain_3D->bx);
		
		cell_indx_global = *cell_indx_local+(*block_idx)*domain_3D->by*domain_3D->bx*domain_3D->bz;

	}else
	{
		cell_indx_global = INDEX_MAX;
	}

	return cell_indx_global;
}

void cel_domain_3D_get_node_coordinate_indicies(  cel_domain_3D* domain_3D, INDEX cell_indx, INDEX* node_x_indx, INDEX* node_y_indx, INDEX* node_z_indx)
{
	INDEX block_indx, temp_indx;
	INDEX plane_x, plane_y, plane_z;
	INDEX b_plane_x, b_plane_y, b_plane_z;
	
	block_indx = cell_indx / domain_3D->bxyz;
	temp_indx = cell_indx-block_indx*domain_3D->bxyz;
	
	plane_z = (INDEX) block_indx/(domain_3D->dx*domain_3D->dy);
	plane_y = (INDEX)(block_indx-(plane_z)*domain_3D->dx*domain_3D->dy)/domain_3D->dx;
	plane_x =(INDEX)(block_indx-(plane_z)*domain_3D->dx*domain_3D->dy-(plane_y)*domain_3D->dx);
	
	b_plane_z = (INDEX) temp_indx/(domain_3D->bx*domain_3D->by);
	b_plane_y = (INDEX)(temp_indx-(b_plane_z)*domain_3D->bx*domain_3D->by)/domain_3D->bx;
	b_plane_x =(INDEX)(temp_indx-(b_plane_z)*domain_3D->bx*domain_3D->by-(b_plane_y)*domain_3D->bx);
	
	*node_x_indx = plane_x*domain_3D->bx+b_plane_x;
	*node_y_indx = plane_y*domain_3D->by+b_plane_y;
	*node_z_indx = plane_z*domain_3D->bz+b_plane_z;
}

void cel_domain_3D_get_block_boundary(  cel_domain_3D* domain_3D, INDEX block_indx, INDEX* xx_start, INDEX* yy_start, INDEX* zz_start,INDEX* xx_end, INDEX* yy_end, INDEX* zz_end)
{
	INDEX plane_x, plane_y, plane_z;
	
	plane_z = (INDEX) block_indx/(domain_3D->dx*domain_3D->dy);
	plane_y = (INDEX)(block_indx-plane_z*domain_3D->dx*domain_3D->dy)/domain_3D->dx;
	plane_x =(INDEX)(block_indx-plane_z*domain_3D->dx*domain_3D->dy-plane_y*domain_3D->dx);
	
	*zz_start = plane_z*domain_3D->bz;
	*yy_start = plane_y*domain_3D->by;
	*xx_start = plane_x*domain_3D->bx;
	
	*zz_end = *zz_start+domain_3D->bz;
	*yy_end = *yy_start+domain_3D->by;
	*xx_end = *xx_start+domain_3D->bx;
}

