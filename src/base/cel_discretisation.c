/*     
 * File:   cel_discretisation.c -- see headers for details and function's comments
 * Project CRESTA (see details on https://cresta-project.eu)
 * CRESTA Exascale library
 * Send you email about the bugs to
 * The High Performance Computing Center Stuttgart (HLRS) of the University of Stuttgart
 * Dmitry Khabi,  (email: khabi@hlrs.de)
 * Created on March 25, 2013
*/

#include "cel_discretisation.h"
#include "cel_discretisation_poisson_3d.h"
#include <stdio.h>
int cel_discretisation(cel_discretisation_task* desc_task, cel_domain_3D* domain_3D, INDEX block_indx, cel_sparse_matrix_coo* sparse_matrix_coo, int output_on)
{
        int err_code;
        err_code = 0;
        switch (desc_task->problem)
        {
                case POISSON_3D:
                        err_code = cel_discretisation_poisson_3D(desc_task, domain_3D, block_indx, sparse_matrix_coo, output_on);
                        break;
                default:
                        err_code = 2;
						if(output_on>0)
							fprintf(stderr,"ERROR in cel_discretisation - disretisation problem is not defined!\n");
                        break;
        }
	
        if(err_code > 0 )
		{
			if(output_on>0)
				fprintf(stderr,"ERROR in cel_discretisation\n");
		}
        return err_code;
}

int cel_discretisation_boundary_conditions(cel_discretisation_task* desc_task, cel_domain_3D* domain_3D, INDEX block_indx, cel_sparse_matrix_coo* sparse_matrix_coo, cel_vector* boundary_vector, cel_vector* x_vector, cel_vector* solution_vector, int output_on)
{
         int err_code;
        err_code = 0;
        switch (desc_task->problem)
        {
                case POISSON_3D:
                {
                        switch (desc_task->boundary_cond)
                        {
                                case X_2_PLUS_Y_2_PLUS_Z_2:
                                        err_code += cel_discretisation_poisson_3D_boundary_conditions(desc_task,domain_3D,block_indx,boundary_vector,x_vector,solution_vector,output_on);
                                        err_code += cel_discretisation_poisson_3D_eliminate_dirichlet_values(desc_task,domain_3D,boundary_vector,sparse_matrix_coo,output_on);
                                        break;
                                default:
                                        err_code = 3;
										if(output_on>0)
											fprintf(stderr,"ERROR in cel_discretisation_boundary_conditions - boundary conditions are not defined!\n");
										break;
                        }
                        break;
                }
                default:
                        err_code = 2;
		        if(output_on>0)
			        fprintf(stderr,"ERROR in cel_discretisation_boundary_conditions - disretisation problem is not defined !\n");
                        break;
        }
	
        if(err_code > 0 )
		{
			if(output_on>0)
				fprintf(stderr,"ERROR in el_discretisation_boundary_conditions\n");
		}
        return err_code;
}

