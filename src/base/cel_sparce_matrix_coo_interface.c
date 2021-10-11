/*     
 * File:   cel_test_petsc.c
 * Project CRESTA (see details on https://cresta-project.eu)
 * Send you email about the bugs to
 * The High Performance Computing Center Stuttgart (HLRS) of the University of Stuttgart
 * Dmitry Khabi,  (email: khabi@hlrs.de)
 * Created on March 25, 2013
*/

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include "cel_types.h"
#include "cel_vector.h"
#include "cel_sparse_matrix_coo.h"
#include "cel_sparse_matrix_csr.h"
#include "cel_discretisation.h"
#include "cel_timer.h"

int cel_sparce_matrix_coo_interface(int* rank, cel_discretisation_task* desc_algorithm, cel_domain_3D* domain_3D,\
 cel_sparse_matrix_coo* sparse_matrix_coo, cel_vector* boundary_vector,\
 cel_vector* th_solution_vector, int* verbose)
{

	INDEX block_indx;
	cel_vector* cel_vectors;
	cel_vector* x_vector;
	int my_rank = *rank;
        int err;
  int verbose_level = *verbose;
  //if(verbose_level>0)
  //  printf("rank = %d\n",my_rank);

	cel_vectors = cel_vector_new(1,0);
	x_vector = cel_vectors;

	if(cel_domain_3D_init(domain_3D->nx, domain_3D->ny, domain_3D->nz, domain_3D->dx, domain_3D->dy, domain_3D->dz, domain_3D, verbose_level)!=0)
		return 3;

	block_indx = my_rank;
  
  if(verbose_level>0)
    cel_domain_3D_print(domain_3D,1);

	if((err=cel_discretisation(desc_algorithm, domain_3D, block_indx, sparse_matrix_coo, verbose_level))!=0)
	{
		return err;
	}

	if(cel_discretisation_boundary_conditions(desc_algorithm, domain_3D, block_indx, sparse_matrix_coo, boundary_vector, x_vector, th_solution_vector, verbose_level)!=0)
	{
		return 6;
	}
	cel_vector_free(cel_vectors, 0);
	free(cel_vectors);

	return 0;
}

