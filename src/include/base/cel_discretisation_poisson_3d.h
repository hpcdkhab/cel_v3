/*     
 * File:   cel_discretisation_poisson_3d.h
 * Project CRESTA (see details on https://cresta-project.eu) Exascale library
 * Send you email about the bugs to
 * The High Performance Computing Center Stuttgart (HLRS) of the University of Stuttgart
 * Dmitry Khabi,  (email: khabi@hlrs.de)
 * Created on March 25, 2013
*/
#ifndef CEL_DISCRETISATION_POISSON_3D_H
#define CEL_DISCRETISATION_POISSON_3D_H
#include "cel_types.h"
#include "cel_domain_3d.h"
#include "cel_discretisation.h"
#include "cel_sparse_matrix_coo.h"
#include "cel_vector.h"

//! Do finite difference discretisation for poisson problem in regular grid on domain [0;1]x[0;1]x[0;1]
//! @param (in) desc_task pointer to structure with task description
//! @param (in) domain_3D pointer to structure with decomposition's parameter of structured grid on 3D domain
//! @param (in) block_indx index of sub block to discretisate
//! @param (in) sparse_matrix_coo pointer to \struct sparse matrix in coordinate format (COO)
//! @param (in) output_on 1 - print error to sdterr; 0 - no output
//! @return  0 - no error ; -1 - not enought memory
int cel_discretisation_poisson_3D(cel_discretisation_task* disc_task, cel_domain_3D* domain_3D, INDEX block_indx, cel_sparse_matrix_coo* sparse_matrix_coo, int output_on);

//! initialize vector b with boundary conditions
//! @param (in) desc_task pointer to structure with task description
//! @param (in) domain_3D pointer to structure of decomposition's parameter set of structured grid on 3D domain
//! @param (in) block_indx index of block to discretisate
//! @param (out) boundary_vector pointer to \struct of vector with boundary conditions (RHS of linear system)
//! @param (out) x_vector pointer to \struct of start vector for the solver
//! @param (in) output_on 1 - print error to sdterr; 0 - no output
//! @return  0 - no error ; -1 - not enought memory
int cel_discretisation_poisson_3D_boundary_conditions(cel_discretisation_task* desc_task, cel_domain_3D* domain_3D, INDEX block_indx, cel_vector* boundary_vector, cel_vector* x_vector, cel_vector* solution_vector, int output_on);

//! Eliminate dirichlet values
//! @param (in) desc_task pointer to structure with task description
//! @param (in) domain_3D pointer to structure with decomposition's parameter of structured grid on 3D domain
//! @param (out) boundary_vector pointer to rhs
//! @param (out) matrix_coo pointer to sparse matrix in coo format
//! @param (in) output_on 1 - print error to sdterr; 0 - no output
//! @return  0 - no error ; -1 - not enought memory
int cel_discretisation_poisson_3D_eliminate_dirichlet_values(cel_discretisation_task* desc_task, cel_domain_3D* domain_3D, cel_vector* boundary_vector,cel_sparse_matrix_coo* sparse_matrix_coo, int output_on);

//! get discretization error (integral of error over domain)
//! @param (in) desc_task pointer to /struct with task description
//! @param (in) domain_3D pointer to /struct with decomposition's parameter of structured grid on 3D domain
//! @param (in) x_vector pointer to solution
//! @return  error of discretisation
REAL cel_discretisation_poisson_3D_error(cel_discretisation_task* desc_task, cel_domain_3D* domain_3D, cel_vector* x_vector);

void _poisson_3D_node_discr_7_stencil(cel_discretisation_task* desc_task, cel_domain_3D* domain_3D, INDEX ii, INDEX jj, INDEX kk, INDEX* rows, INDEX* cols, REAL* values, INDEX *current_indx, INDEX* rows_num);
void _poisson_3D_node_discr_19_stencil(cel_discretisation_task* desc_task, cel_domain_3D* domain_3D, INDEX ii, INDEX jj, INDEX kk, INDEX* rows, INDEX* cols, REAL* values, INDEX *current_indx, INDEX* rows_num);
void _poisson_3D_node_discr_27_stencil(cel_discretisation_task* desc_task, cel_domain_3D* domain_3D, INDEX ii, INDEX jj, INDEX kk, INDEX* rows, INDEX* cols, REAL* values, INDEX *current_indx, INDEX* rows_num);
REAL _poisson_3D_solution(cel_discretisation_task*  desc_task, REAL xx, REAL yy, REAL zz);
REAL _poisson_3D_rhs_boundary(cel_discretisation_task*  desc_task, REAL xx, REAL yy, REAL zz);
REAL _poisson_3D_rhs_inner(cel_discretisation_task*  desc_task, REAL xx, REAL yy, REAL zz);

#endif

