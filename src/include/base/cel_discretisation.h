/*     
 * File:   cel_discretisation.h
 * Project CRESTA (see details on https://cresta-project.eu) Exascale library
 * Send you email about the bugs to
 * The High Performance Computing Center Stuttgart (HLRS) of the University of Stuttgart
 * Dmitry Khabi,  (email: khabi@hlrs.de)
 * Created on March 25, 2013
*/
#ifndef CEL_DISCRETISATION_H
#define CEL_DISCRETISATION_H
#include "cel_types.h"
#include "cel_domain_3d.h"
#include "cel_vector.h"
#include "cel_sparse_matrix_coo.h"
//! \enum problems
enum cel_task_problem
{
	POISSON_3D = 1
};
//! \enum discretisation stencils
enum cel_task_stencil
{
	STENCIL_7 = 1,
	STENCIL_19 = 2,
	STENCIL_27 = 3
};
//! \enum boundary conditions
enum cel_task_boundary_cond
{
        USER_DEFINED = 0,
	X_2_PLUS_Y_2_PLUS_Z_2 = 1
};


//! \struct choose of decretisation
typedef struct 
{
	enum cel_task_problem problem;//! problem
	enum cel_task_stencil stencil;//! discretisation stencil
	enum cel_task_boundary_cond boundary_cond;//! boundary conditions
} cel_discretisation_task;

//! Do finite difference discretisation in regular grid on domain [0;1]x[0;1]x[0;1]
//! @param (in) desc_algorithm pointer to structure with task description
//! @param (in) domain_3D pointer to \struct with decomposition's parameter of structured grid on 3D domain
//! @param (in) block_indx index of sub block to discretisate
//! @param (in) sparse_matrix_coo pointer to \struct sparse matrix in coordinate format (COO)
//! @param (in) output_on 1 - print error to sdterr; 0 - no output
//! @return  0 - no error ; >1 - error
int cel_discretisation(cel_discretisation_task* desc_task, cel_domain_3D* domain_3D, INDEX block_indx, cel_sparse_matrix_coo* sparse_matrix_coo, int output_on);

//! Initialize vector b with boundary conditions
//! @param (in) desc_task pointer to \struct with task description
//! @param (in) domain_3D pointer to \struct of decomposition's parameter set of structured grid on 3D domain
//! @param (in) block_indx index of sub block to discretisate
//! @param (in) sparse_matrix_coo pointer to \struct sparse matrix in coordinate format (COO)
//! @param (out) boundary_vector pointer to \struct of vector with boundary conditions (RHS of linear system)
//! @param (out) x_vector pointer to \struct of start vector for the solver
//! @param (in) output_on 1 - print error to sdterr; 0 - no output
//! @return  0 - no error ; >1 - error
int cel_discretisation_boundary_conditions(cel_discretisation_task* desc_task, cel_domain_3D* domain_3D, INDEX block_indx, cel_sparse_matrix_coo* sparse_matrix_coo, cel_vector* boundary_vector, cel_vector* x_vector, cel_vector* solution_vector, int output_on);
#endif

