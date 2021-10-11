!     
! File:   cel_c_sparce_matrix_coo_interface_module.f90
! Project CRESTA (see details on https://cresta-project.eu) Exascale library
! Send you email about the bugs to
! The High Performance Computing Center Stuttgart (HLRS) of the University of Stuttgart
! Dmitry Khabi,  (email: khabi@hlrs.de)
! Created on May 29, 2013
!

module cel_c_sp_mat_interface_module

use, intrinsic :: ISO_C_BINDING
implicit none

interface
 integer(kind=C_INT) function cel_c_sp_mat_interface &
        (rank, discretisation_task, domain_3D, sparse_matrix_coo, boundary_vector, &
         th_solution_vector, verbose_level) &
   bind(C, name='cel_sparce_matrix_coo_interface')
   use, intrinsic :: iso_c_binding
use cel_types_module
use cel_discretisation_task_module
use cel_domain_3d_module
use cel_sp_mat_module
use cel_vec_module
   implicit none
   integer(kind=c_int) :: rank
   type(cel_c_discretisation_task_type) :: discretisation_task
   type(cel_c_domain_3D_type) :: domain_3D
   type(cel_c_sp_mat_type) :: sparse_matrix_coo
   type(cel_c_vec_type) :: boundary_vector
   type(cel_c_vec_type) :: th_solution_vector
   integer(kind=c_int) :: verbose_level

 end function cel_c_sp_mat_interface
end interface


end module cel_c_sp_mat_interface_module
