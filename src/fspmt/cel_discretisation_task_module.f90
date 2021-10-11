!     
! File:   cel_domain_3d_module.f90
! Project CRESTA (see details on https://cresta-project.eu) Exascale library
! Send you email about the bugs to
! The High Performance Computing Center Stuttgart (HLRS) of the University of Stuttgart
! Dmitry Khabi,  (email: khabi@hlrs.de)
! Created on May 29, 2013
!

module cel_discretisation_task_module
use cel_types_module
use, intrinsic :: ISO_C_BINDING
implicit none

integer(c_int), parameter :: POISSON_3D = 1
integer(c_int), parameter :: STENCIL_7 = 1
integer(c_int), parameter :: STENCIL_19 = 2
integer(c_int), parameter :: STENCIL_27 = 3
integer(c_int), parameter :: USER_DEFINED = 0
integer(c_int), parameter :: X_2_PLUS_Y_2_PLUS_Z_2 = 1

! struct choose of decretisation, interface to c
type, bind(c) :: cel_c_discretisation_task_type
  integer(kind=c_int) :: problem! problem
  integer(kind=c_int) :: stencil! problem
  integer(kind=c_int) :: boundary_cond! problem
end type cel_c_discretisation_task_type

! struct choose of decretisation
type :: cel_discretisation_task_type
  integer(kind=c_enum) :: problem! problem
  integer(kind=c_enum) :: stencil! problem
  integer(kind=c_enum) :: boundary_cond! problem
end type cel_discretisation_task_type

interface cel_convert_c
  module procedure cel_discretisation_task_c_to_f
end interface cel_convert_c

contains

subroutine cel_discretisation_task_c_to_f(cel_c_discretisation_task, cel_discretisation_task, to_deallocate_c)
  type(cel_c_discretisation_task_type), intent(inout) :: cel_c_discretisation_task
  type(cel_discretisation_task_type), intent(inout) :: cel_discretisation_task
  logical(kind=lk), intent(in) :: to_deallocate_c
  
  cel_discretisation_task%problem = INT(cel_c_discretisation_task%problem,kind=c_enum)
  cel_discretisation_task%stencil = INT(cel_c_discretisation_task%stencil,kind=c_enum)
  cel_discretisation_task%boundary_cond = INT(cel_c_discretisation_task%boundary_cond,kind=c_enum)
  
end subroutine cel_discretisation_task_c_to_f

end module cel_discretisation_task_module
