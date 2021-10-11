!     
! File:   cel_cg_module.f90
! Project CRESTA (see details on https://www.cresta-project.eu) Exascale library
! Send you email about the bugs to
! The High Performance Computing Center Stuttgart (HLRS) of the University of Stuttgart
! Dmitry Khabi,  (email: khabi@hlrs.de)
! Created on Jul 5, 2013
!

module cel_cg_module
use cel_types_module
use cel_error_module
use cel_vec_module
use OMP_LIB
implicit none

type cel_cg_type
  type(cel_vec_type)                            :: vector_x
  type(cel_vec_type)                            :: vector_r
  type(cel_vec_type)                            :: vector_d
  type(cel_vec_type)                            :: vector_q
  type(cel_vec_type)                            :: vector_tmp
  real(kind=rk)                                 :: sigma_new
  real(kind=rk)                                 :: sigma_old
  real(kind=rk)                                 :: betta
  real(kind=rk)                                 :: alpha
  real(kind=rk)                                 :: scalar_d_q
  real(kind=rk)                                 :: real_tmp
  real(kind=rk)                                 :: norm2
  real(kind=rk)                                 :: eps
  real(kind=rk)                                 :: residuum
  integer(kind=ik)                              :: iter_max
  integer(kind=ik)                              :: iter_done
  integer(kind=mpik)                            :: mpi_win_vector_x
  integer(kind=mpik)                            :: mpi_win_vector_d
  integer(kind=ik)                              :: verbosity
end type cel_cg_type
 
interface del
    module procedure cel_cg_del
end interface del

contains

subroutine cel_cg_del(cel_cg)
  type(cel_cg_type), intent(inout)  :: cel_cg
  call del(cel_cg%vector_r)
  call del(cel_cg%vector_x)
  call del(cel_cg%vector_d)
  call del(cel_cg%vector_q)
end subroutine cel_cg_del

end module cel_cg_module
