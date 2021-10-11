!     
! File:   cel_perf_cg_module.f90
! Project CRESTA (see details on https://cresta-project.eu) Exascale library
! Send you email about the bugs to
! The High Performance Computing Center Stuttgart (HLRS) of the University of Stuttgart
! Autor: Dmitry Khabi,  (email: khabi@hlrs.de)
! Created on Okt 29, 2013
!

module  cel_perf_cg_module
!!use constants_module
use cel_types_module
use cel_cpu_param_module
use cel_perf_module
use OMP_LIB
implicit none

type(cel_perf_sum_type)   :: cel_perf_cg_sum
type(cel_perf_type)   :: cel_perf_cg
type(cel_perf_type)   :: cel_perf_bcast_comm

logical(kind=lk),parameter :: cel_perf_on = .false.

integer(kind=ik) :: cel_perf_mv_test_counter = 0_ik
integer(kind=ik) :: cel_perf_vect_test_counter = 0_ik
integer(kind=ik) :: cel_perf_mv_mpi_test_counter = 0_ik
integer(kind=ik) :: cel_perf_vec_mpi_test_counter = 0_ik


integer(kind=omp_lock_kind) :: cel_perf_counter_master_lock
integer(kind=omp_lock_kind) :: cel_perf_counter_worker_lock

integer(kind=ik),parameter :: i_type_mv=1
character(len=64),parameter   :: i_type_mv_label=TRIM("matrix_vector_mult")
character(len=256),parameter   :: i_type_mv_descr=TRIM("matrix_vector_mult") // &
  TRIM( "; platform: ")// TRIM(cel_cpu_param_description)
integer(kind=ik),parameter :: i_type_vect=2
character(len=64),parameter   :: i_type_vect_label=TRIM("vector_op")
character(len=256),parameter   :: i_type_vect_descr=TRIM("vector_op")// &
  TRIM( "; platform: ")// TRIM(cel_cpu_param_description)
integer(kind=ik),parameter :: i_type_mv_mpi=3
character(len=64),parameter   :: i_type_mv_mpi_label=TRIM("mpi_matrix_vector_mult")
character(len=256),parameter   :: i_type_mv_mpi_descr=TRIM("mpi_matrix_vector_mult")// &
  TRIM( "; platform: ")// TRIM(cel_cpu_param_description)
integer(kind=ik),parameter :: i_type_vec_mpi=4
character(len=64),parameter   :: i_type_vec_mpi_label=TRIM("mpi_vector_op")
character(len=256),parameter   :: i_type_vec_mpi_descr=TRIM("mpi_vector_op")// &
  TRIM( "; platform: ")// TRIM(cel_cpu_param_description)
integer(kind=ik),parameter :: i_type_max=4


end module  cel_perf_cg_module
!***********************************************************************

