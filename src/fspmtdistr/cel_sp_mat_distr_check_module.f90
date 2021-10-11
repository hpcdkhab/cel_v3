!     
! File:  cel_blspmt_cel.f90
! Project CRESTA (see details on https://cresta-project.eu) Exascale library
! Send you email about the bugs to
! The High Performance Computing Center Stuttgart (HLRS) of the University of Stuttgart
! Dmitry Khabi,  (email: khabi@hlrs.de)
! Created on Aug 27, 2013
!

module cel_sp_mat_distr_check_module
!library fbase
use cel_types_module
use cel_input_parameters_module
use cel_omp_module
use cel_comm_module
use cel_timer_interface_module
use cel_sp_mat_distr_module
use cel_sp_mat_distr_vec_module
use cel_sp_mat_check_module
!use cel_cg_module
!library fspmt
use cel_vec_module
use cel_sp_mat_module
!library fblspmt

!external libraries
use MPI
implicit none

contains

subroutine cel_sp_mat_distr_check_diagonal_dominance(a_sp_mats,distr_diagonal,&
                                distr_sum_off_diagonal,distr_dominance,&
                                cel_comm,vtxdist,cel_omp,cel_omp_shared_work,&
                                cel_perf_counter,input_parameters,err_code)
  type(cel_sp_mat_type), dimension(:), allocatable,intent(inout)    :: a_sp_mats
  type(cel_sp_mat_distr_vec_type), intent(out)                      :: distr_diagonal
  type(cel_sp_mat_distr_vec_type), intent(out)                      :: distr_sum_off_diagonal
  type(cel_sp_mat_distr_vec_type), intent(out)                      :: distr_dominance
  type(cel_comm_type), intent(inout)                                :: cel_comm
  type(cel_omp_shared_work_type), intent(inout)                     :: cel_omp_shared_work
  integer(kind=ik),dimension(:), intent(in)                         :: vtxdist
  type(cel_omp_type), intent(in)                                    :: cel_omp
  type(cel_perf_counter_type), intent(inout)                        :: cel_perf_counter
  type(cel_input_parameters_type)  , intent(in)                     :: input_parameters
  integer(kind=ik), intent(out)                                     :: err_code
  integer(kind=ik) :: ii
  integer(kind=mpik) :: ierr
  integer(kind=ik) :: length
  
  err_code = 0_ik
  !allocate distributed vector for matrix diagonals
  call cel_omp_shared_work_start(cel_omp_shared_work,cel_omp%num_threads)
  call cel_sp_mat_distr_vec_new_from_procs&
    (distr_diagonal,cel_comm%get_from,vtxdist,cel_omp)
  call cel_sp_mat_distr_vec_new_from_procs&
    (distr_sum_off_diagonal,cel_comm%get_from,vtxdist,cel_omp)
  call cel_sp_mat_distr_vec_new_from_procs&
    (distr_dominance,cel_comm%get_from,vtxdist,cel_omp)
  call cel_omp_shared_work_end_and_wait(cel_omp_shared_work,cel_omp,LOOP_SLEEP_DEFLENGTH)
  length = distr_diagonal%num_elements
  !collect diagonal for each matrix block
  !diagonal of the first matrix block is the global matrix diagonal
  do ii=1_ik, min(1,size(a_sp_mats))
    call cel_sp_mat_check_diagonal_dominance(&
       a_sp_mats(ii),distr_diagonal%values(1_ik:length,ii),&
       distr_sum_off_diagonal%values(1_ik:length,ii),&
       distr_dominance%values(1_ik:length,ii),&
       vtxdist,cel_omp,cel_omp_shared_work,&
       cel_perf_counter,input_parameters,err_code) 
  end do

end subroutine cel_sp_mat_distr_check_diagonal_dominance


end module cel_sp_mat_distr_check_module
