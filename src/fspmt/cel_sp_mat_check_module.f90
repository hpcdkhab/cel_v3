!     
! File:  cel_blspmt_cel.f90
! Project CRESTA (see details on https://cresta-project.eu) Exascale library
! Send you email about the bugs to
! The High Performance Computing Center Stuttgart (HLRS) of the University of Stuttgart
! Dmitry Khabi,  (email: khabi@hlrs.de)
! Created on Aug 27, 2013
!

module cel_sp_mat_check_module
!library fbase
use cel_types_module
use cel_input_parameters_module
use cel_omp_shared_work_module
use cel_omp_module
use cel_perf_module
use cel_timer_interface_module
use cel_sp_mat_module
!use cel_cg_module
!library fspmt
use cel_vec_module
use cel_sp_mat_module
!library fblspmt

implicit none

contains


subroutine cel_sp_mat_check_diagonal_dominance(a_sp_mat,diagonal,&
                                sum_off_diagonal,diagonal_dominance,&
                                vtxdist,cel_omp,cel_omp_shared_work,&
                                cel_perf_counter,input_parameters,err_code)
  type(cel_sp_mat_type), intent(inout)                              :: a_sp_mat
  real(kind=rk),dimension(:), intent(inout)                         :: diagonal
  real(kind=rk),dimension(:), intent(inout)                         :: sum_off_diagonal
  real(kind=rk),dimension(:), intent(inout)                         :: diagonal_dominance
  type(cel_omp_shared_work_type), intent(inout)                     :: cel_omp_shared_work
  integer(kind=ik),dimension(:), intent(in)                         :: vtxdist
  type(cel_omp_type), intent(in)                                    :: cel_omp
  type(cel_perf_counter_type), intent(inout)                        :: cel_perf_counter
  type(cel_input_parameters_type), intent(in)                       :: input_parameters
  integer(kind=ik), intent(out)                                     :: err_code
  integer(kind=ik) :: ii
  integer(kind=mpik) :: ierr
  integer(kind=ik) :: omp_start_idx,omp_end_idx,num_threads
  integer(kind=ik) :: total_number,col_offset,row_offset, vec_length
  real(kind=rk) :: value
  integer(kind=ik) :: col_idx, row_idx

  err_code = 0_ik
  !if(cel_omp%master_num .eq. 2) then
  !initialize diagonal with zero
  !divide work between threads
  vec_length = size(diagonal,dim=1)
  call cel_omp_shared_work_start(cel_omp_shared_work,cel_omp%num_threads)
  call cel_omp_shared_work_distr(cel_omp_shared_work,&
                                 omp_start_idx,omp_end_idx,num_threads,&
                                 vec_length,vec_length,&
                                 cel_omp,input_parameters%verbosity,err_code)
  if(cel_omp%worker_num .le. num_threads) then
    diagonal(omp_start_idx:omp_end_idx)=0.0_rk
    sum_off_diagonal(omp_start_idx:omp_end_idx)=0.0_rk
  end if
  call cel_omp_shared_work_end_and_wait(cel_omp_shared_work,cel_omp,LOOP_SLEEP_DEFLENGTH)
  
    
  !collect diagonal elements
  !divide work between threads
  call cel_omp_shared_work_start(cel_omp_shared_work,cel_omp%num_threads)
  call cel_omp_shared_work_distr(cel_omp_shared_work,&
                                 omp_start_idx,omp_end_idx,num_threads,&
                                 a_sp_mat%elements_num,a_sp_mat%elements_num,&
                                 cel_omp,input_parameters%verbosity,err_code)
  if(cel_omp%worker_num .le. num_threads) then
    total_number = (omp_end_idx-omp_start_idx+1_ik)
    col_offset =  a_sp_mat%col_offset
    row_offset =  a_sp_mat%row_offset
    !write(*,'(I0,A,I0,A,I0,A,I0, A,I0,A,I0)') cel_omp%worker_num, ": elems:",a_sp_mat%elements_num,&
    !  "; omp_start_idx:",omp_start_idx,"; omp_end_idx:",omp_end_idx,&
    !  "; col_offset:",col_offset,"; row_offset:",row_offset 
    do ii = omp_start_idx, omp_end_idx
      col_idx = a_sp_mat%cols(ii) - col_offset
      row_idx = a_sp_mat%rows(ii) - row_offset
      value = a_sp_mat%values(ii)
      if(col_idx .eq. row_idx) then
        diagonal(row_idx) = abs(value)
      else
        sum_off_diagonal(row_idx)=sum_off_diagonal(row_idx)+abs(value)
      end if
    end do
  end if
  call cel_omp_shared_work_end_and_wait(cel_omp_shared_work,cel_omp,LOOP_SLEEP_DEFLENGTH)
  call cel_omp_shared_work_start(cel_omp_shared_work,cel_omp%num_threads)
  call cel_omp_shared_work_distr(cel_omp_shared_work,&
                                 omp_start_idx,omp_end_idx,num_threads,&
                                 a_sp_mat%rows_num,a_sp_mat%rows_num,&
                                 cel_omp,input_parameters%verbosity,err_code)
   if(cel_omp%worker_num .le. num_threads) then
     do ii = omp_start_idx,omp_end_idx
       diagonal_dominance(ii)=diagonal(ii) - sum_off_diagonal(ii)
     end do
   end if                              
   call cel_omp_shared_work_end_and_wait(cel_omp_shared_work,cel_omp,LOOP_SLEEP_DEFLENGTH)
  !end if
end subroutine cel_sp_mat_check_diagonal_dominance


end module cel_sp_mat_check_module
