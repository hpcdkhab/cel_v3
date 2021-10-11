!
! File:   cel_omp_shared_work_module.f90
! Project CRESTA (see details on https://www.cresta-project.eu) Exascale library
! Send you email about the bugs to
! The High Performance Computing Center Stuttgart (HLRS) of the University of Stuttgart
! Dmitry Khabi,  (email: khabi@hlrs.de)
! Created on Nov 28, 2013
!

module cel_omp_shared_work_module
use cel_types_module
use cel_error_module
use cel_cpu_param_module
use cel_omp_module
use cel_timer_interface_module
use OMP_LIB
implicit none

  

  interface new
    module procedure cel_omp_shared_work_new
  end interface new

  interface del
    module procedure cel_omp_shared_work_del
  end interface del
 
  type cel_omp_shared_work_type
    integer(kind=ik)                                     :: num_working_threads
    integer(kind=tik)                                    :: loop_sleep_counter
    integer(kind=omp_lock_kind)                          :: omp_lock
    integer(kind=omp_lock_kind)                          :: counter
    integer(kind=omp_lock_kind),allocatable,dimension(:) :: go_local
    integer(kind=omp_lock_kind),allocatable,dimension(:) :: go
  end type cel_omp_shared_work_type

contains

subroutine cel_omp_shared_work_new(cel_omp_shared_work,cel_omp)
  type(cel_omp_shared_work_type), intent(inout) :: cel_omp_shared_work
  type(cel_omp_type), intent(in)                :: cel_omp

  if(cel_omp%is_master) then
    call omp_init_lock(cel_omp_shared_work%omp_lock)
    cel_omp_shared_work%num_working_threads = 0_ik
    cel_omp_shared_work%loop_sleep_counter = 0_tik
    cel_omp_shared_work%counter = 0
    allocate(cel_omp_shared_work%go_local(cel_omp%num_threads))
    allocate(cel_omp_shared_work%go(cel_omp%num_threads))
    cel_omp_shared_work%go_local = 1
    cel_omp_shared_work%go = 0
  end if
  
end subroutine cel_omp_shared_work_new

subroutine cel_omp_shared_work_del(cel_omp_shared_work)
  type(cel_omp_shared_work_type), intent(inout) :: cel_omp_shared_work
  integer(kind=ik) :: ii

  call omp_destroy_lock(cel_omp_shared_work%omp_lock)
  cel_omp_shared_work%num_working_threads = 0_ik
  cel_omp_shared_work%loop_sleep_counter = 0_tik
  cel_omp_shared_work%counter = 0
  if(allocated(cel_omp_shared_work%go_local))&
    deallocate(cel_omp_shared_work%go_local)
  if(allocated(cel_omp_shared_work%go))&
    deallocate(cel_omp_shared_work%go)
end subroutine cel_omp_shared_work_del

subroutine cel_omp_shared_work_start(cel_omp_shared_work,num_threads)
  type(cel_omp_shared_work_type), intent(inout) :: cel_omp_shared_work
  integer(kind=ik), intent(in)                  :: num_threads
  integer(kind=ik) :: num_locks

  call omp_set_lock (cel_omp_shared_work%omp_lock)
  if(cel_omp_shared_work%num_working_threads .ne. num_threads) then
    cel_omp_shared_work%num_working_threads = num_threads
  endif
  call omp_unset_lock (cel_omp_shared_work%omp_lock)
  
end subroutine cel_omp_shared_work_start

subroutine cel_omp_shared_work_end_and_wait(cel_omp_shared_work,cel_omp,loop_sleep_length)
  type(cel_omp_shared_work_type), intent(inout) :: cel_omp_shared_work
  type(cel_omp_type), intent(in)                :: cel_omp
  integer(kind=ik), intent(in)                  :: loop_sleep_length
  integer(kind=ik)  :: counter, go, ii
  
  call cel_omp_shared_work_go_(go,cel_omp_shared_work,cel_omp)
  cel_omp_shared_work%go_local(cel_omp%profile_indx) = go
  call cel_omp_shared_work_counter_inc_(counter,cel_omp_shared_work)
  if(counter .eq. cel_omp_shared_work%num_working_threads) then
    call omp_set_lock (cel_omp_shared_work%omp_lock)
    cel_omp_shared_work%counter = 0_ik
    do ii=1_ik,cel_omp%num_threads
      cel_omp_shared_work%go(ii) = 1_ik - cel_omp_shared_work%go(ii)
    end do
    call omp_unset_lock (cel_omp_shared_work%omp_lock)    
  else
   call cel_omp_shared_work_wait_go_(cel_omp_shared_work,cel_omp,loop_sleep_length)
  endif

end subroutine cel_omp_shared_work_end_and_wait

subroutine cel_omp_shared_work_go_(go, cel_omp_shared_work,cel_omp)
  integer(kind=ik), intent(out)                 :: go
  type(cel_omp_shared_work_type), intent(inout) :: cel_omp_shared_work
  type(cel_omp_type), intent(in)                :: cel_omp

  call omp_set_lock (cel_omp_shared_work%omp_lock)
  go = cel_omp_shared_work%go(cel_omp%profile_indx)
  call omp_unset_lock (cel_omp_shared_work%omp_lock)
  
end subroutine cel_omp_shared_work_go_

subroutine cel_omp_shared_work_go_go_(go,go_local, cel_omp_shared_work,cel_omp)
  integer(kind=ik), intent(out)                 :: go
  integer(kind=ik), intent(out)                 :: go_local
  type(cel_omp_shared_work_type), intent(inout) :: cel_omp_shared_work
  type(cel_omp_type), intent(in)                :: cel_omp

  call omp_set_lock (cel_omp_shared_work%omp_lock)
  go = cel_omp_shared_work%go(cel_omp%profile_indx)
  go_local = cel_omp_shared_work%go_local(cel_omp%profile_indx)
  call omp_unset_lock (cel_omp_shared_work%omp_lock)
  
end subroutine cel_omp_shared_work_go_go_

subroutine cel_omp_shared_work_counter_(counter, cel_omp_shared_work)
  integer(kind=ik), intent(out)                 :: counter
  type(cel_omp_shared_work_type), intent(inout) :: cel_omp_shared_work

  call omp_set_lock (cel_omp_shared_work%omp_lock)
  counter = cel_omp_shared_work%counter
  call omp_unset_lock (cel_omp_shared_work%omp_lock)
  
end subroutine cel_omp_shared_work_counter_

subroutine cel_omp_shared_work_counter_inc_(counter,cel_omp_shared_work)
  integer(kind=ik),intent(out)                   :: counter
  type(cel_omp_shared_work_type), intent(inout) :: cel_omp_shared_work

  call omp_set_lock (cel_omp_shared_work%omp_lock)
  cel_omp_shared_work%counter = cel_omp_shared_work%counter + 1_ik
  counter = cel_omp_shared_work%counter
  call omp_unset_lock (cel_omp_shared_work%omp_lock)
  
end subroutine cel_omp_shared_work_counter_inc_

subroutine cel_omp_shared_work_wait_go_(cel_omp_shared_work,cel_omp,loop_sleep_length)
  type(cel_omp_shared_work_type), intent(inout) :: cel_omp_shared_work
  type(cel_omp_type), intent(in)                :: cel_omp
  integer(kind=ik), intent(in)                  :: loop_sleep_length
  integer(kind=ik) :: ii, go,go_local, length
  integer(kind=itk) :: count_wait
  
  count_wait = 0_ik
  call cel_omp_shared_work_go_go_(go,go_local,cel_omp_shared_work,cel_omp)
  do while(go .eq. go_local)
    count_wait = count_wait + cel_omp_shared_work_loop_sleep(cel_omp_shared_work, loop_sleep_length)
    call cel_omp_shared_work_go_go_(go,go_local,cel_omp_shared_work,cel_omp)
  end do
  !if(count_wait .gt. 100000000_ik) then
  !   call cel_warning("bad thread load balancing(1)",count_wait, cel_is_in_debug)
  !end if
end subroutine cel_omp_shared_work_wait_go_

function cel_omp_shared_work_loop_sleep(cel_omp_shared_work,loop_sleep_length)&
  result(num_nops)
  type(cel_omp_shared_work_type), intent(inout) :: cel_omp_shared_work
  integer(kind=ik), intent(in)                  :: loop_sleep_length
  integer(kind=itk):: num_nops
  integer(kind=ik) :: ii

  num_nops = 0_ik
  do ii=1_ik, loop_sleep_length
    num_nops = num_nops + nops()
  end do

end function cel_omp_shared_work_loop_sleep

subroutine cel_omp_shared_work_distr(cel_omp_shared_work,&
                                     start_idx,end_indx,num_threads,&
                                     shared_work_length,shared_amount_of_work,&
                                     cel_omp,output_on,err_code)
  type(cel_omp_shared_work_type), intent(inout) :: cel_omp_shared_work
  integer(kind=ik), intent(out)                 ::start_idx
  integer(kind=ik), intent(out)                 ::end_indx
  integer(kind=ik), intent(out)                 ::num_threads
  integer(kind=ik), intent(in)                  ::shared_work_length
  integer(kind=ik), intent(in)                  ::shared_amount_of_work
  type(cel_omp_type), intent(in)                ::cel_omp
  logical(kind=lk), intent(in)                  ::output_on
  integer(kind=ik), intent(out)                 ::err_code
  integer(kind=ik) :: aligned_work_length
  integer(kind=ik) :: delta_work_length, aligned_num_blocks
  integer(kind=ik) :: num_free_threads, aligned_num_blocks_thread
  integer(kind=ik) :: delta_align_blocks
  err_code = 0_ik

  aligned_work_length = shared_work_length / AW_SMALLEST_OMP_TASK
  aligned_work_length = aligned_work_length * AW_SMALLEST_OMP_TASK
  delta_work_length = shared_work_length - aligned_work_length
  aligned_num_blocks = aligned_work_length / AW_SMALLEST_OMP_TASK

  num_free_threads = cel_omp%num_workers
  aligned_num_blocks_thread = aligned_num_blocks/num_free_threads
  aligned_num_blocks_thread = aligned_num_blocks_thread
  delta_align_blocks = aligned_num_blocks-aligned_num_blocks_thread*num_free_threads

  start_idx = aligned_num_blocks_thread*AW_SMALLEST_OMP_TASK*(cel_omp%worker_num-1_ik) + 1_ik
  end_indx = start_idx + aligned_num_blocks_thread*AW_SMALLEST_OMP_TASK - 1_ik
  end_indx = end_indx + cel_omp%worker_num/num_free_threads*&
        (delta_work_length + delta_align_blocks*AW_SMALLEST_OMP_TASK)

  num_threads = num_free_threads

  start_idx = max(0_ik,start_idx)
  end_indx = max(-1_ik,end_indx)
  !write(cel_output_unit,'(9(A,I0))') "work_length:",shared_work_length,&
  !"; num_threads:",num_threads,"; start_idx:",start_idx,&
  !": end_indx:",end_indx,&
  !"; w_num:", cel_omp%worker_num,&
  !"; num_bl:",aligned_num_blocks,&
  !"; num_bl_thread:",aligned_num_blocks_thread,&
  !": delta_length:",delta_work_length,&
  !": delta_ablocks:",delta_align_blocks

end subroutine cel_omp_shared_work_distr

end module cel_omp_shared_work_module
