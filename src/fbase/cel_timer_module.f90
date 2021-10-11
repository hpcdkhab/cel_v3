!
! File:   cel_timer_module.f90
! Project CRESTA (see details on https://cresta-project.eu) Exascale library
! Send you email about the bugs to
! The High Performance Computing Center Stuttgart (HLRS) of the University of Stuttgart
! Dmitry Khabi,  (email: khabi@hlrs.de)
! Created on May 29, 2013
!

module  cel_timer_module
!!use constants_module

implicit none

real(kind=8)      :: machine_cycle_time
real(kind=8)      :: time_per_clock_tic_in_seconds
real(kind=8)      :: old_time_stamp=0.

integer(kind=8),external  :: getrdtsc
integer(kind=8)           :: clock_start
!!real(kind=8),external  :: clock_tic  ! if clock_tic done via F95 cpu_time


contains
!***************************************************************************************************************/
!***********************************************************************
!!function get_machine_cycle_time_for(hostname) result(result_value)
!!character(len=*)   :: hostname
!!real(kind=8)       :: result_value
!!
!!select case(hostname)
!!
!!case('node0','node1','node2','node3','node4','node5')
!!   result_value=1.d0/565.d6  ! NEC SX-6
!!case('nfhpc32')
!!   result_value=1.d0/1400.d6  ! Pentium M
!!case('asama')
!!   result_value=1.d0/1500.d6  ! asama     
!!case('ia64')
!!   result_value=1.d0/1300.d6  ! ia64     
!!case default
!!   print*,'get_machine_cycle_time_for: no machine_cycle_time identified for hostname "',trim(hostname),'"'
!!   stop
!!end select
!!
!!end function get_machine_cycle_time_for
!***********************************************************************
function clock_tics_per_second()
real(kind=8) :: start,end
!!real         :: start,end
!!real(kind=8) :: start_counter,end_counter
integer(kind=8) :: start_counter,end_counter
real(kind=8) :: clock_tics_per_second
real(kind=8),external   :: omp_get_wtime


! Initialization of the global reference clock_start; in this way we decrease the possibility of overflow
clock_start =getrdtsc()

start=0.
end=0.
   !!call cpu_time(start)
   start=omp_get_wtime()
do 
   !!call cpu_time(end)
   end=omp_get_wtime()
   if( end > start ) exit
enddo

start_counter=getrdtsc()
start=end
! we measure a time intervall of 1.0 seconds by the fortran95 call cpu_time in seconds
! and at the same time by the call clock_tic
   do 
      !!call cpu_time(end)
      end=omp_get_wtime()
   if( end - start > 1. ) exit
   enddo
end_counter=getrdtsc()


clock_tics_per_second=real(end_counter-start_counter,kind=8)/real(end - start,kind=8)

time_per_clock_tic_in_seconds=1.d0/clock_tics_per_second

!!print*,'end_counter,start_counter=',end_counter,start_counter
!!print*,'end_counter-start_counter=',end_counter-start_counter
!!print*,'end-start=',end-start
!!print*,'clock_tics_per_second=',clock_tics_per_second

end function clock_tics_per_second
!***********************************************************************
function clock_tic_increment_in_seconds()
! calculates the increment of my_second in seconds
real(kind=8) :: start,end
real(kind=8) :: clock_tic_increment_in_seconds
real(kind=8) :: clock_tics_per_second_0

!!print*,'clock_tic_increment_in_seconds: begin'

! to ensure that the parameter time_per_clock_tic_in_seconds will be set
clock_tics_per_second_0=clock_tics_per_second()

!!print*,'clock_tics_per_second=',clock_tics_per_second_0

call my_second(start)
do 
call my_second(end)
if( end > start ) exit
enddo

start=end
do 
call my_second(end)
if( end > start ) exit
enddo

clock_tic_increment_in_seconds=end - start


end function clock_tic_increment_in_seconds
!***********************************************************************
subroutine my_second(second)
real(kind=8)      :: second

!!integer(kind=ik)   :: counter
real(kind=8)   :: counter


counter=getrdtsc()
!!print*,'counter=',counter

!!second=machine_cycle_time*(counter-clock_start)
second=time_per_clock_tic_in_seconds*(counter-clock_start)

if(old_time_stamp > second) then
   print*,'my_second: error'
   print*,'old_time_stamp=',old_time_stamp
   print*,'second=',second
   stop
endif
   old_time_stamp=second


end subroutine my_second
!***********************************************************************
subroutine test_routine
real(kind=8) :: gg

gg=clock_tics_per_second()
print*,'clock_tics_per_second=',gg

end subroutine test_routine
!***********************************************************************
end module  cel_timer_module
!***********************************************************************
!!program test_timer
!!use timer_module
!!call test_routine
!!end program test_timer
