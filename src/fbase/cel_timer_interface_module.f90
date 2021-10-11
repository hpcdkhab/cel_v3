!
! File:   cel_timer_interface_module.f90
! Project CRESTA (see details on https://cresta-project.eu) Exascale library
! Send you email about the bugs to
! The High Performance Computing Center Stuttgart (HLRS) of the University of Stuttgart
! Dmitry Khabi,  (email: khabi@hlrs.de)
! Created on May 29, 2013
!

!===============
!Module with measure tools routines
!===============
MODULE cel_timer_interface_module
 
 implicit none
  integer, parameter :: tik = 8
  integer, parameter :: trk = 8
 interface
!  function read_tsc_register() result(cycles) bind(c, name = 'read_tsc_register')
!   use, intrinsic :: iso_c_binding
!   integer(kind = c_long) :: cycles
!  end function read_tsc_register

  function getrdtsc() result(cycles) bind(c, name = 'getrdtsc')
   use, intrinsic :: iso_c_binding
   integer(kind = c_long) :: cycles
  end function getrdtsc
  
  function nops() result(num_nops) bind(c, name = 'nops')
   use, intrinsic :: iso_c_binding
   integer(kind = c_long) :: num_nops
  end function nops

  function get_time() result(seconds) bind(c, name = 'get_time')
   use, intrinsic :: iso_c_binding
   real(kind = c_double) :: seconds
  end function get_time


  subroutine micro_sleep( micro_seconds) BIND(c,name="micro_sleep")
    use, INTRINSIC :: ISO_C_BINDING
    integer(kind = c_long) :: micro_seconds
  end subroutine micro_sleep
  
  subroutine get_time_stamp(seconds,nanoseconds) BIND(c,name="get_time_stamp")
    use, INTRINSIC :: ISO_C_BINDING
    integer(kind = c_long_long) :: seconds
    integer(kind = c_long_long) :: nanoseconds
  end subroutine get_time_stamp
  
  function get_tsc_cycle(accuracy,maximumDuration,output_on) &
  result(tsc_cycles) bind(c, name = 'get_tsc_cycle')
   use, intrinsic :: iso_c_binding
   real(kind = c_double),intent(in) :: accuracy
   real(kind = c_double),intent(in) :: maximumDuration
   integer(kind = c_int), intent(in) :: output_on
   real(kind = c_double) :: tsc_cycles
  end function get_tsc_cycle

  function craypat_off() result(ierr) bind(c, name = 'craypat_off')
   use, intrinsic :: iso_c_binding
   integer(kind = c_long) :: ierr
  end function craypat_off 

  function craypat_on() result(ierr) bind(c, name = 'craypat_on')
   use, intrinsic :: iso_c_binding
   integer(kind = c_long) :: ierr
  end function craypat_on

  function vtrace_off() result(ierr) bind(c, name = 'vtrace_off')
   use, intrinsic :: iso_c_binding
   integer(kind = c_long) :: ierr
  end function vtrace_off 

  function vtrace_on() result(ierr) bind(c, name = 'vtrace_on')
   use, intrinsic :: iso_c_binding
   integer(kind = c_long) :: ierr
  end function vtrace_on


  subroutine vt_start_1() BIND(c,name="vt_start_1")
    use, INTRINSIC :: ISO_C_BINDING
  end subroutine vt_start_1

  subroutine vt_end_1() BIND(c,name="vt_end_1")
    use, INTRINSIC :: ISO_C_BINDING
  end subroutine vt_end_1

   integer(kind=C_INT) function cel_free_double_c(array) &
     bind(C, name='cel_free_double_c')
     use, intrinsic :: iso_c_binding
     implicit none
     real(kind = c_double) ::array
   end function cel_free_double_c

   integer(kind=C_INT) function cel_free_int_c(array) &
     bind(C, name='cel_free_int_c')
     use, intrinsic :: iso_c_binding
     implicit none
     integer(kind = C_INT) ::array
   end function cel_free_int_c

   integer(kind=C_INT) function setThreadPolicy_sched(thread_num) &
     bind(C, name='setThreadPolicy_sched')
     use, intrinsic :: iso_c_binding
     implicit none
     integer(kind = C_INT) ::thread_num
   end function setThreadPolicy_sched

 end interface

end module cel_timer_interface_module
