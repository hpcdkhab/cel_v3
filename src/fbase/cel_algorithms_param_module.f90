!
! File:   cel_algorithms_param_module.f90
! Project CRESTA (see details on https://cresta-project.eu) Exascale library
! Send you email about the bugs to
! The High Performance Computing Center Stuttgart (HLRS) of the University of Stuttgart
! Dmitry Khabi,  (email: khabi@hlrs.de)
! Created on May 29, 2013
!

module cel_algorithms_param_module
use cel_types_module
implicit none
  integer(kind=ik), parameter :: cel_algorithms_param_merge_max_array_length = 32 ! the array length to switch merge sort to buble sort
  integer(kind=ik), parameter :: cel_algorithms_param_merge_max_r_level = 32 ! maximum recursion level of the merge sort

end module cel_algorithms_param_module
