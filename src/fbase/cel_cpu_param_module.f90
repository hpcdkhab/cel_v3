!
! File:   cel_cpu_param_module.f90
! Project CRESTA (see details on https://cresta-project.eu) Exascale library
! Send you email about the bugs to
! The High Performance Computing Center Stuttgart (HLRS) of the University of Stuttgart
! Dmitry Khabi,  (email: khabi@hlrs.de)
! Created on May 29, 2013
!

module cel_cpu_param_module
use cel_types_module
implicit none
  character(LEN=1024), parameter :: cel_cpu_param_description="Numtest4: Intel(R) Xeon(R) CPU E5-2687W 0 @ 3.10GHz"
  !character(LEN=1024) :: cel_cpu_param_description="Intel(R) Xeon(R) CPU E5-2670 0 @ 2.60GHz"
  integer(kind=ik), parameter :: cel_cpu_param_cores = 12 ! number of cores on cpu
  integer(kind=ik), parameter :: cel_cpu_param_hwthreads = 24 ! number of hadrware threads
  integer(kind=ik), parameter :: cel_cpu_param_max_num_threads = 128 ! number of hadrware threads
  integer(kind=ik), parameter :: cel_cpu_param_l1_cache_size = (32)*1024 ! L1 Cache
  integer(kind=ik), parameter :: cel_cpu_param_core_cache_size = (256-32)*1024 ! L2 Cache
  integer(kind=ik), parameter :: cel_cpu_param_shared_cache_size = 30*1024*1024-&
   cel_cpu_param_core_cache_size*cel_cpu_param_cores ! L3 Cache
  integer(kind=ik), parameter :: cel_cpu_param_num_sockets = 2 ! number of the sockets on the board
  integer(kind=ik), parameter :: cel_cpu_param_cache_line_bytes = 64 ! number of cores on cpu
  real(kind=rk), parameter ,dimension(16):: cel_cpu_frequencies = &
     (/ 3.11,3.1,3.0,2.8,2.7,2.5,2.4,2.3,2.1,2.0,1.9,1.7,1.6,1.5,1.3,1.2/)
end module cel_cpu_param_module
