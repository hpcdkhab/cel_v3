!     
! File:   cel_omp_module.f90
! Project CRESTA (see details on https://www.cresta-project.eu) Exascale library
! Send you email about the bugs to
! The High Performance Computing Center Stuttgart (HLRS) of the University of Stuttgart
! Dmitry Khabi,  (email: khabi@hlrs.de)
! Created on Jul 5, 2013
!

module cel_omp_module
use cel_types_module
use cel_error_module
use OMP_LIB
implicit none

  integer, parameter :: cel_master_out_thread = 0
  integer, parameter :: cel_worker_out_thread = 0

interface init
  module procedure cel_omp_init
end interface init

interface print
  module procedure cel_omp_print
end interface print
 
  type cel_omp_type
    integer(kind=ik)   :: num_threads
    integer(kind=ik)   :: thread_num
    integer(kind=ik)   :: num_masters
    integer(kind=ik)   :: master_num
    integer(kind=mpik) :: master_comm
    integer(kind=mpik) :: mpi_type_ik
    integer(kind=mpik) :: mpi_type_rk
    integer(kind=mpik) :: mpi_type_perf_ik
    integer(kind=ik)   :: num_workers
    integer(kind=ik)   :: worker_num
    integer(kind=ik)   :: profile_indx
    logical(kind=lk)   :: verbosity
    logical(kind=lk)   :: is_master
    logical(kind=lk)   :: is_worker
  end type cel_omp_type


contains

subroutine cel_omp_init(cel_omp, num_threads, thread_num,&
  mpi_size, mpi_rank, mpi_comm, mpi_type_ik, mpi_type_rk,&
  mpi_type_perf_ik, output_on)
  type(cel_omp_type), intent(inout) :: cel_omp
  integer(kind=ik), intent(in) :: num_threads
  integer(kind=ik), intent(in) :: thread_num
  integer(kind=mpik), intent(in) :: mpi_size
  integer(kind=mpik), intent(in) :: mpi_rank
  integer(kind=mpik), intent(in) :: mpi_comm
  integer(kind=mpik), intent(in) :: mpi_type_rk
  integer(kind=mpik), intent(in) :: mpi_type_ik
  integer(kind=mpik), intent(in) :: mpi_type_perf_ik
  logical(kind=lk), intent(in) :: output_on

  cel_omp%verbosity = output_on
  cel_omp%num_threads = num_threads
  cel_omp%thread_num = thread_num
  cel_omp%mpi_type_rk = mpi_type_rk
  cel_omp%mpi_type_ik = mpi_type_ik
  cel_omp%mpi_type_perf_ik = mpi_type_perf_ik
  !if(cel_omp%num_workers + cel_omp%num_masters .NE. cel_omp%num_threads) then
  !  call cel_error("Error cel_omp_init: check master and worker number",&
  !   1_ik, output_on, .TRUE._lk)
  !end if
  if (cel_omp%thread_num .LT. 1_ik) then
    cel_omp%is_master=.TRUE._lk
    cel_omp%is_worker=.FALSE._lk
    cel_omp%num_workers=num_threads-1
    cel_omp%worker_num=cel_index_undef
    cel_omp%verbosity = .TRUE._lk
  else
    cel_omp%is_master=.FALSE._lk
    cel_omp%is_worker=.TRUE._lk
    cel_omp%num_workers=num_threads-1
    cel_omp%worker_num=thread_num
    cel_omp%verbosity = .FALSE._lk
  endif
  if (cel_omp%thread_num .LT. 1_ik) then
  end if
  cel_omp%profile_indx = thread_num + 1_ik
  cel_omp%num_masters = mpi_size
  cel_omp%master_num = mpi_rank 
  cel_omp%master_comm = mpi_comm
  
end subroutine cel_omp_init


subroutine cel_omp_print(cel_omp, output_on)
  type(cel_omp_type) :: cel_omp
  logical(kind=lk) :: output_on
  
  if(output_on) then
    write(*,'(A,I0,A)') "cel_omp_print(",cel_omp%thread_num,"):"
    write(*,'(A,I0)') "cel_omp%num_threads: ", cel_omp%num_threads
    write(*,'(A,I0)') "cel_omp%thread_num: ", cel_omp%thread_num
    write(*,'(A,I0)') "cel_omp%num_masters: ", cel_omp%num_masters
    write(*,'(A,I0)') "cel_omp%master_num: ", cel_omp%master_num
    write(*,'(A,I0)') "cel_omp%num_workers: ", cel_omp%num_workers
    write(*,'(A,I0)') "cel_omp%worker_num: ", cel_omp%worker_num
    write(*,'(A,I0)') "cel_omp%profile_idx: ", cel_omp%profile_indx
    write(*,'(A,L5)') "cel_omp%verbosity: ", cel_omp%verbosity
    write(*,'(A,L5)') "cel_omp%is_master: ", cel_omp%is_master
    write(*,'(A,L5)') "cel_omp%is_worker: ", cel_omp%is_worker
  end if
end subroutine cel_omp_print

end module cel_omp_module
