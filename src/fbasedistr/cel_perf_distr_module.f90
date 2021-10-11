!     
! File:   cel_perf_distr_module.f90
! Project CRESTA (see details on https://www.cresta-project.eu) Exascale library
! Send you email about the bugs to
! The High Performance Computing Center Stuttgart (HLRS) of the University of Stuttgart
! Dmitry Khabi,  (email: khabi@hlrs.de)
! Created on Dezember 03, 2013
!

module cel_perf_distr_module
use cel_types_module
use cel_timer_interface_module
use cel_error_module
use cel_base_print_module
use cel_omp_module
use cel_perf_module
use OMP_LIB
use MPI
implicit none

contains

subroutine cel_perf_distr_collect(cel_perf_total,cel_perf_local,cel_omp,collect_by_all)
  type(cel_perf_counter_type), intent(inout)  :: cel_perf_total
  type(cel_perf_counter_type), intent(in)     :: cel_perf_local
  type(cel_omp_type), intent(in)              :: cel_omp
  logical(kind=lk), intent(in) :: collect_by_all
  integer(kind=perf_ik),dimension(15) :: int_buffer
  integer(kind=perf_ik),dimension(15) :: int_buffer_recv
  real(kind=rk),dimension(7) :: real_buffer
  real(kind=rk),dimension(7,5) :: real_buffer_recv
  integer(kind=mpik) :: int_count
  integer(kind=mpik) :: real_count
  integer(kind=mpik) :: ierr
  
  if(cel_omp%is_master) then
    int_buffer(1_ik)=int(cel_perf_local%add,kind=perf_ik)
    int_buffer(2_ik)=int(cel_perf_local%mult,kind=perf_ik)
    int_buffer(3_ik)=int(cel_perf_local%divide,kind=perf_ik)
    int_buffer(4_ik)=int(cel_perf_local%sqr,kind=perf_ik)
    int_buffer(5_ik)=int(cel_perf_local%abso,kind=perf_ik)
    int_buffer(6_ik)=int(cel_perf_local%load,kind=perf_ik)
    int_buffer(7_ik)=int(cel_perf_local%store,kind=perf_ik)
    int_buffer(8_ik)=int(cel_perf_local%omp_locks,kind=perf_ik)
    int_buffer(9_ik)=int(cel_perf_local%number_of_repetitions,kind=perf_ik)
    int_buffer(10_ik)=int(cel_perf_local%size,kind=perf_ik)
    int_buffer(11_ik)=int(cel_perf_local%send_bytes,kind=perf_ik)
    int_buffer(12_ik)=int(cel_perf_local%get_bytes,kind=perf_ik)
    int_buffer(13_ik)=int(cel_perf_local%send_messages_num,kind=perf_ik)
    int_buffer(14_ik)=int(cel_perf_local%get_messages_num,kind=perf_ik)
    int_buffer(15_ik)=int(cel_perf_local%num_elements,kind=perf_ik)
    int_count = 15_mpik
    real_buffer(1_ik)=real(cel_perf_local%time(1_ik),kind=rk)
    real_buffer(2_ik)=real(cel_perf_local%time_comm(1_ik),kind=rk)
    real_buffer(3_ik)=real(cel_perf_local%time_local_calc(1_ik),kind=rk)
    real_buffer(4_ik)=real(cel_perf_local%time_halo_calc(1_ik),kind=rk)
    real_buffer(5_ik)=real(cel_perf_local%time_jad_mv_phase_1(1_ik),kind=rk)
    real_buffer(6_ik)=real(cel_perf_local%time_jad_mv_phase_2(1_ik),kind=rk)
    real_buffer(7_ik)=real(cel_perf_local%time_jad_mv_phase_3(1_ik),kind=rk)
    
    real_count= 7_mpik
    if(collect_by_all) then
      call MPI_Allreduce(int_buffer(:),int_buffer_recv(:),int_count,&
                        cel_omp%mpi_type_perf_ik,MPI_SUM,cel_omp%master_comm,ierr)
      call MPI_Allreduce(real_buffer(:),real_buffer_recv(:,1),real_count,&
                        cel_omp%mpi_type_rk,MPI_SUM,cel_omp%master_comm,ierr)
      call MPI_Allreduce(real_buffer(:),real_buffer_recv(:,2),real_count,&
                        cel_omp%mpi_type_rk,MPI_MIN,cel_omp%master_comm,ierr)
      call MPI_Allreduce(real_buffer(:),real_buffer_recv(:,3),real_count,&
                        cel_omp%mpi_type_rk,MPI_MAX,cel_omp%master_comm,ierr)

      real_buffer_recv(:,4) = real_buffer_recv(:,1) / real(cel_omp%num_masters,kind=8)

    end if
    cel_perf_total%add = cel_perf_total%add + int(int_buffer_recv(1_ik),kind=perf_ik)
    cel_perf_total%mult = cel_perf_total%mult + int(int_buffer_recv(2_ik),kind=perf_ik)
    cel_perf_total%divide = cel_perf_total%divide + int(int_buffer_recv(3_ik),kind=perf_ik)
    cel_perf_total%sqr = cel_perf_total%sqr + int(int_buffer_recv(4_ik),kind=perf_ik)
    cel_perf_total%abso = cel_perf_total%abso + int(int_buffer_recv(5_ik),kind=perf_ik)
    cel_perf_total%load = cel_perf_total%load + int(int_buffer_recv(6_ik),kind=perf_ik)
    cel_perf_total%store = cel_perf_total%store + int(int_buffer_recv(7_ik),kind=perf_ik)
    cel_perf_total%omp_locks = cel_perf_total%omp_locks + int(int_buffer_recv(8_ik),kind=perf_ik)
    cel_perf_total%number_of_repetitions = cel_perf_total%number_of_repetitions +&
     int(int_buffer_recv(9_ik),kind=perf_ik)
    cel_perf_total%size = cel_perf_total%size + int(int_buffer_recv(10_ik),kind=perf_ik)
    cel_perf_total%send_bytes = cel_perf_total%send_bytes +&
     int(int_buffer_recv(11_ik),kind=perf_ik)
    cel_perf_total%get_bytes = cel_perf_total%get_bytes +&
     int(int_buffer_recv(12_ik),kind=perf_ik)
    cel_perf_total%send_messages_num = cel_perf_total%send_messages_num +&
     int(int_buffer_recv(13_ik),kind=perf_ik)
    cel_perf_total%get_messages_num = cel_perf_total%get_messages_num +&
     int(int_buffer_recv(14_ik),kind=perf_ik)
    cel_perf_total%num_elements = cel_perf_total%num_elements +&
     int(int_buffer_recv(15_ik),kind=perf_ik)
    cel_perf_total%time(1_ik) =&
     cel_perf_total%time(1_ik) +&
     real(cel_perf_local%time(1_ik),kind=8)
    cel_perf_total%time_comm(1_ik) =&
     cel_perf_total%time_comm(1_ik) + &
      real(cel_perf_local%time_comm(1_ik),kind=8)
    cel_perf_total%time_local_calc(1_ik) = &
     cel_perf_total%time_local_calc(1_ik) +&
      real(cel_perf_local%time_local_calc(1_ik),kind=8)
    cel_perf_total%time_jad_mv_phase_1(1_ik) = cel_perf_total%time_jad_mv_phase_1(1_ik) +&
     real(cel_perf_local%time_jad_mv_phase_1(1_ik),kind=8)
    cel_perf_total%time_jad_mv_phase_2(1_ik) = &
    cel_perf_total%time_jad_mv_phase_1(1_ik) +&
     real(cel_perf_local%time_jad_mv_phase_2(1_ik),kind=8)
    cel_perf_total%time_jad_mv_phase_3(1_ik) =&
     cel_perf_total%time_jad_mv_phase_1(1_ik) +&
     real(cel_perf_local%time_jad_mv_phase_3(1_ik),kind=8)

    cel_perf_total%time(2_ik:5_ik) = &
    cel_perf_total%time(2_ik:5_ik) +&
     real(real_buffer_recv(1_ik,1:4),kind=8)
    cel_perf_total%time_comm(2_ik:5_ik) =&
     cel_perf_total%time_comm (2_ik:5_ik)+&
     real(real_buffer_recv(2_ik,1:4),kind=8)
    cel_perf_total%time_local_calc(2_ik:5_ik) =&
     cel_perf_total%time_local_calc(2_ik:5_ik) +&
     real(real_buffer_recv(3_ik,1:4),kind=8)
    cel_perf_total%time_halo_calc(2_ik:5_ik) =&
     cel_perf_total%time_halo_calc(2_ik:5_ik)+&
     real(real_buffer_recv(4_ik,1:4),kind=8)
    cel_perf_total%time_jad_mv_phase_1(2_ik:5_ik) =&
     cel_perf_total%time_jad_mv_phase_1(2_ik:5_ik)+&
     real(real_buffer_recv(5_ik,1:4),kind=8)
    cel_perf_total%time_jad_mv_phase_2(2_ik:5_ik) =&
     cel_perf_total%time_jad_mv_phase_2(2_ik:5_ik)+&
     real(real_buffer_recv(6_ik,1:4),kind=8)
    cel_perf_total%time_jad_mv_phase_3(2_ik:5_ik) =&
     cel_perf_total%time_jad_mv_phase_3(2_ik:5_ik)+&
     real(real_buffer_recv(7_ik,1:4),kind=8)
    

  end if
end subroutine cel_perf_distr_collect

end module cel_perf_distr_module
