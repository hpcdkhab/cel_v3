! Definition of the distributed matrix vector multiplication
! File:   cel_fop_mv_module.f90
! Project CRESTA Exascale library (see details on https://www.cresta-project.eu)
! Send you email about the bugs to
! The High Performance Computing Center Stuttgart (HLRS) of the University of Stuttgart
! Dmitry Khabi,  (email: khabi@hlrs.de)
! Created on Jun 14, 2013
!

module cel_fop_mv_module
use cel_types_module
use cel_timer_interface_module
use cel_cpu_param_module
use cel_sp_mat_distr_opt_module
use cel_sp_mat_distr_vec_module
use cel_sp_mat_module
use cel_sp_mat_distr_module
use cel_sp_mat_distr_gl_module
use cel_perf_module
use cel_base_print_module
use cel_omp_module
use cel_omp_shared_work_module
use MPI
use OMP_LIB
implicit none
contains
!Distributed matrix vector multiplication y=y+A*x
subroutine cel_fop_mv_axy_distr(a_sp_mats,x_distr_vectors,y_distr_vectors,&
                                distr_vect_temp,cel_comm,vtxdist,cel_omp,cel_omp_shared_work,&
                                cel_perf_counter,output_on,err_code)
  type(cel_sp_mat_type), dimension(:), allocatable,intent(inout)    :: a_sp_mats
  type(cel_sp_mat_distr_vec_type), intent(inout)                    :: x_distr_vectors
  type(cel_sp_mat_distr_vec_type), intent(inout)                    :: y_distr_vectors
  type(cel_sp_mat_distr_vec_type), intent(inout)                    :: distr_vect_temp
  type(cel_comm_type), intent(inout)                                :: cel_comm
  type(cel_omp_shared_work_type), intent(inout)                     :: cel_omp_shared_work
  integer(kind=ik),dimension(:), intent(in)                         :: vtxdist
  type(cel_omp_type), intent(in)                                    :: cel_omp
  type(cel_perf_counter_type), intent(inout)                        :: cel_perf_counter
  logical(kind=lk)  , intent(in)                                    :: output_on
  integer(kind=ik), intent(out)                                     :: err_code
  integer(kind=ik) :: ii
  integer(kind=mpik) :: ierr
  real(kind=trk) :: time, time_hallo, time_local
  err_code = 0_ik
  !call reset(cel_perf_counter)
  !start communication
  if(cel_omp%is_master) time = get_time()
  if(cel_omp%is_master .and. cel_omp%num_masters .GT. 1_ik) then
    call cel_fop_mv_axy_distr_master_comm(x_distr_vectors,cel_comm,&
                                          vtxdist,a_sp_mats(1_ik)%mv_communicator,&
                                          cel_omp,cel_perf_counter,&
                                          output_on,err_code)
    cel_perf_counter%time_comm(1_ik) = cel_perf_counter%time_comm(1) +  (get_time() - time)
  elseif(cel_omp%is_master) then
    cel_perf_counter%time_comm(1_ik) = cel_perf_counter%time_comm(1) +  (get_time() - time)
  end if
  if(cel_omp%is_worker) then
    if(cel_omp%worker_num .eq. cel_omp%num_workers) time_local = get_time()
    call cel_fop_mv_axy_distr_worker_local(a_sp_mats(1_ik),x_distr_vectors%values(:,1_ik),&
                                           y_distr_vectors%values(:,1_ik),&
                                           cel_omp,cel_omp_shared_work,&
                                           cel_perf_counter,output_on,err_code)
    if(cel_omp%worker_num .eq. cel_omp%num_workers) then
      cel_perf_counter%time_local_calc(1_ik) = cel_perf_counter%time_local_calc(1_ik) +  (get_time() - time_local)
    end if
  end if
  !$omp barrier
  !calculation on received data
  if(cel_omp%is_worker) then
    if(cel_omp%worker_num .eq. cel_omp%num_workers)  time_hallo = get_time()
    do ii=2_ik, size(a_sp_mats)
      call cel_fop_mv_axy_distr_worker_local(a_sp_mats(ii),x_distr_vectors%values(:,ii),&
                                          y_distr_vectors%values(:,1_ik),&
                                          cel_omp,cel_omp_shared_work,&
                                          cel_perf_counter,output_on,err_code)
    end do
    if(cel_omp%worker_num .eq. cel_omp%num_workers) then
      cel_perf_counter%time_halo_calc(1_ik) = cel_perf_counter%time_halo_calc(1) +  (get_time() - time_hallo)
    end if
  end if
  !$omp barrier
  if(cel_omp%is_master) then
    cel_perf_counter%time(1_ik) = cel_perf_counter%time(1_ik) + (get_time() - time)
  end if
  
end subroutine cel_fop_mv_axy_distr

!Communication part of distributed matrix vector multiplication
subroutine cel_fop_mv_axy_distr_master_comm(x_distr_vectors,cel_comm,vtxdist,mv_communicator,&
                                            cel_omp, cel_perf_counter,output_on,err_code)
  type(cel_sp_mat_distr_vec_type), intent(inout)                    :: x_distr_vectors
  type(cel_comm_type), intent(inout)                                :: cel_comm
  integer(kind=ik),dimension(:), intent(in)                         :: vtxdist
  integer(kind=ik), intent(in)                                      :: mv_communicator
  type(cel_omp_type), intent(in)                                    :: cel_omp
  type(cel_perf_counter_type), intent(inout)                        :: cel_perf_counter
  logical(kind=lk)  , intent(in)                                    :: output_on
  integer(kind=ik), intent(out)                                     :: err_code

  err_code = 0_ik
  select case (mv_communicator)
    case(MV_COMM_ISEND)
      call cel_sp_mat_distr_vec_isend_start_all(x_distr_vectors,cel_comm,vtxdist,cel_omp,&
                                      cel_perf_counter,err_code)
      !call cel_error("Error cel_fop_mv_axy_distr_master_comm: vec_start_all", &
      !               err_code,.TRUE._lk,.FALSE._lk)
      call cel_sp_mat_distr_vec_end_all(cel_comm,cel_omp,err_code)
            !call cel_error("Error cel_fop_mv_axy_distr_master_comm: vec_end_all", &
      !           err_code,.TRUE._lk,.FALSE._lk)
    case(MV_COMM_ISEND_IDX)
      call cel_sp_mat_distr_vec_isend_idx_start_all(x_distr_vectors,cel_comm,vtxdist,cel_omp,&
                                      cel_perf_counter,err_code)
      !call cel_error("Error cel_fop_mv_axy_distr_master_comm: vec_start_all", &
      !               err_code,.TRUE._lk,.FALSE._lk)
      call cel_sp_mat_distr_vec_end_all(cel_comm,cel_omp,err_code)
            !call cel_error("Error cel_fop_mv_axy_distr_master_comm: vec_end_all", &
      !           err_code,.TRUE._lk,.FALSE._lk)
    case(MV_COMM_ISEND_GROUPS)
      call cel_sp_mat_distr_vec_isend_groups_start_all(x_distr_vectors,cel_comm,vtxdist,cel_omp,&
                                      cel_perf_counter,err_code)
      !call cel_error("Error cel_fop_mv_axy_distr_master_comm: vec_start_all", &
      !               err_code,.TRUE._lk,.FALSE._lk)
      call cel_sp_mat_distr_vec_end_all(cel_comm,cel_omp,err_code)
            !call cel_error("Error cel_fop_mv_axy_distr_master_comm: vec_end_all", &
      !           err_code,.TRUE._lk,.FALSE._lk)
    case(MV_COMM_IBCAST)
      call cel_sp_mat_distr_vec_ibcast_start_all(x_distr_vectors,cel_comm,vtxdist,cel_omp,&
                                          cel_perf_counter,err_code)
      !call cel_error("Error cel_fop_mv_axy_distr_master_comm: ibcast_start_all", &
      !               err_code,.TRUE._lk,.FALSE._lk)
      call cel_sp_mat_distr_vec_ibcast_end_all(cel_comm,cel_omp,err_code)
      !call cel_error("Error cel_fop_mv_axy_distr_master_comm: cel_sp_mat_distr_vec_ibcast_end_all", &
      !               err_code,.TRUE._lk,.FALSE._lk)
    case(MV_COMM_ISEND_COMM_BUFFER)
      call cel_sp_mat_distr_vec_isend_comm_buffer_start_all(x_distr_vectors,cel_comm,&
         vtxdist,cel_omp,cel_perf_counter,err_code)
      !call cel_error("Error cel_fop_mv_axy_distr_master_comm: ibcast_start_all", &
      !               err_code,.TRUE._lk,.FALSE._lk)
      call cel_sp_mat_distr_vec_end_all(cel_comm,cel_omp,err_code)
            !call cel_error("Error cel_fop_mv_axy_distr_master_comm: vec_end_all", &
      !           err_code,.TRUE._lk,.FALSE._lk)
      call cel_sp_mat_distr_vec_synchronize_comm_buffer_master(x_distr_vectors,cel_comm,&
             cel_omp,cel_perf_counter,err_code)
    case default
        call cel_error("Error cel_fop_mv_axy_distr_worker_local: no mv communicator of the choised type defined", &
                   err_code,.TRUE._lk,.TRUE._lk)
    end select

end subroutine cel_fop_mv_axy_distr_master_comm

!Calculation part of distributed matrix vector multiplication
subroutine cel_fop_mv_axy_distr_worker_local(a_sp_mat,x_local_vector,y_local_vector,&
                                             cel_omp,cel_omp_shared_work,&
                                             cel_perf_counter,output_on,err_code)
  type(cel_sp_mat_type) ,intent(inout)                            :: a_sp_mat
  real(kind=rk), dimension(:), intent(in)                           :: x_local_vector
  real(kind=rk), dimension(:), intent(inout)                        :: y_local_vector
  type(cel_omp_type), intent(in)                                    :: cel_omp
  type(cel_omp_shared_work_type), intent(inout)                     :: cel_omp_shared_work
  type(cel_perf_counter_type), intent(inout)                            :: cel_perf_counter
  logical(kind=lk)  , intent(in)                                    :: output_on
  integer(kind=ik), intent(out)                                     :: err_code
  err_code = 0_ik

  select case (a_sp_mat%mv_algorithm) 
    case(MV_CSR)
      call cel_fop_mv_axy_distr_worker_local_csr(a_sp_mat,x_local_vector,y_local_vector,&
                                             cel_omp,cel_omp_shared_work,&
                                             cel_perf_counter,output_on,err_code)
    case(MV_JAD)
      call cel_fop_mv_axy_distr_worker_local_jad(a_sp_mat,x_local_vector,y_local_vector,&
                                             cel_omp,cel_omp_shared_work,&
                                             cel_perf_counter,output_on,err_code)
    case(MV_COO)
      call cel_fop_mv_axy_distr_worker_local_coo(a_sp_mat,x_local_vector,y_local_vector,&
                                             cel_omp,cel_omp_shared_work,&
                                             cel_perf_counter,output_on,err_code)
    case default
      call cel_error("Error cel_fop_mv_axy_distr_worker_local: no mv algorithm defined", &
                 err_code,.TRUE._lk,.TRUE._lk)
  end select
  call cel_error("Error cel_fop_mv_axy_distr_worker_local: mv", &
                 err_code,cel_is_in_debug,.FALSE._lk)

end subroutine cel_fop_mv_axy_distr_worker_local

!Calculation part of distributed matrix vector multiplication in format csr
subroutine cel_fop_mv_axy_distr_worker_local_csr(a_sp_mat,x_local_vector,y_local_vector,&
                                                 cel_omp,cel_omp_shared_work,&
                                                 cel_perf_counter, output_on,err_code)
  type(cel_sp_mat_type) ,intent(inout)                            :: a_sp_mat
  real(kind=rk),dimension(:), intent(in)                            :: x_local_vector
  real(kind=rk),dimension(:), intent(inout)                         :: y_local_vector
  type(cel_omp_type), intent(in)                                    :: cel_omp
  type(cel_omp_shared_work_type), intent(inout)                     :: cel_omp_shared_work
  type(cel_perf_counter_type), intent(inout)                        :: cel_perf_counter
  logical(kind=lk)  , intent(in)                                    :: output_on
  integer(kind=ik), intent(out)                                     :: err_code
  integer(kind=ik) :: ii, jj, nrow, col_idx
  integer(kind=ik) :: start_idx,end_idx,omp_start_idx,omp_end_idx,num_threads
  integer(kind=ik) :: total_number,col_offset
  real(kind=rk) ::value_tmp, value, value_x

  err_code = 0_ik
  
  nrow = size(a_sp_mat%xadj) - 1_ik
  call cel_omp_shared_work_distr(cel_omp_shared_work,&
                                 omp_start_idx,omp_end_idx,num_threads,&
                                 nrow,a_sp_mat%elements_num,&
                                 cel_omp,output_on,err_code)
  if(cel_omp%worker_num .le. num_threads) then
    total_number = 0_perf_ik
    col_offset = a_sp_mat%col_offset
    do ii=omp_start_idx, omp_end_idx
      value_tmp = 0.0_rk
      start_idx = a_sp_mat%xadj(ii)
      end_idx = a_sp_mat%xadj(ii+1_ik) - 1_ik
      total_number = total_number + (end_idx-start_idx+1_ik)
      !dir$ ivdep
      !dir$ vector always
      do jj = start_idx, end_idx
        col_idx = a_sp_mat%cols(jj) - col_offset
        value = a_sp_mat%values(jj)
        value_x = x_local_vector(col_idx)
        value_tmp = value_tmp + value * value_x
      end do
      y_local_vector(ii) =y_local_vector(ii) + value_tmp
    end do
    if(cel_omp%worker_num .eq. 1_ik) then
      cel_perf_counter%add = cel_perf_counter%add+&
        a_sp_mat%elements_num+nrow
      cel_perf_counter%mult = cel_perf_counter%mult+&
        a_sp_mat%elements_num 
      cel_perf_counter%load = cel_perf_counter%load+&
        nrow*2_ik*ik+a_sp_mat%elements_num*1_ik*ik+a_sp_mat%elements_num*2_ik*rk
      cel_perf_counter%store = cel_perf_counter%store+&
       nrow*rk
     end if
  else
    !free worker thread
    !write(*,'(I0,A,I0)') cel_omp%master_num,": worker_num:",  cel_omp%worker_num,"; num_threads:", num_threads
    !call cel_error("Error cel_fop_mv_axy_distr_worker_local_csr:  not enought threads", &
    !             1_ik,.TRUE._lk,.FALSE._lk)
  end if

end subroutine cel_fop_mv_axy_distr_worker_local_csr

!Calculation part of distributed matrix vector multiplication in format csr
subroutine cel_fop_mv_axy_distr_worker_local_jad(a_sp_mat,x_local_vector,y_local_vector,&
                                                 cel_omp,cel_omp_shared_work,&
                                                 cel_perf_counter, output_on,err_code)
  type(cel_sp_mat_type) ,intent(inout)                              :: a_sp_mat
  real(kind=rk),dimension(:), intent(in)                            :: x_local_vector
  real(kind=rk),dimension(:), intent(inout)                         :: y_local_vector
  type(cel_omp_type), intent(in)                                    :: cel_omp
  type(cel_omp_shared_work_type), intent(inout)                     :: cel_omp_shared_work
  type(cel_perf_counter_type), intent(inout)                        :: cel_perf_counter
  logical(kind=lk)  , intent(in)                                    :: output_on
  integer(kind=ik), intent(out)                                     :: err_code
  integer(kind=ik) :: ii, jj, length,k1, length_max
  real(kind=rk)    :: real_1, time, time_2, time_3
  integer(kind=ik) :: el_index_1
  integer(kind=ik) :: omp_start_idx, omp_end_idx, num_threads
  integer(kind=ik) :: col_offset, row_offset

  err_code = 0_ik
  length_max = size(a_sp_mat%jad_tmp_values,dim=1)
  if(cel_omp%worker_num .eq. cel_omp%num_workers) time = get_time()
  !reset the temporaly vector
  call cel_omp_shared_work_start(cel_omp_shared_work,cel_omp%num_workers)
  call cel_omp_shared_work_distr(cel_omp_shared_work,&
                                 omp_start_idx,omp_end_idx,num_threads,&
                                 length_max,length_max*VRL,&
                                 cel_omp,output_on,err_code)
  if(cel_omp%worker_num .le. num_threads) then
    a_sp_mat%jad_tmp_values(omp_start_idx:omp_end_idx) = 0.0_rk
  end if
  call cel_omp_shared_work_end_and_wait(cel_omp_shared_work,cel_omp,LOOP_SLEEP_DEFLENGTH)
  !write(cel_output_unit,'(2(A,I0))') "lost work_end:", cel_omp%worker_num,&
  !  "; length_max:",length_max
  if(cel_omp%worker_num .eq. cel_omp%num_workers) then
      cel_perf_counter%time_jad_mv_phase_1(1) = cel_perf_counter%time_jad_mv_phase_1(1)+&
        (get_time() - time)
      time_2 = get_time()
  end if
  !perform matrix vector multiplication
  call cel_omp_shared_work_start(cel_omp_shared_work,cel_omp%num_workers)
  call cel_omp_shared_work_distr(cel_omp_shared_work,&
                               omp_start_idx,omp_end_idx,num_threads,&
                               a_sp_mat%jad_mmax,a_sp_mat%jad_mmax*VRL,&
                               cel_omp,output_on,err_code)
  if(cel_omp%worker_num .le. num_threads) then
    do ii=omp_start_idx, omp_end_idx
      k1=a_sp_mat%jad_begin(ii)-1_ik
      length = a_sp_mat%jad_begin(ii+1)-k1 - 1_ik
      do jj=1_ik, length
        el_index_1 = a_sp_mat%jad_cols(k1+jj)
        real_1 = a_sp_mat%jad_values(k1+jj)*x_local_vector(el_index_1)
        a_sp_mat%jad_tmp_values(jj)=&
          a_sp_mat%jad_tmp_values(jj)+real_1
      end do
    end do
  end if

  call cel_omp_shared_work_end_and_wait(cel_omp_shared_work,cel_omp,LOOP_SLEEP_DEFLENGTH)
!  write(cel_output_unit,'(4(A,I0))') "end:", cel_omp%worker_num,&
!    "; from:",cel_omp%num_workers,"; omp_start_idx:", omp_start_idx, ": omp_end_idx:",omp_end_idx
  if(cel_omp%worker_num .eq. cel_omp%num_workers) then
      cel_perf_counter%time_jad_mv_phase_2(1) = cel_perf_counter%time_jad_mv_phase_2(1)+&
        (get_time() - time_2)
      time_3 = get_time()
  end if
  !copy and permute temporaly values to the solution vector
  call cel_omp_shared_work_start(cel_omp_shared_work,cel_omp%num_workers)
  call cel_omp_shared_work_distr(cel_omp_shared_work,&
                                 omp_start_idx,omp_end_idx,num_threads,&
                                 a_sp_mat%jad_nmax,a_sp_mat%jad_nmax,&
                                 cel_omp,output_on,err_code)
 ! write(cel_output_unit,'(4(A,I0))') "3 shared start:", cel_omp%thread_num,&
 !   "; from:",cel_omp%num_workers,"; omp_start_idx:", omp_start_idx, ": omp_end_idx:",omp_end_idx
  if(cel_omp%worker_num .le. num_threads) then
    do ii=omp_start_idx, omp_end_idx
      el_index_1 = a_sp_mat%jad_perm(ii)
      y_local_vector(el_index_1) = a_sp_mat%jad_tmp_values(ii)
    end do
  end if
 ! write(cel_output_unit,'(4(A,I0))') "3 shared end:", cel_omp%thread_num,&
  !  "; from:",cel_omp%num_workers,"; omp_start_idx:", omp_start_idx, ": omp_end_idx:",omp_end_idx
  call cel_omp_shared_work_end_and_wait(cel_omp_shared_work,cel_omp,LOOP_SLEEP_DEFLENGTH)
  !write(cel_output_unit,'(4(A,I0))') "go 3 shared end:", cel_omp%thread_num,&
 !   "; from:",cel_omp%num_workers,"; omp_start_idx:", omp_start_idx, ": omp_end_idx:",omp_end_idx
  if(cel_omp%worker_num .eq. cel_omp%num_workers) then
      cel_perf_counter%time_jad_mv_phase_3(1) = cel_perf_counter%time_jad_mv_phase_3(1)+&
        (get_time() - time_3)
  end if
!    write(cel_output_unit,'(4(A,I0))') "work_end 2:", cel_omp%worker_num,&
!    "; from:",cel_omp%num_workers,"; omp_start_idx:", omp_start_idx, ": omp_end_idx:",omp_end_idx
  if(cel_omp%worker_num .eq. 1_ik) then
    cel_perf_counter%add = cel_perf_counter%add + a_sp_mat%elements_num +&
      length_max*num_threads
    cel_perf_counter%mult = cel_perf_counter%mult + a_sp_mat%elements_num
    cel_perf_counter%load = cel_perf_counter%load +&
      a_sp_mat%jad_mmax*2_ik*ik + a_sp_mat%elements_num*1_ik*ik + &
      a_sp_mat%elements_num*3_ik*rk+length_max*num_threads*2_ik*rk + &
      a_sp_mat%jad_nmax*ik+a_sp_mat%jad_nmax*rk
    cel_perf_counter%store = cel_perf_counter%store + &
      length_max*rk + a_sp_mat%elements_num*rk + &
      length_max*rk + a_sp_mat%jad_nmax*ik+a_sp_mat%jad_nmax*rk
  end if

end subroutine cel_fop_mv_axy_distr_worker_local_jad



!Calculation part of distributed matrix vector multiplication in format coo
subroutine cel_fop_mv_axy_distr_worker_local_coo(a_sp_mat,x_local_vector,y_local_vector,&
                                                 cel_omp,cel_omp_shared_work,&
                                                 cel_perf_counter, output_on,err_code)
  type(cel_sp_mat_type) ,intent(inout)                            :: a_sp_mat
  real(kind=rk),dimension(:), intent(in)                            :: x_local_vector
  real(kind=rk),dimension(:), intent(inout)                         :: y_local_vector
  type(cel_omp_type), intent(in)                                    :: cel_omp
  type(cel_omp_shared_work_type), intent(inout)                     :: cel_omp_shared_work
  type(cel_perf_counter_type), intent(inout)                            :: cel_perf_counter
  logical(kind=lk)  , intent(in)                                    :: output_on
  integer(kind=ik), intent(out)                                     :: err_code
  integer(kind=ik) :: ii, col_idx
  integer(kind=ik) :: omp_start_idx,omp_end_idx,num_threads
  integer(kind=ik) :: total_number, col_offset, row_offset, row_idx
  real(kind=rk) :: value, value_tmp

  err_code = 0_ik

  call cel_omp_shared_work_distr(cel_omp_shared_work,&
                                 omp_start_idx,omp_end_idx,num_threads,&
                                 a_sp_mat%elements_num,a_sp_mat%elements_num,&
                                 cel_omp,output_on,err_code)
  if(cel_omp%worker_num .le. num_threads) then
    total_number = (omp_end_idx-omp_start_idx+1_ik)
    col_offset =  a_sp_mat%col_offset
    row_offset =  a_sp_mat%row_offset
   ! write(*,'(I0,A,I0,A,I0,A,I0, A,I0,A,I0)') cel_omp%worker_num, ": elems:",a_sp_mat%elements_num,&
   !   "; omp_start_idx:",omp_start_idx,"; omp_end_idx:",omp_end_idx,&
   !   "; col_offset:",col_offset,"; row_offset:",row_offset 
    do ii = omp_start_idx, omp_end_idx
      col_idx = a_sp_mat%cols(ii) - col_offset
      row_idx = a_sp_mat%rows(ii) - row_offset
      value = a_sp_mat%values(ii)
      value_tmp = value * x_local_vector(col_idx)
      !$omp atomic update
      y_local_vector(row_idx) = y_local_vector(row_idx) + value_tmp
    end do
    if(cel_omp%worker_num .eq. 1_ik) then
      cel_perf_counter%add = cel_perf_counter%add + a_sp_mat%elements_num
      cel_perf_counter%mult = cel_perf_counter%add + a_sp_mat%elements_num
      cel_perf_counter%load = cel_perf_counter%load +a_sp_mat%elements_num*2_ik*ik+&
         a_sp_mat%elements_num*2_ik*rk
      cel_perf_counter%store = cel_perf_counter%mult + a_sp_mat%elements_num*rk
    end if
  end if

end subroutine cel_fop_mv_axy_distr_worker_local_coo

end module cel_fop_mv_module

