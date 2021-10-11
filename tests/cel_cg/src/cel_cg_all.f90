
program cel_cg_program
!library fbase
use cel_types_module
use cel_input_parameters_module
use cel_omp_module
use cel_timer_interface_module
use cel_perf_module
use cel_perf_cg_module
use cel_cpu_param_module
!library fspmt
use cel_sp_mat_generate_module
use cel_cg_module

!external libraries
use OMP_LIB
use MPI
use, intrinsic :: iso_c_binding

implicit none
  type(cel_input_parameters_type)               :: input_parameters !user input
  integer(kind = mpik)                          :: ierr, mpi_comm, mpi_rank !mpi parameters
  integer(kind = mpik)                          :: mpi_size, mpi_type_ik, ik_mpi !mpi parameters
  integer(kind = mpik)                          :: mpi_type_rk, rk_mpi !mpi parameters
  logical(kind=lk)                              :: output_on !output parameters
  logical(kind=lk)                              :: to_stop !debuger parameters
  integer(kind=ik)                              :: test_idx !perf. evaluation
  integer (kind=ik)                             :: curr_num_threads !threads number
  type(cel_omp_type), dimension(:), allocatable :: omp !parameters for threads; thread num, master or worker thread and so on
  integer(kind=ik)                              :: num_threads, omp_indx !OpenMP / threads control
  integer(kind=ik)                              :: block_indx !block_id of coo matrix (distibution parameter)
  type(cel_vec_type)                            :: boundary_vector !boundary vector of system to solve
  type(cel_vec_type)                            :: th_solution_vector !theoretical solution of the system
  type(cel_vec_type)                            :: solution_vector !solution of the system
  real(kind=rk)                                 :: time_start, time_end !perf. evaluation
  integer(kind=ik)                              :: num_rows, num_elements, global_num_elements, global_num_rows !matrix parameters
  integer(kind=ik)                              :: err_code !error variable
  integer(kind=c_int)                           :: num_threads_c !matrix properties
  integer(kind=ik)                              :: ii, num_nodes
  logical(kind=lk)                              :: flag_error
  type(cel_sp_mat_type)                       :: sparse_matrix_coo !sparse matrix in coordinate format to solve 
  type(cel_sp_mat_type)                       :: sparse_matrix_coo_global!filled by proc 0 (before the distribution of the matrix)
  real(kind=rk)                                 :: max_time, average_time, min_time, tmp_time !perf. evaluation
  real(kind=rk)                                 :: mv_max_time, mv_average_time, mv_min_time, mv_tmp_time !perf. evaluation
  real(kind=rk)                                 :: mv_mflops, mv_gb_per_sec, gl_mv_mflops, gl_mv_gb_per_sec
  real(kind=rk)                                 :: mv_operations
  character(len=256)                            :: perf_filename
  integer(kind=ik)                              :: socket_id !0 or 2 in two sockets system
  integer(kind=ik)                              :: num_procs_per_socket
  integer (kind=kmp_affinity_mask_kind)         :: thread_mask
  type(cel_cg_type)                             :: cel_cg
  integer(kind=ik), dimension(5)                :: buffer


  call cel_input_parameters_read(input_parameters)!read user input parameters from command line

  !obtain and set number of threads
  curr_num_threads = input_parameters%omp%num_threads
  num_threads = INT(curr_num_threads,kind=c_ik)  
  num_threads_c = INT(num_threads,kind=c_ik)  
  to_stop = .TRUE._lk
  
  call omp_set_dynamic(.false._lk)
  call omp_set_nested(.false._lk)
  call omp_set_num_threads(num_threads_c)
  
  mpi_comm = MPI_COMM_WORLD
  ik_mpi = INT(ik, kind = mpik)
  rk_mpi = INT(rk, kind = mpik)
  call MPI_INIT(ierr)
  call cel_error_proc("MPI_INIT", mpi_rank, ierr, output_on, to_stop)
  call MPI_COMM_RANK(mpi_comm, mpi_rank, ierr)
  call cel_error_proc("MPI_COMM_RANK", mpi_rank, ierr, output_on, to_stop)
  call MPI_COMM_SIZE(mpi_comm, mpi_size, ierr)
  call cel_error_proc("MPI_COMM_SIZE", mpi_rank, ierr, output_on, to_stop)
  call MPI_TYPE_MATCH_SIZE(MPI_TYPECLASS_INTEGER, ik_mpi, mpi_type_ik, ierr)
  call cel_error_proc("MPI_TYPE_MATCH_SIZE INTEGER", mpi_rank, ierr, output_on, to_stop)
  call MPI_TYPE_MATCH_SIZE(MPI_TYPECLASS_REAL, rk_mpi, mpi_type_rk, ierr)
  call cel_error_proc("MPI_TYPE_MATCH_SIZE REAL", mpi_rank, ierr, output_on, to_stop)
  err_code = craypat_off()
  if (mpi_rank .EQ. 0) then
    cel_is_in_debug = .TRUE._lk
  else
    cel_is_in_debug = .FALSE._lk
  endif

  if (mpi_rank .EQ. 0) then
    output_on = input_parameters%verbosity
    if(input_parameters%verbosity)  cel_cg%verbosity = 1
  else
    output_on = .FALSE._lk
  endif
  !some checks of the input parameters
  if(input_parameters%generate_matrix .GT. 1) then
    if(.not. input_parameters%dx*input_parameters%dy*input_parameters%dz .EQ. mpi_size) then
      if (mpi_rank .EQ. 0) then
        write(*,'(A,A)') &
          "Error: please, check the input parameters:",& 
          " dx*dy*dz must be equal to the number of processes !"
      end if
      CALL MPI_FINALIZE(ierr)
      stop
    end if
    num_nodes = input_parameters%nx*input_parameters%ny*input_parameters%nz/&
    input_parameters%dx*input_parameters%dy*input_parameters%dz
  end if
  if(input_parameters%omp%num_threads .LT. 2) then
    if (mpi_rank .EQ. 0) then
      write(*,'(A,A)') &
        "Error: please, check the input parameters:",&
        " the parameter num_threads must be greater than one!"
    end if
    CALL MPI_FINALIZE(ierr)
    stop
  end if
  
  !allocate memory for thread and processes topology omp
  if(allocated(omp)) deallocate(omp)
  allocate(omp(num_threads))
  
  num_procs_per_socket = 1!cel_cpu_param_num_sockets
  socket_id = mod(mpi_rank,num_procs_per_socket)
  
  omp_indx = cel_master_out_thread+1
  !start omp threads
  !$omp parallel default(private)&
    !$omp& shared(mpi_comm, mpi_size,mpi_rank,mpi_type_ik,&
      !$omp& mpi_type_rk,&
      !$omp& num_elements,global_num_elements,&
      !$omp& num_rows,global_num_rows,&
      !$omp& cel_cg,omp,sparse_matrix_coo,sparse_matrix_coo_global,&
      !$omp& boundary_vector,solution_vector,th_solution_vector,&
      !$omp& num_threads,input_parameters)
  
  omp_indx = omp_get_thread_num() + 1
  !call  kmp_create_affinity_mask(thread_mask)
  !if(socket_id .eq. 0) then
  !  if(kmp_set_affinity_mask_proc( &
  !     INT(omp_indx-1), thread_mask) .NE. 0) then
  !     write (*,'(4(A,I0))') "rank: ", mpi_rank,"; thread_num: ",omp_indx-1, "; socket: ",socket_id, "; core:", omp_indx-1
  !     call cel_error("Error in Main kmp_set_affinity_mask_proc",&
  !     1_ik, .TRUE._lk, .TRUE._lk)
  !  endif
  !else
  !  if(kmp_set_affinity_mask_proc( &
  !     INT(omp_indx-1)+cel_cpu_param_cores, thread_mask) .NE. 0) then
  !     write (*,'(4(A,I0))') "rank: ", mpi_rank,"; thread_num: ",omp_indx-1, "; socket: ",socket_id, "; core:", omp_indx-1+num_procs_per_socket
  !     call cel_error("Error in Main kmp_set_affinity_mask_proc",&
  !     1_ik, .TRUE._lk, .TRUE._lk)
  !  endif
  !endif
  !if(kmp_set_affinity(thread_mask) .NE. 0) then
  !      call cel_error("Error in Main kmp_set_affinity",&
  !      1_ik, .TRUE._lk, .TRUE._lk)
  !end if
  !initialize thread and processes topology omp
  call init(omp(omp_indx), num_threads, omp_indx-1,&
      mpi_size, mpi_rank, mpi_comm, mpi_type_ik, mpi_type_rk, input_parameters%verbosity)
  !$omp barrier

  block_indx = cel_index_undef
  if(.not. (omp(omp_indx)%num_workers .EQ. omp_get_num_threads()-1) ) then
      call cel_error("MAIN check number of threads: could not run the request number of threads", &
        &1_ik,.TRUE._lk,.TRUE._lk)
  end if
  if(omp(omp_indx)%is_master) then
    output_on = input_parameters%verbosity .and. omp(omp_indx)%is_master .and. mpi_rank .EQ. 0
    to_stop = .TRUE._lk
    test_idx = 1_ik
  end if
  !$omp barrier
  
  !make coo matrix with poisson
  time_start = get_time()
  if(input_parameters%generate_matrix .EQ. 1) then
    if(omp(omp_indx)%is_master) then
      block_indx = omp(omp_indx)%master_num
      call cel_sp_mat_generate(sparse_matrix_coo, &
        boundary_vector, th_solution_vector, &
        input_parameters, omp(omp_indx), test_idx, err_code, &
        block_indx=block_indx)
      num_rows = sparse_matrix_coo%rows_num
      num_elements = sparse_matrix_coo%elements_num
      call MPI_ALLREDUCE(num_rows, global_num_rows, 1_mpik, mpi_type_ik, MPI_SUM, mpi_comm, ierr)
      call MPI_ALLREDUCE(num_elements, global_num_elements, 1_mpik, mpi_type_ik, MPI_SUM, mpi_comm, ierr)
      call cel_error_proc("MPI_ALLREDUCE", mpi_rank, ierr, output_on, to_stop)
    endif
  else if(input_parameters%generate_matrix .EQ. 2) then
     if(omp(omp_indx)%is_master) then
      block_indx = omp(omp_indx)%master_num
      call cel_sp_mat_generate_simple_test(sparse_matrix_coo, &
        boundary_vector, th_solution_vector, &
        input_parameters, omp(omp_indx), test_idx, err_code)
      num_rows = sparse_matrix_coo%rows_num
      num_elements = sparse_matrix_coo%elements_num
      call MPI_ALLREDUCE(num_rows, global_num_rows, 1_mpik, mpi_type_ik, MPI_SUM, mpi_comm, ierr)
      call MPI_ALLREDUCE(num_elements, global_num_elements, 1_mpik, mpi_type_ik, MPI_SUM, mpi_comm, ierr)
      call cel_error_proc("MPI_ALLREDUCE", mpi_rank, ierr, output_on, to_stop)
    end if
  else if(input_parameters%read_matrix .EQ. 1) then
      if(mpi_rank .eq. 0_mpik) then
       if(omp(omp_indx)%is_master) then
        block_indx = omp(omp_indx)%master_num
        call cel_sp_mat_read_matrix_market_file(sparse_matrix_coo_global, &
            input_parameters%matrix_data_dir%directory_name, &
            input_parameters%matrix_data_dir%file_name, &
            input_parameters%block_rows_num,&
            input_parameters%block_cols_num,&
            omp(omp_indx),&
            output_on, err_code)
        end if
      end if
      !$omp barrier
      if(omp(omp_indx)%is_master) then
        buffer(1)=sparse_matrix_coo_global%elements_num
        buffer(2)=sparse_matrix_coo_global%rows_num
        buffer(3)=sparse_matrix_coo_global%cols_num
        buffer(4)=sparse_matrix_coo_global%block_rows_num
        buffer(5)=sparse_matrix_coo_global%block_cols_num
        call MPI_Bcast(buffer,5_ik,omp(omp_indx)%mpi_type_ik,0_mpik,omp(omp_indx)%master_comm,ierr)
        sparse_matrix_coo_global%elements_num=buffer(1)
        sparse_matrix_coo_global%rows_num=buffer(2)
        sparse_matrix_coo_global%cols_num=buffer(3)
        sparse_matrix_coo_global%block_rows_num=buffer(4)
        sparse_matrix_coo_global%block_cols_num=buffer(5)
       ! call cel_sp_mat_distribute(sparse_matrix_coo, sparse_matrix_coo_global,&
       !   omp(omp_indx), output_on, err_code)
      end if
      if(omp(omp_indx)%is_master) then
        call MPI_BARRIER(mpi_comm, ierr)
      end if
      !$omp barrier
  else
    !$omp barrier
    if(omp(omp_indx)%is_master) then
      call cel_error_proc("Error matrix can not be generated; check the command line", mpi_rank, 1_mpik, output_on, .FALSE._lk)
      call MPI_FINALIZE(ierr);
      stop
    endif
  end if
  !fill vtxdist
 ! if(input_parameters%generate_matrix .GT. 0) then
 !   if(omp(omp_indx)%is_master) then
 !     if(allocated(sparse_matrix_coo%vtxdist)) deallocate(sparse_matrix_coo%vtxdist)
 !     allocate(sparse_matrix_coo%vtxdist(omp(omp_indx)%num_masters+1_ik))
 !     sparse_matrix_coo%vtxdist(1)=1_ik
 !     call MPI_ALLGATHER(num_rows,1_mpik,omp(omp_indx)%mpi_type_ik,&
 !       sparse_matrix_coo%vtxdist(2),1_mpik,omp(omp_indx)%mpi_type_ik,&
 !       omp(omp_indx)%master_comm,ierr)
 !     call cel_error_proc("MPI_ALLREDUCE", mpi_rank, ierr, output_on, to_stop)
 !     do ii=1,omp(omp_indx)%num_masters
 !       sparse_matrix_coo%vtxdist(ii+1_ik)=sparse_matrix_coo%vtxdist(ii+1_ik)+sparse_matrix_coo%vtxdist(ii)
 !     end do
 !   end if
 ! endif
 ! if(output_on) call print("vtxdist",sparse_matrix_coo%vtxdist)
  
  !coo matrix with poisson is done
  !$omp barrier
  if(omp(omp_indx)%is_master) time_end  = get_time()
  if(output_on) then
      write (*,'(A,I0,1(A,E15.6),4(A,I0))') "mt generated(",omp(omp_indx)%master_num, &
      "); time:",time_end-time_start,&
      "; num_elements:",num_elements,&
      "; gl_num_elements:",global_num_elements,&
      "; num_rows:",num_rows,&
      "; gl_num_rows:",global_num_rows
      !"; iter_max:",input_parameters%iter_max
  end if
  !$omp barrier
  if(omp(omp_indx)%is_master) then
    cel_cg%iter_max = input_parameters%iter_max
    cel_cg%eps = input_parameters%eps
    if(cel_perf_on) then
      cel_perf_mv_test_counter = 0_ik
      cel_perf_vect_test_counter = 0_ik
      cel_perf_mv_mpi_test_counter = 0_ik
      cel_perf_vec_mpi_test_counter = 0_ik
      call reset(cel_perf_cg_sum)
    endif
  end if 

  if(omp(omp_indx)%is_master) call new(solution_vector,num_rows)
  if(omp(omp_indx)%is_master) err_code = craypat_on()
  !$omp barrier
  if(omp(omp_indx)%is_master) time_start = get_time()
  !CG START
  do ii=1, input_parameters%attempts_num

   ! call cel_sp_mat_distr_cg(solution_vector%values,&
   !   boundary_vector%values,&
   !   th_solution_vector%values,&
   !   cel_cel,&
   !   cel_cg,&
   !   omp(omp_indx),&
   !   output_on,&
   !   err_code)
  end do
  !CG END
  if(omp(omp_indx)%is_master) time_end  = get_time()
  if(omp(omp_indx)%is_master) err_code = craypat_off()
  !collect some performance statistics
  !$omp barrier
  if(omp(omp_indx)%is_master) then
    
      tmp_time = time_end-time_start
      call MPI_ALLREDUCE(tmp_time, max_time, 1_mpik, mpi_type_rk, MPI_MAX, mpi_comm, ierr)
      call MPI_ALLREDUCE(tmp_time, min_time, 1_mpik, mpi_type_rk, MPI_MIN, mpi_comm, ierr)
      call MPI_ALLREDUCE(tmp_time, average_time, 1_mpik, mpi_type_rk, MPI_SUM, mpi_comm, ierr)
      average_time = average_time / omp(omp_indx)%num_masters

    if(cel_perf_on) then
      call cel_perf_incr_sum_perf(cel_perf_cg,i_type_mv,cel_perf_cg_sum)
      call cel_perf_incr_sum_perf(cel_perf_cg,i_type_vect,cel_perf_cg_sum)
      mv_operations = real(cel_perf_cg_sum%add,kind=perf_rk)+&
        &real(cel_perf_cg_sum%mult,kind=perf_rk)+&
        &real(cel_perf_cg_sum%divide,kind=perf_rk)+&
        &real(cel_perf_cg_sum%sqr,kind=perf_rk)
      mv_mflops=(mv_operations/tmp_time)*1.d-6
      mv_gb_per_sec=((real(cel_perf_cg_sum%load,kind=perf_rk)+real(cel_perf_cg_sum%store,kind=perf_rk)) / &
        tmp_time)*1.d-9
      call MPI_ALLREDUCE(mv_mflops, gl_mv_mflops, 1_mpik, mpi_type_rk, MPI_SUM, mpi_comm, ierr)
      call MPI_ALLREDUCE(mv_gb_per_sec, gl_mv_gb_per_sec, 1_mpik, mpi_type_rk, MPI_SUM, mpi_comm, ierr)
     ! if(omp(omp_indx)%master_num .eq. 0) then
      if(output_on) then
        write(*,*) "GLOBAL PERFORMANCE:"
        write(*,'(3(A,E13.6),(1(A,I10)),(A,E13.6))')  &
         "time:",tmp_time,&
         "; gl MFLOPS:", gl_mv_mflops,&
         "; gl GB/sec:", gl_mv_gb_per_sec,&
         "; iter: ", cel_cg%iter_done,&
         "; res.: ",cel_cg%residuum
        write (*,'(3(A,E13.6))') "av. time: ", average_time,&
         "; max time: ",max_time,"; min_time: ", min_time

      endif
      call MPI_BARRIER(mpi_comm, ierr)
      if(output_on) then
        write(*,'(A,I0,5(A,E13.6),(1(A,I0)),(A,E13.6))') "performance(",omp(omp_indx)%master_num, &
         ");time:",tmp_time,&
         "; MFLOPS:", mv_mflops,&
         "; GB/sec:", mv_gb_per_sec
      end if
    else if(output_on) then
      write(*,*) "GLOBAL PERFORMANCE:"
      write (*,'(3(A,E13.6))') "av. time: ", average_time,&
       "; max time: ",max_time,"; min_time: ", min_time
    endif
    
  endif
 !$omp barrier
  !if(omp(omp_indx)%is_master) then
  !  flag_error = .false._lk
  !  do ii=1, num_rows
  !    if(ABS(solution_vector%values(ii)-boundary_vector%values(ii)) .GE. 0.0000001) then
  !      flag_error = .true._lk
  !      !write (*,'(I0,A,I0,A,E16.6,A,E16.6)') omp(omp_indx)%master_num,&
  !      ! ": error in ii:",ii,"; calc:", vector_y%values(ii), "; th:",boundary_vector%values(ii)
  !    end if
  !  end do
!    if(flag_error) then
!      write (*,*) "MV with ERROR"
!    else
!      write (*,*) "NO ERROR"
!    endif
 ! end if
!  if(omp(omp_indx)%is_master) then
!          call MPI_BARRIER(mpi_comm, ierr); call cel_error_proc("MPI_BARRIER", mpi_rank, ierr, output_on, to_stop)
! endif

 !$omp end parallel
  !if(cel_perf_on .and. output_on) then
  !  write(perf_filename,'(A,I0)') trim("cel_cg_proc_"), omp(omp_indx)%master_num
  !  call cel_perf_write_performance(cel_perf_cg,perf_filename)
  !end if
 ! call del(sparse_matrix_coo)
 ! call del(boundary_vector)
 ! call del(th_solution_vector)
 ! call del(solution_vector)
  !deallocate(omp)
  CALL MPI_FINALIZE(ierr);
  
end program cel_cg_program

