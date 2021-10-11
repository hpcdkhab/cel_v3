
program cel_cg_program
!library fbase
use cel_types_module
use cel_input_parameters_module
use cel_omp_module
use cel_timer_interface_module
use cel_perf_module
use cel_perf_distr_module
use cel_cpu_param_module
use cel_vtkio_module
use cel_omp_shared_work_module
use cel_profile_module
!library fcomm
use cel_comm_module
!library fspmt
use cel_sp_mat_generate_module
use cel_sp_mat_module
!library fsp_mat_distr
use cel_sp_mat_distr_vec_module
use cel_sp_mat_distr_vec_module
use cel_sp_mat_distr_check_module
!library fspmtconverter
use cel_sp_mat_jad_converter_module
!library fop
use cel_fop_mv_module
use cel_fop_vec_distr_module
use cel_alg_module
use cel_sp_mat_distr_gl_module
use cel_cgalg_module
use cel_gmresalg_module
use cel_sp_mat_distr_format_module
!external libraries
!use io_benchmark_main_module, only :io_get_matrix_in_format_sparse_31_one_proc
!mpi and cel_omps libraries
use OMP_LIB
use MPI
use, intrinsic :: iso_c_binding
implicit none
  
  type(cel_input_parameters_type)               :: input_parameters !user input
  integer(kind = mpik)                          :: ierr, mpi_rank,cel_mpi_comm !mpi parameters
  integer(kind = mpik)                          :: mpi_size, mpi_type_ik, ik_mpi !mpi parameters
  integer(kind = mpik)                          :: mpi_type_rk, rk_mpi, mpi_type_perf_ik, perf_ik_mpi, mpi_provided !mpi parameters
  logical(kind=lk)                              :: output_on !output parameters
  logical(kind=lk)                              :: to_stop !debuger parameters
  integer(kind=ik)                              :: test_idx !perf. evaluation
  integer (kind=ik)                             :: curr_num_threads !threads number
  type(cel_omp_type), dimension(:), allocatable :: cel_omps !parameters for threads; thread num, master or worker thread and so on
  type(cel_omp_shared_work_type)                :: cel_omp_shared_work !manages the number of the worker threads, which are involved in computational tasks
  integer(kind=ik)                              :: num_threads, omp_indx !OpenMP / threads control
  integer(kind=ik)                              :: block_indx !block_id of coo matrix (distibution parameter)
  type(cel_vec_type)                            :: boundary_vector !boundary vector of system to solve
  type(cel_vec_type)                            :: th_solution_vector !theoretical solution of the system
  type(cel_vec_type)                            :: global_calc_y_vector !theoretical solution of the system
  real(kind=rk)                                 :: time_start, time_end !perf. evaluation
  integer(kind=ik)                              :: num_rows, num_elements, global_num_elements, global_num_rows !matrix parameters
  integer(kind=ik)                              :: err_code !error variable
  integer(kind=c_int)                           :: num_threads_c, num_threads_check 
  integer(kind=ik)                              :: ii, num_nodes, core_id
  logical(kind=lk)                              :: flag_error
  type(cel_sp_mat_type), dimension(:), allocatable  :: sparse_matrix_a !sparse matrix to solve 
  type(cel_sp_mat_type)                         :: sparse_matrix_coo_global!filled by proc 0 (before the distribution of the matrix)
  type(cel_sp_mat_type)                         :: sparse_matrix_coo_local!filled by proc 0 (before the distribution of the matrix)
  real(kind=rk)                                 :: max_time, average_time, min_time, tmp_time !perf. evaluation
  real(kind=rk)                                 :: mv_max_time, mv_average_time, mv_min_time, mv_tmp_time !perf. evaluation
  real(kind=rk)                                 :: gl_bandwidth, gl_flops
  real(kind=rk)                                 :: mv_operations
  character(len=256)                            :: perf_filename
  integer(kind=ik)                              :: socket_id !0 or 2 in two sockets system
  integer(kind=ik)                              :: num_procs_per_socket
  !integer (kind=kmp_affinity_mask_kind)         :: thread_mask
  integer(kind=ik), dimension(5)                :: buffer
  integer(kind=ik), dimension(:), allocatable   :: vtxdist ! array of matrix distribution (i-th el. is the index of the first row on i-th process)
  integer(kind=ik), dimension(:), allocatable   :: procs_array! array with proc_id for the communication elements in sp_mat_local_array
  integer(kind=ik), dimension(:), allocatable   :: procs_array_first_idx ! array with first idx of data needed from the process in procs_array
  integer(kind=ik), dimension(:), allocatable   :: procs_array_last_idx  ! array with last idx of data needed from the process in procs_array
  integer(kind=ik), dimension(:,:), allocatable :: com_buffer_perm
  integer(kind=ik), dimension(:), allocatable   :: length_comm_buffer
  type(cel_comm_type)                           :: cel_comm
  !real(kind=rk) , allocatable, dimension(:,:)   :: vtk_nodes
  type(cel_sp_mat_distr_vec_type)               :: x_distr_vector,r_distr_vector,boundary_distr_vector,solution_distr_vector
  type(cel_sp_mat_distr_vec_type)               :: d_distr_vector, x_distr_vector_th
  REAL(kind = rk), ALLOCATABLE, DIMENSION(:)    :: matrix_diagonal
  REAL(kind = rk), ALLOCATABLE, DIMENSION(:)    :: matrix_up_triangl
  INTEGER(kind = ik), ALLOCATABLE, DIMENSION(:) :: column_indices
  INTEGER(kind = ik), ALLOCATABLE, DIMENSION(:) :: row_ptr
  INTEGER(kind=ik) ::row_ptr_length
  character(len=1024) :: temp_string
  character(len=1024) :: std_output_filename, err_output_filename
  type(cel_perf_counter_type)                   :: local_perf, total_perf, converter_perf
  integer :: iostat
  integer(kind=ik) :: start_time_sec,start_time_nanosec,end_time_sec,end_time_nanosec
  type(cel_profile_type), dimension(:), allocatable :: profiles
  real(kind=rk) :: norm2, h3, disc_err
  logical(kind=4) :: is_thread_main,to_scale
  type(cel_gmresalg_parameter_type)                              :: cel_gmres
  type(cel_cgalg_parameter_type)                                 :: cel_cg
  type(cel_sp_mat_distr_vec_type)                                :: cg_vector_x
  type(cel_sp_mat_distr_vec_type)                                :: cg_vector_r
  type(cel_sp_mat_distr_vec_type)                                :: cg_vector_d
  type(cel_sp_mat_distr_vec_type)                                :: cg_vector_q
  type(cel_sp_mat_distr_vec_type)                                :: cg_vector_v
  type(cel_sp_mat_distr_vec_type),dimension(:), allocatable      :: gmres_vector_v
  type(cel_sp_mat_distr_vec_type),dimension(:), allocatable      :: gmres_vector_z
  type(cel_sp_mat_distr_vec_type),dimension(:), allocatable      :: gmres_vector_w
  type(cel_sp_mat_distr_vec_type)                                :: cg_vector_tmp
  type(cel_perf_counter_type)                                    :: cg_perf_counter_mv
  type(cel_perf_counter_type)                                    :: cg_perf_counter_vec
  type(cel_profile_type),dimension(:), allocatable               :: cg_cel_profiles
  real(kind=rk), allocatable, dimension(:,:)        :: cel_gmres_hh
  real(kind=rk), allocatable, dimension(:)          :: cel_gmres_cc
  real(kind=rk), allocatable, dimension(:)          :: cel_gmres_ss
  real(kind=rk), allocatable, dimension(:)          :: cel_gmres_gamma
  real(kind=rk), allocatable, dimension(:)          :: cel_gmres_alpha
  type(cel_sp_mat_distr_vec_type)                                :: distr_diagonal
  type(cel_sp_mat_distr_vec_type)                                :: distr_off_diagonal
  type(cel_sp_mat_distr_vec_type)                                :: distr_diagonal_dominance
  integer(kind=ik)                                            :: iter_count
  real(kind=rk)                                               :: residuum
  call cel_input_parameters_read(input_parameters)!read user input parameters from command line
  
  !obtain and set number of threads
  curr_num_threads = input_parameters%cel_omp%num_threads
  num_threads = INT(curr_num_threads,kind=c_ik)  
  num_threads_c = INT(num_threads,kind=c_ik)  
  to_stop = .TRUE._lk
  
  call omp_set_dynamic(.false._lk)
  call omp_set_nested(.false._lk)
  call omp_set_num_threads(num_threads_c)
  
  cel_mpi_comm = MPI_COMM_WORLD
  ik_mpi = INT(ik, kind = mpik)
  rk_mpi = INT(rk, kind = mpik)
  perf_ik_mpi = INT(perf_ik, kind = mpik)
  mpi_provided=0
  call MPI_INIT_THREAD(MPI_THREAD_FUNNELED,mpi_provided,ierr)
!  call MPI_INIT(ierr)
  call MPI_COMM_RANK(cel_mpi_comm, mpi_rank, ierr)
  call MPI_COMM_SIZE(cel_mpi_comm, mpi_size, ierr)
  call MPI_TYPE_MATCH_SIZE(MPI_TYPECLASS_INTEGER, ik_mpi, mpi_type_ik, ierr)
  call MPI_TYPE_MATCH_SIZE(MPI_TYPECLASS_INTEGER, perf_ik_mpi, mpi_type_perf_ik, ierr)
  call MPI_TYPE_MATCH_SIZE(MPI_TYPECLASS_REAL, rk_mpi, mpi_type_rk, ierr)
  if(mpi_rank .eq. 0) write(*,'(A,I0)') "mpi_provided:",mpi_provided
  if(mpi_rank .eq. 0) write(*,'(A,I0)') "MPI_THREAD_FUNNELED:",MPI_THREAD_FUNNELED
  err_code = craypat_off()
  if (mpi_rank .EQ. 0) then
    cel_is_in_debug = .TRUE._lk
  else
    cel_is_in_debug = .TRUE._lk
  endif

  if(mpi_rank .EQ. 0 .and. .false._lk) then
    output_on = input_parameters%verbosity
    if(input_parameters%verbosity)  cel_cg%verbosity = .true._lk
    Write(std_output_filename,'(A)') "out"
    Write(err_output_filename,'(A)') "err"
    cel_output_unit = 3000
    open(unit=cel_output_unit,iostat=iostat,file=std_output_filename,action='READWRITE',status='REPLACE')
    if(iostat .gt. 0) then
      write(cel_output_unit,'(A)') &
          "Error: could not open file for the output:"
    end if
    cel_error_unit = 30001
    open(unit=cel_output_unit,iostat=iostat,file=err_output_filename,action='READWRITE',status='REPLACE')
    if(iostat .gt. 0) then
      write(cel_error_unit,'(A)') &
          "Error: could not open file for the error output:"
    end if
    call print(input_parameters)
  else
    output_on = .FALSE._lk
  endif
  !some checks of the input parameters
  if(input_parameters%generate_matrix .GT. 1) then
    if(.not. input_parameters%dx*input_parameters%dy*input_parameters%dz .EQ. mpi_size) then
      if (mpi_rank .EQ. 0) then
        write(cel_error_unit,'(A,A)') &
          "Error: please, check the input parameters:",& 
          " dx*dy*dz must be equal to the number of processes !"
      end if
      CALL MPI_FINALIZE(ierr)
      stop
    end if
    num_nodes = input_parameters%nx*input_parameters%ny*input_parameters%nz/&
    input_parameters%dx*input_parameters%dy*input_parameters%dz
  end if
  if(input_parameters%cel_omp%num_threads .LT. 2) then
    if (mpi_rank .EQ. 0) then
      write(cel_error_unit,'(A,A)') &
        "Error: please, check the input parameters:",&
        " the parameter num_threads must be greater than one!"
    end if
    CALL MPI_FINALIZE(ierr)
    stop
  end if
  
  !allocate memory for thread and processes topology cel_omps
  if(allocated(cel_omps)) deallocate(cel_omps)
  allocate(cel_omps(num_threads))
  
  num_procs_per_socket = cel_cpu_param_num_sockets
    
  omp_indx = cel_master_out_thread+1
  !start omp threads
  !$omp parallel default(none) private(output_on,omp_indx,err_code,ii,core_id,&
      !$omp& converter_perf)&
    !$omp& shared(cel_mpi_comm, mpi_size,mpi_rank,mpi_type_ik,&
      !$omp& mpi_type_rk,&
      !$omp& num_elements,global_num_elements,&
      !$omp& num_rows,global_num_rows,&
      !$omp& cel_omps,sparse_matrix_a,sparse_matrix_coo_local, sparse_matrix_coo_global,&
      !$omp& boundary_distr_vector,solution_distr_vector,th_solution_vector,&
      !$omp& global_calc_y_vector,&
      !$omp& x_distr_vector, r_distr_vector,d_distr_vector,&
      !$omp& num_threads,input_parameters,buffer,vtxdist,&
      !$omp& procs_array,cel_comm,&
      !$omp& boundary_vector,&
      !$omp& matrix_diagonal,  matrix_up_triangl,&
      !$omp& column_indices, row_ptr, row_ptr_length,&
      !$omp& cel_omp_shared_work,cel_cg,&
      !$omp& cel_output_unit, cel_error_unit,&
      !$omp& local_perf, total_perf, socket_id,&
      !$omp& profiles, norm2,com_buffer_perm,length_comm_buffer,&
      !$omp& procs_array_first_idx,procs_array_last_idx,&
      !$omp& ierr, ik_mpi, perf_ik_mpi,&
      !$omp& rk_mpi, mpi_provided,mpi_type_perf_ik,to_stop,test_idx,&
      !$omp& curr_num_threads,block_indx, time_start, time_end,&
      !$omp& num_threads_check,num_threads_c,flag_error,&
      !$omp& max_time, average_time, min_time, tmp_time,&
      !$omp& mv_max_time, mv_average_time, mv_min_time, mv_tmp_time,&
      !$omp& gl_bandwidth, gl_flops,mv_operations,perf_filename,&
      !$omp& num_procs_per_socket, temp_string,std_output_filename,&
      !$omp& err_output_filename,start_time_sec,start_time_nanosec,&
      !$omp& end_time_sec,end_time_nanosec,is_thread_main,cel_is_in_debug,&
      !$omp& cg_vector_x,cg_vector_r,cg_vector_d,cg_vector_q,cg_vector_tmp,&
      !$omp& cg_perf_counter_mv,cg_perf_counter_vec,cg_cel_profiles,cg_vector_v,&
      !$omp& gmres_vector_v,gmres_vector_z,gmres_vector_w,cel_gmres,h3,disc_err,to_scale,&
      !$omp& cel_gmres_hh,cel_gmres_cc,cel_gmres_ss,cel_gmres_gamma,cel_gmres_alpha,&
      !$omp& distr_diagonal,distr_off_diagonal,distr_diagonal_dominance)
  core_id = omp_get_thread_num()
  omp_indx = omp_get_thread_num() + 1
  ii = int(mpi_rank,kind=ik)/num_procs_per_socket
    !initialize thread and processes topology cel_omps
  call init(cel_omps(omp_indx), num_threads, omp_indx-1,&
      mpi_size, mpi_rank, cel_mpi_comm, mpi_type_ik, mpi_type_rk, mpi_type_perf_ik, input_parameters%verbosity)
  !even
  !if( 2_ik*ii .eq. mpi_rank) then
  !  err_code = setThreadPolicy_sched(int(core_id,kind=4))
  !  write(cel_output_unit,'(4(A,I0))')&
  !    "rank:", mpi_rank,"; worker:",cel_omps(omp_indx)%worker_num,"; core_id:",core_id,"; err_code=",err_code 
  !else
  !  err_code = setThreadPolicy_sched(int(core_id+cel_cpu_param_cores,kind=4))
  !  write(cel_output_unit,'(4(A,I0))')&
  !    "rank:", mpi_rank,"; worker:",cel_omps(omp_indx)%worker_num,"; core_id:",core_id,"; err_code=",err_code 
  !end if

  ! PINNING - ONLY INTEl COMPILER
  !call  kmp_create_affinity_mask(thread_mask)
  !if(socket_id .eq. 0) then
  !  if(kmp_set_affinity_mask_proc( &
  !     INT(omp_indx-1), thread_mask) .NE. 0) then
  !     call cel_error("Error in Main kmp_set_affinity_mask_proc",&
  !     1_ik, .TRUE._lk, .TRUE._lk)
  !  endif
  !  write (cel_output_unit,'(4(A,I0))') "rank: ", mpi_rank,"; thread_num: ",omp_indx-1, "; socket: ",socket_id, "; core:", omp_indx-1
  !else
  !  if(kmp_set_affinity_mask_proc( &
  !     INT(omp_indx-1)+cel_cpu_param_cores, thread_mask) .NE. 0) then
  !     call cel_error("Error in Main kmp_set_affinity_mask_proc",&
  !     1_ik, .TRUE._lk, .TRUE._lk)
  !  endif
  !  write (cel_output_unit,'(4(A,I0))') "rank: ", mpi_rank,"; thread_num: ",omp_indx-1, "; socket: ",socket_id, "; core:", omp_indx-1+cel_cpu_param_cores
  !endif
  !if(kmp_set_affinity(thread_mask) .NE. 0) then
  !      call cel_error("Error in Main kmp_set_affinity",&
  !      1_ik, .TRUE._lk, .TRUE._lk)
  !end if
  

  !$omp barrier
  if(cel_omps(omp_indx)%is_master) then
     call MPI_IS_THREAD_MAIN(is_thread_main, ierr)
     if(.not. is_thread_main) then
       call cel_error_proc("Warning in main : master thread is not mpi main thread",&
        mpi_rank, ierr, .TRUE._lk, .TRUE._lk)
     end if
  end if
  call new(cel_omp_shared_work,cel_omps(omp_indx))
  !$omp barrier
  block_indx = cel_index_undef
  if(.not. (cel_omps(omp_indx)%num_workers .EQ. omp_get_num_threads()-1) ) then
      call cel_error("MAIN check number of threads: could not run the request number of threads", &
        &1_ik,.TRUE._lk,.TRUE._lk)
  end if
  if(cel_omps(omp_indx)%is_master) then
    output_on = input_parameters%verbosity .and. cel_omps(omp_indx)%is_master .and. mpi_rank .EQ. 0
    to_stop = .TRUE._lk
    test_idx = 1_ik
    if(cel_omps(omp_indx)%master_num .eq. 0_ik) then
      num_threads_check =  omp_get_num_threads()
      write (cel_output_unit,'(A,I0)') "num threads:", num_threads_check
    end if
  end if
  !$omp barrier
  
  !generate coo matrix 
  if(input_parameters%generate_matrix .GT. 0) then
    if(cel_omps(omp_indx)%is_master) time_start = get_time()
    if(input_parameters%generate_matrix .EQ. 1) then
      if(cel_omps(omp_indx)%is_master) then
        block_indx = cel_omps(omp_indx)%master_num
        call cel_sp_mat_generate(sparse_matrix_coo_local, &
          boundary_vector, th_solution_vector, &
          input_parameters, cel_omps(omp_indx), test_idx, err_code, &
          block_indx=block_indx)
      endif
    else if(input_parameters%generate_matrix .EQ. 2) then
       if(cel_omps(omp_indx)%is_master) then
        block_indx = cel_omps(omp_indx)%master_num
        call cel_sp_mat_generate_simple_test(sparse_matrix_coo_local, &
          boundary_vector, th_solution_vector, &
          input_parameters, cel_omps(omp_indx), test_idx, err_code)
      end if
    end if
    if(cel_omps(omp_indx)%is_master) then
      !fill vtxdist
       num_rows = sparse_matrix_coo_local%rows_num
      if(allocated(vtxdist)) deallocate(vtxdist)
      allocate(vtxdist(cel_omps(omp_indx)%num_masters+1_ik))
      vtxdist(1)=1_ik
      call MPI_ALLGATHER(num_rows,1_mpik,cel_omps(omp_indx)%mpi_type_ik,&
        vtxdist(2:size(vtxdist)),1_mpik,cel_omps(omp_indx)%mpi_type_ik,&
        cel_omps(omp_indx)%master_comm,ierr)
      call cel_error_proc("Error in main MPI_ALLREDUCE", mpi_rank, ierr, output_on, to_stop)
      do ii=1,cel_omps(omp_indx)%num_masters
        vtxdist(ii+1_ik)=vtxdist(ii+1_ik)+vtxdist(ii)
      end do
      !write the distributed matrix in vtk format
      !Warning large IO output; default: off
      if(cel_omps(omp_indx)%num_masters .lt. 33) then
        if(cel_omps(omp_indx)%is_master .and. input_parameters%write_matrix) then
          write(temp_string,'(A,A,I0,A,I0,A)') trim(input_parameters%matrix_data_dir%directory_name),&
          trim("/"),mpi_size,trim("_matrix_"),mpi_rank,trim(".vtk")
          call write_vtk_coo_matrix(sparse_matrix_coo_local%cols,&
            sparse_matrix_coo_local%rows,&
            sparse_matrix_coo_local%values, temp_string, .true.)
        end if
        if(cel_omps(omp_indx)%is_master .and. input_parameters%write_comm_graph) then
          !write(temp_string,'(I0,A,I0,A)') mpi_size,"_matrix_",mpi_rank,".vtk"
          !call write_vtk_coo_matrix(sparse_matrix_coo_local%cols,&
          !  sparse_matrix_coo_local%rows,&
          !  sparse_matrix_coo_local%values, temp_string, .true.)
          write(temp_string,'(A,I0)') trim(input_parameters%comm_graph_dir%directory_name),&
            trim("/"),cel_omps(omp_indx)%num_masters
          call cel_comm_show_send_group_graphvizio(temp_string,&
            "comm_group_of_proc_", "grapviz", cel_comm, cel_omps(omp_indx), err_code)
        end if
      else
        call cel_error_proc("Warning: The comm group graph was not exported: too many processors",&
         mpi_rank, ierr, output_on, .FALSE._lk)
      end if
    end if
  endif
  !$omp barrier


  if(input_parameters%read_matrix .GT. 0_ik) then
      if(mpi_rank .eq. 0_mpik) then
        if(cel_omps(omp_indx)%is_master) then
          if(input_parameters%read_matrix .EQ. 1_ik) then
           !read coo matrix in market format
            block_indx = cel_omps(omp_indx)%master_num
            call cel_sp_mat_read_matrix_market_file(sparse_matrix_coo_global, &
                input_parameters%matrix_data_dir%directory_name, &
                input_parameters%matrix_data_dir%file_name, &
                input_parameters%block_rows_num,&
                input_parameters%block_cols_num,&
                cel_omps(omp_indx),&
                err_code)
            allocate(th_solution_vector%values(sparse_matrix_coo_global%rows_num))
            allocate(global_calc_y_vector%values(sparse_matrix_coo_global%rows_num))
            th_solution_vector%values = 0.0_rk
            allocate(boundary_vector%values(sparse_matrix_coo_global%rows_num))
            boundary_vector%values = 0.0_rk
          elseif(input_parameters%read_matrix .EQ. 2) then
            call cel_error("Error in main: to read fmps matrix, please link the wp4 libraries!", &
             &1_ik,cel_is_in_debug,.TRUE._lk)
              !read coo matrix in fmps format
           !   call io_get_matrix_in_format_sparse_31_one_proc(input_parameters%matrix_data_dir%directory_name, .FALSE._lk, &
           !     &matrix_diagonal, matrix_up_triangl, column_indices, row_ptr, &
           !     &row_ptr_length, boundary_vector%values)
           !   call cel_sp_mat_from_format_sparse_31_full(sparse_matrix_coo_global, matrix_diagonal,&
           !     matrix_up_triangl, column_indices, row_ptr, &
           !     input_parameters%block_rows_num, input_parameters%block_cols_num,output_on, err_code)
           !   call cel_sp_mat_format_sparse_31_sort(sparse_matrix_coo_global, output_on, err_code)
           !   allocate(th_solution_vector%values(sparse_matrix_coo_global%rows_num))
           !   th_solution_vector%values = 1.0_rk
          end if
        end if
      endif
      if(cel_omps(omp_indx)%is_master) then
        call MPI_BARRIER(cel_mpi_comm, ierr)
      end if
      !$omp barrier
      !fill sp_mat_local with the local part of the global matrix sp_mat_global
      if(cel_omps(omp_indx)%is_master) then
        buffer(1)=sparse_matrix_coo_global%elements_num
        buffer(2)=sparse_matrix_coo_global%rows_num
        buffer(3)=sparse_matrix_coo_global%cols_num
        buffer(4)=sparse_matrix_coo_global%block_rows_num
        buffer(5)=sparse_matrix_coo_global%block_cols_num
        call MPI_Bcast(buffer,5_ik,cel_omps(omp_indx)%mpi_type_ik,0_mpik,cel_omps(omp_indx)%master_comm,ierr)
        sparse_matrix_coo_global%elements_num=buffer(1)
        sparse_matrix_coo_global%rows_num=buffer(2)
        sparse_matrix_coo_global%cols_num=buffer(3)
        sparse_matrix_coo_global%block_rows_num=buffer(4)
        sparse_matrix_coo_global%block_cols_num=buffer(5)
      end if
      call cel_sp_mat_distr_local(sparse_matrix_coo_local, vtxdist, sparse_matrix_coo_global,&
        cel_omps(omp_indx), err_code)
      call cel_error("Error in main: cel_sp_mat_distr_local!", &
        &err_code,cel_is_in_debug,.FALSE._lk)
      !write the global matrix in vtk format
      !Warning lagre IO output; default: off
      if(cel_omps(omp_indx)%num_masters .lt. 33) then
        if(cel_omps(omp_indx)%is_master .and. input_parameters%write_matrix) then
          write(temp_string,'(I0,A,I0,A)') mpi_size,"_matrix_",mpi_rank,".vtk"
    !      if(mpi_rank .eq. 0_mpik) then
            call write_vtk_coo_matrix(sparse_matrix_coo_global%cols,&
              sparse_matrix_coo_global%rows,&
              sparse_matrix_coo_global%values, temp_string, .true.)
    !     end if
        end if
      end if
  end if



  !calculate global number of rows and non-zero elements
  if(cel_omps(omp_indx)%is_master) then
    num_rows = sparse_matrix_coo_local%rows_num
    num_elements = sparse_matrix_coo_local%elements_num
    call MPI_ALLREDUCE(num_rows, global_num_rows, 1_mpik, mpi_type_ik, MPI_SUM, cel_mpi_comm, ierr)
    call MPI_ALLREDUCE(num_elements, global_num_elements, 1_mpik, mpi_type_ik, MPI_SUM, cel_mpi_comm, ierr)
    call cel_error_proc("MPI_ALLREDUCE", mpi_rank, ierr, output_on, to_stop)
  endif
  !$omp barrier
  !calculate who is the neighbours processes (involved in the communication)
  !and convert matrix in csr format
  call cel_sp_mat_distr_from_sp_mat_local(sparse_matrix_a, procs_array, &
         procs_array_first_idx, procs_array_last_idx,&
         com_buffer_perm,length_comm_buffer,&
         vtxdist, sparse_matrix_coo_local,&
         cel_omps(omp_indx), input_parameters, err_code)
  call cel_error("Error in main: cel_sp_mat_distr_from_sp_mat_local!", &
         &err_code,cel_is_in_debug,.FALSE._lk)

  call cel_comm_get_connections(cel_comm, procs_array,&
         procs_array_first_idx, procs_array_last_idx,&
         com_buffer_perm,length_comm_buffer,&
         cel_omps(omp_indx), err_code)
  !if(cel_omps(omp_indx)%master_num .eq. 0 ) then
  !  call cel_comm_print("proc 1",cel_comm,cel_omps(omp_indx))
  !end if
  call cel_error("Error in main: cel_comm_get_connections!", &
         &err_code,cel_is_in_debug,.FALSE._lk)
  call cel_comm_init_vtxdist_counts(cel_comm, vtxdist, cel_omps(omp_indx), err_code)
  call cel_error("Error in main: cel_comm_init_vtxdist_counts!", &
         &err_code,cel_is_in_debug,.FALSE._lk)
  !create the communicators for the halo data transfer
  !if(cel_omps(omp_indx)%is_master) then
  !write(*,"(3(A,I0))") "input_parameters%mv_communicator:",input_parameters%mv_communicator,&
  !      "; MV_COMM_ISEND_GROUPS:",MV_COMM_ISEND_GROUPS,&  
  !      "; MV_COMM_ISEND_GROUPS:",MV_COMM_IBCAST
  if(input_parameters%mv_communicator .eq. MV_COMM_ISEND_GROUPS .or. &
    input_parameters%mv_communicator .eq. MV_COMM_IBCAST) then
    call cel_comm_get_communicators(cel_comm, cel_omps(omp_indx), err_code)
    call cel_error("Error in main: cel_comm_get_communicators!", &
         &err_code,cel_is_in_debug,.FALSE._lk)
  else
    if(cel_omps(omp_indx)%is_master) then
      allocate(cel_comm%get_comms(0))
    end if
  endif
  !end if
  !$omp barrier
  if(cel_omps(omp_indx)%is_master) time_end  = get_time()
  if(output_on .and. cel_omps(omp_indx)%is_master) then
      write (cel_output_unit,'(A,I0,A,E15.6)',advance="no") "mt generated(",cel_omps(omp_indx)%master_num, &
      "); time:",time_end-time_start
      write (cel_output_unit,'(6(A,I0))') "; num_elements:",num_elements,&
      "; gl_num_elements:",global_num_elements,&
      "; num_rows:",num_rows,&
      "; gl_num_rows:",global_num_rows,&
      "; procs:",cel_omps(omp_indx)%num_masters,&
      "; iter_max:",input_parameters%iter_max
  end if
  !write the comm groups in graphvizio format
  !Warning large IO output; default: off
  if(cel_omps(omp_indx)%num_masters .lt. 33) then
    if(cel_omps(omp_indx)%is_master .and. input_parameters%write_comm_graph) then
       write(temp_string,'(A,A,I0)') trim(input_parameters%comm_graph_dir%directory_name),&
         trim("/"),cel_omps(omp_indx)%num_masters
      call cel_comm_show_send_group_graphvizio(temp_string,&
         "comm_group_of_proc_", "grapviz", cel_comm, cel_omps(omp_indx), err_code)
   end if
  else
    if(cel_omps(omp_indx)%is_master .and. input_parameters%write_comm_graph) then
      call cel_error_proc("Warning: The comm group graph was not exported: too many processors",&
        mpi_rank, ierr, output_on, .FALSE._lk)
    end if
  end if
  
  !$omp barrier
   if(cel_omps(omp_indx)%is_master) time_start = get_time()
  call reset(converter_perf)
   
  call cel_sp_mat_distr_set_mv_parameters(sparse_matrix_a,cel_omps(omp_indx),&
                               input_parameters,err_code)
  if(cel_omps(omp_indx)%is_master) call MPI_BARRIER(cel_mpi_comm, ierr)
  if(cel_omps(omp_indx)%is_master .and. cel_omps(omp_indx)%master_num .eq. 0_ik) then
   write(cel_output_unit,'(A,I0)') "cel_sp_mat_distr_set_mv_parameters"
  end if
  if(input_parameters%matrix_scaler .eq. MATRIX_SCALER_DIAGONAL) then
      call cel_sp_mat_distr_scale(sparse_matrix_a,&
                            cel_omps(omp_indx),cel_omp_shared_work,local_perf, output_on,err_code)
  end if
  if(cel_omps(omp_indx)%is_master) call MPI_BARRIER(cel_mpi_comm, ierr)
  if(cel_omps(omp_indx)%is_master .and. cel_omps(omp_indx)%master_num .eq. 0_ik) then
   write(cel_output_unit,'(A,I0)') "matrix_row_diagonal_scaling"
  end if
  call cel_sp_mat_distr_format(sparse_matrix_a,vtxdist,cel_omp_shared_work,&
                             cel_omps(omp_indx),local_perf,output_on,err_code)
  if(cel_omps(omp_indx)%is_master) call MPI_BARRIER(cel_mpi_comm, ierr)
  if(cel_omps(omp_indx)%is_master .and. cel_omps(omp_indx)%master_num .eq. 0_ik) then
   write(cel_output_unit,'(A,I0)') "cel_sp_mat_distr_format"
  end if
  if(input_parameters%preconditioner .eq. PRECONDITIONER_JACOBI) then
     to_scale = input_parameters%matrix_scaler .eq. MATRIX_SCALER_DIAGONAL
     call cel_sp_mat_distr_set_jacobi(sparse_matrix_a,to_scale,&
                              cel_omps(omp_indx),cel_omp_shared_work,&
                              local_perf, output_on,err_code)
  end if
  if(cel_omps(omp_indx)%is_master) call MPI_BARRIER(cel_mpi_comm, ierr)
  if(cel_omps(omp_indx)%is_master .and. cel_omps(omp_indx)%master_num .eq. 0_ik) then
   write(cel_output_unit,'(A,I0)') "cel_sp_mat_distr_set_jacobi"
  end if
  !=============================================
  !=======distr. sp. matrix is done==========
  !=============================================
  if(cel_omps(omp_indx)%is_master) time_end  = get_time()
 ! if(cel_omps(omp_indx)%master_num .eq. 1 .and. cel_omps(omp_indx)%is_master) then
 !   call print("sparse matrix",sparse_matrix_a(1))
 ! end if

  if(input_parameters%generate_matrix .GT. 0) then
   !if(cel_omps(omp_indx)%is_master) then
   !     call new(y_distr_vector,1_ik,sparse_matrix_a(1_ik)%rows_num)
   !     y_distr_vector%values(:,1_ik) = 0.0_rk
   ! end if
    call new(x_distr_vector,cel_comm%get_from,vtxdist, cel_omps(omp_indx))
    call new(r_distr_vector,cel_comm%get_from,vtxdist, cel_omps(omp_indx))
    call new(d_distr_vector,cel_comm%get_from,vtxdist, cel_omps(omp_indx))
    call new(boundary_distr_vector,cel_comm%get_from,vtxdist, cel_omps(omp_indx))
    call new(solution_distr_vector,cel_comm%get_from,vtxdist, cel_omps(omp_indx))
    call cel_sp_mat_distr_local_vec(boundary_distr_vector,boundary_vector%values,vtxdist,cel_omps(omp_indx))
    call cel_sp_mat_distr_local_vec(solution_distr_vector,th_solution_vector%values,vtxdist,cel_omps(omp_indx))
  else if(input_parameters%read_matrix .GT. 0) then
  ! if(cel_omps(omp_indx)%is_master) then
  !      call new(y_distr_vector,1_ik,sparse_matrix_a(1_ik)%rows_num)
  !      y_distr_vector%values(:,1_ik) = 0.0_rk
  !  end if
    call new(x_distr_vector,cel_comm%get_from,vtxdist, cel_omps(omp_indx))
    call new(r_distr_vector,cel_comm%get_from,vtxdist, cel_omps(omp_indx))
    call new(d_distr_vector,cel_comm%get_from,vtxdist, cel_omps(omp_indx))
    call new(boundary_distr_vector,cel_comm%get_from,vtxdist, cel_omps(omp_indx))
    call new(solution_distr_vector,cel_comm%get_from,vtxdist, cel_omps(omp_indx))
    call cel_sp_mat_distr_global_vec(boundary_distr_vector,boundary_vector%values,cel_comm,vtxdist,cel_omps(omp_indx))
    call cel_sp_mat_distr_global_vec(solution_distr_vector,th_solution_vector%values,cel_comm,vtxdist,cel_omps(omp_indx))
  end if
  !call print("START ",distr_vectors,cel_omps(omp_indx),output_on)
  !call cel_sp_mat_distr_vec_start_all(distr_vectors,cel_comm,vtxdist,cel_omps(omp_indx),err_code)
  !call cel_sp_mat_distr_vec_end_all(cel_comm,cel_omps(omp_indx),err_code)
  !call print("END ",distr_vectors,cel_omps(omp_indx),cel_omps(omp_indx)%master_num .eq. 1_ik)
  !=============================================
  !========distributed vectors are created======
  !=============================================
  if(cel_omps(omp_indx)%is_master .and. input_parameters%write_performance) then
    call reset(local_perf)
    call reset(total_perf)
    total_perf%cpu_freq_idx = input_parameters%cpu_freq_idx
    local_perf%size = num_rows
    local_perf%num_elements = num_elements
  end if
  if(cel_omps(omp_indx)%is_master) then
    allocate(profiles(1))
  end if

    if(input_parameters%benchmark_id .eq. BENCHMARK_CG) then
      call new(cel_cg,vtxdist,&
        cg_vector_x,cg_vector_r,cg_vector_d,cg_vector_q,cg_vector_v,cg_vector_tmp,&
        local_perf,cg_perf_counter_mv, cg_perf_counter_vec,cg_cel_profiles,&
        cel_omps(omp_indx),cel_comm, input_parameters,err_code)
    elseif(input_parameters%benchmark_id .eq. BENCHMARK_GMRES) then
      call new(cel_gmres,vtxdist,&
        cg_vector_x,cg_vector_r,gmres_vector_v,gmres_vector_z,gmres_vector_w,cg_vector_tmp,&
        cel_gmres_hh,cel_gmres_cc,cel_gmres_ss,cel_gmres_gamma,cel_gmres_alpha,&
        local_perf,cg_perf_counter_mv, cg_perf_counter_vec,cg_cel_profiles,&
        cel_omps(omp_indx),cel_comm, input_parameters,err_code)
    end if

  !$omp barrier
  if(input_parameters%matrix_scaler .eq. MATRIX_SCALER_DIAGONAL) then
     call cel_sp_mat_distr_boundary_diag_scale(boundary_distr_vector, sparse_matrix_a(1_ik)%scaled_diagonal,&
                            cel_omps(omp_indx),cel_omp_shared_work,local_perf, output_on,err_code)
     if(input_parameters%benchmark_id .eq. BENCHMARK_MV) then
       call cel_sp_mat_distr_vector_diag_scale_inv(solution_distr_vector, sparse_matrix_a(1_ik)%scaled_diagonal,&
                            cel_omps(omp_indx),cel_omp_shared_work,local_perf, output_on,err_code)
     end if
  end if

  if(cel_omps(omp_indx)%is_master) call MPI_BARRIER(cel_mpi_comm, ierr)
  !===================================
  !==============CG START=============
  !===================================
  if(cel_omps(omp_indx)%is_master) then
    cel_cg%iter_max = input_parameters%iter_max
    cel_cg%eps = input_parameters%eps
  end if 
  !$omp barrier
  if(cel_omps(omp_indx)%is_master) tmp_time = 0.0_rk
  if(cel_omps(omp_indx)%is_master .and. cel_omps(omp_indx)%master_num .eq. 0_ik) then
   write(cel_output_unit,'(A,I0)') "start benchmark(1-MV;2-R_D;3-NORM2;4-CG;):",&
    input_parameters%benchmark_id
  end if
  if(cel_omps(omp_indx)%is_master) call MPI_BARRIER(cel_mpi_comm, ierr)
  if(cel_omps(omp_indx)%is_master) then
   profiles(1)%id = cel_omps(omp_indx)%num_threads
   call get_time_stamp(profiles(1)%start_time_sec,profiles(1)%start_time_nanosec)
  end if
  !$omp barrier
  !!!!!!!!!!!!!!!!!
  !call cel_sp_mat_distr_check_diagonal_dominance(sparse_matrix_a,distr_diagonal,&
  !                              distr_off_diagonal,distr_diagonal_dominance,&
  !                              cel_comm,vtxdist,cel_omps(omp_indx),&
  !                              cel_omp_shared_work,&
  !                              local_perf,input_parameters,err_code)
  !call cel_sp_mat_distr_vec_print_sorted("Diagonal",distr_diagonal,&
   !                             cel_omps(omp_indx),cel_comm,&
   !                             input_parameters%verbosity,err_code)
 ! call cel_sp_mat_distr_vec_print_sorted("Off Diagonal",distr_off_diagonal,&
  !                              cel_omps(omp_indx),cel_comm,&
  !                              input_parameters%verbosity,err_code)
  !call cel_sp_mat_distr_vec_print_sorted("Diagonal dominance",distr_diagonal_dominance,&
   !                             cel_omps(omp_indx),cel_comm,&
   !                             input_parameters%verbosity,err_code)
  !!!!!!!!!!!!!!!!!
  if(cel_omps(omp_indx)%is_master) ierr=craypat_on()
  do ii=1, input_parameters%attempts_num
   !$omp barrier
    if(cel_omps(omp_indx)%is_master) time_start = get_time()
    select case (input_parameters%benchmark_id)
      case (BENCHMARK_MV)
        call cel_fop_vec_distr_to_zero(r_distr_vector,cel_omps(omp_indx),cel_omp_shared_work,&
                                local_perf, output_on,err_code)
        call cel_fop_mv_axy_distr(sparse_matrix_a,solution_distr_vector,r_distr_vector,&
                             d_distr_vector,cel_comm,vtxdist,cel_omps(omp_indx),cel_omp_shared_work,&
                             local_perf, output_on, err_code)
         call cel_error("Error in main  cel_fop_mv_axy_distr!", &
                         &err_code,cel_is_in_debug,.FALSE._lk)
      case (BENCHMARK_R_D)
        call cel_fop_vec_distr_calc_r_and_d(r_distr_vector,&
                                d_distr_vector,boundary_distr_vector,&
                                cel_omps(omp_indx),cel_omp_shared_work,&
                                local_perf, output_on,err_code)
      case (BENCHMARK_NORM2)
        call cel_fop_vec_distr_norm2_allreduce(norm2,r_distr_vector, &
                                cel_omps(omp_indx),cel_omp_shared_work,&
                                local_perf, output_on,err_code)
      case (BENCHMARK_CG)

        if(input_parameters%preconditioner .eq. PRECONDITIONER_JACOBI) then
         call cel_cgalg_jacobi(cel_cg,sparse_matrix_a,vtxdist, boundary_distr_vector, &
           x_distr_vector,cg_vector_r,cg_vector_d,cg_vector_q,cg_vector_v,cg_vector_tmp,&
           cel_omps(omp_indx),cel_comm,cel_omp_shared_work, &
           local_perf,cg_perf_counter_mv,cg_perf_counter_vec,cg_cel_profiles,&
           iter_count, residuum,&
           output_on,err_code)
       else
          call cel_cgalg(cel_cg,sparse_matrix_a,vtxdist, boundary_distr_vector, &
           x_distr_vector,cg_vector_r,cg_vector_d,cg_vector_q,&
           cg_vector_tmp,cel_omps(omp_indx),cel_comm,cel_omp_shared_work, &
           local_perf,cg_perf_counter_mv,cg_perf_counter_vec,cg_cel_profiles,&
           iter_count, residuum,&
           output_on,err_code)
       endif
      case (BENCHMARK_GMRES)
          if(input_parameters%preconditioner .eq. PRECONDITIONER_JACOBI) then
            call cel_gmresalg_2jacobi(cel_gmres,sparse_matrix_a,vtxdist,boundary_distr_vector,&
              x_distr_vector,cg_vector_r,gmres_vector_v,gmres_vector_z,gmres_vector_w,cg_vector_tmp,&
              cel_gmres_hh,cel_gmres_cc,cel_gmres_ss,cel_gmres_gamma,cel_gmres_alpha,&
              cel_omps(omp_indx),cel_comm,cel_omp_shared_work,&
              local_perf,cg_perf_counter_mv,cg_perf_counter_vec,cg_cel_profiles,&
              iter_count, residuum,&
              output_on,err_code)
          else              
            call cel_gmresalg_2(cel_gmres,sparse_matrix_a,vtxdist,boundary_distr_vector,&
              x_distr_vector,cg_vector_r,gmres_vector_v,gmres_vector_w,cg_vector_tmp,&
              cel_gmres_hh,cel_gmres_cc,cel_gmres_ss,cel_gmres_gamma,cel_gmres_alpha,&
              cel_omps(omp_indx),cel_comm,cel_omp_shared_work,&
              local_perf,cg_perf_counter_mv,cg_perf_counter_vec,cg_cel_profiles,&
              iter_count, residuum,&
              output_on,err_code)
          endif 
          call cel_error("Error in main  BENCHMARK_GMRES !", &
                         &err_code,cel_is_in_debug,.TRUE._lk) 
      case default
        err_code = -4523_ik
        call cel_error("Error in main  benchmark_id is undefined!", &
                         &err_code,cel_is_in_debug,.TRUE._lk)
    end select
    !$omp barrier
    if(cel_omps(omp_indx)%is_master) then
     tmp_time  = tmp_time + (get_time() - time_start)
     call MPI_BARRIER(cel_mpi_comm, ierr)
    end if
   end do
   if(cel_omps(omp_indx)%is_master .and. cel_omps(omp_indx)%master_num .eq. 0_ik) then
     write(cel_output_unit,'(A)') "benchmark finished"
   end if
  if(cel_omps(omp_indx)%is_master) then
    call get_time_stamp(profiles(1)%end_time_sec,profiles(1)%end_time_nanosec)
  end if
  if(cel_omps(omp_indx)%is_master) ierr=craypat_off()
 
  if(input_parameters%write_performance) then
      call cel_perf_distr_collect(total_perf,local_perf,cel_omps(omp_indx),.true._lk)
      total_perf%number_of_repetitions = input_parameters%attempts_num
  end if
  !===================================
  !==============CG END===============
  !===================================

  !collect some performance statistics and check result
  !$omp barrier
  if(cel_omps(omp_indx)%is_master) then
    call MPI_ALLREDUCE(tmp_time, max_time, 1_mpik, mpi_type_rk, MPI_MAX, cel_mpi_comm, ierr)
    call MPI_ALLREDUCE(tmp_time, min_time, 1_mpik, mpi_type_rk, MPI_MIN, cel_mpi_comm, ierr)
    call MPI_ALLREDUCE(tmp_time, average_time, 1_mpik, mpi_type_rk, MPI_SUM, cel_mpi_comm, ierr)
    average_time = average_time / real(mpi_size,kind=8)
    if(input_parameters%write_profile) then
      if(cel_omps(omp_indx)%master_num .EQ. 0_ik) then
       call cel_profile_write(profiles,input_parameters%profile_filename)
      end if
    end if
    if(input_parameters%write_performance) then
      if(cel_omps(omp_indx)%master_num .EQ. 0_ik) then
        mv_operations = real(total_perf%add,kind=perf_rk)+&
          &real(total_perf%mult,kind=perf_rk)+&
          &real(total_perf%divide,kind=perf_rk)+&
          &real(total_perf%sqr,kind=perf_rk)
        gl_flops=(mv_operations/average_time)
        gl_bandwidth=((real(total_perf%load,kind=perf_rk)+real(total_perf%store,kind=perf_rk)) / &
          average_time)
       ! if(cel_omps(omp_indx)%master_num .eq. 0) then
        if(.true._lk) then
          write(cel_output_unit,'(A)') CHAR(13)//CHAR(10)//"==AGREGATED PERFORMANCE COUNTERS=="
          call print(total_perf)
       
          write(cel_output_unit,'(A)') CHAR(13)//CHAR(10)//"==AGREGATED PERFORMANCE=="
          write(cel_output_unit,'(2(A,E13.6),A,I10,A,I10,A,I10,A,E13.6)')  &
                                    "GFlops:", gl_flops*1.d-9,&
               CHAR(13)//CHAR(10)//"GiB/sec:", gl_bandwidth*1.d-9,&
               CHAR(13)//CHAR(10)//"# sends:", total_perf%send_messages_num,&
               CHAR(13)//CHAR(10)//"# gets:", total_perf%get_messages_num,&
               CHAR(13)//CHAR(10)//"iter: ", 0_ik,&
               CHAR(13)//CHAR(10)//"res.: ",0.0_rk
          write (cel_output_unit,'(A,E13.6,A,E13.6,A,E13.6,A,E13.6,A,E13.6)')&
               "av. time: ", average_time,&
               CHAR(13)//CHAR(10)//"max time: ",max_time,&
               CHAR(13)//CHAR(10)//"min_time: ", min_time,&
               CHAR(13)//CHAR(10)//"counter time max: ",total_perf%time(4_ik) / real(mpi_size,kind=8),&
               CHAR(13)//CHAR(10)//"counter comm time max: ",total_perf%time_comm(4_ik) / real(mpi_size,kind=8)
          if(total_perf%time_comm(4) .gt. 0.) then 
              write (cel_output_unit,'(A,E13.6)') &
                     "total time/comm time: ",total_perf%time(4_ik)/total_perf%time_comm(4_ik)
          end if

          call cel_perf_write_performance_counters(total_perf,cel_omps(omp_indx),&
                                                    input_parameters%memory_level,&
                                                    input_parameters%performance_filename)
          if(output_on) then
            write (cel_error_unit,'(A)') CHAR(13)//CHAR(10)//"===TIME==="
            write (cel_error_unit,'(A,E13.6,A,E13.6,A,E13.6)')&
                  "av. time: ", average_time,&
                  CHAR(13)//CHAR(10)//"max time: ",max_time,&
                  CHAR(13)//CHAR(10)//"min_time: ", min_time
          endif
        end if
      end if
    end if
  end if
  !check result
  !$omp barrier
  
  if(input_parameters%check_result) then
    select case (input_parameters%benchmark_id)
      case (BENCHMARK_MV)
        if(cel_omps(omp_indx)%is_master ) then
         if(input_parameters%attempts_num .ne. 1_ik) then
          call cel_error("Error in main  can not check the result, number of attempt must be 1!", &
                             &err_code,cel_is_in_debug,.FALSE._lk)
         else
            flag_error = .false._lk
            do ii=1, sparse_matrix_a(1_ik)%rows_num
              if(ABS(r_distr_vector%values(ii,1_ik)-boundary_vector%values(ii)) .GE. 0.0000001) then
                flag_error = .true._lk
                !Warning large output
                if(.FALSE._lk) then
                  write (cel_output_unit,'(I0,A,I0,A,E16.6,A,E16.6)') cel_omps(omp_indx)%master_num,&
                  ": error in ii:",ii,"; calc:", r_distr_vector%values(ii,1_ik), "; th:",boundary_vector%values(ii)
                else
                    write (cel_output_unit,'(I0,A,I0,A,E16.6,A,E16.6)') cel_omps(omp_indx)%master_num,&
                  ": error in ii:",ii,"; calc:", r_distr_vector%values(ii,1_ik), "; th:",boundary_vector%values(ii)
                  exit
                end if
              end if
            end do
            if(flag_error) then
              write (cel_output_unit,'(A)') "MV with ERROR"
            else
              write (cel_output_unit,'(A)') "NO ERROR"
            endif
          end if
        end if
      case (BENCHMARK_GMRES)
        !d_(i+1)=r_(i+1)+betta*d_(i)
        call cel_fop_vec_distr_scale_1(solution_distr_vector,-1.0_rk, x_distr_vector,&
                                cel_omps(omp_indx),cel_omp_shared_work,&
                                local_perf, output_on,err_code)
        call cel_fop_vec_distr_norm2_allreduce(disc_err,solution_distr_vector, &
                                cel_omps(omp_indx),cel_omp_shared_work,&
                                local_perf, output_on,err_code)
        !$omp  barrier
        if(cel_omps(omp_indx)%is_master  .and. cel_omps(omp_indx)%master_num .EQ. 0_ik) then
          h3 =  sqrt((1.0_rk/(REAL((input_parameters%nx-1_ik)*&
            (input_parameters%nx-1_ik)*&
            (input_parameters%nz-1_ik),kind=rk))))
           write (cel_output_unit,'(A,E16.6)') &
            "CRESTA Linear solver CG disretization error ||x_theory-x_solve||*sqrt(hx*hy*hz):", &
            disc_err*sqrt(h3)
           write (cel_output_unit,'(A,E16.6)') &
            "CRESTA Linear solver CG disretization error per grid node:", &
            disc_err*sqrt(h3)/&
            REAL(input_parameters%nx*input_parameters%ny*input_parameters%nz,kind=rk)
        end if
      end select
    end if
  if(cel_omps(omp_indx)%is_master) then
    call MPI_BARRIER(cel_mpi_comm, ierr); 
    call cel_error_proc("MPI_BARRIER", mpi_rank, ierr, output_on, to_stop)
  endif
 !$omp end parallel


  call del(sparse_matrix_coo_global)
  call del(sparse_matrix_coo_local)
  call del(boundary_vector)
  call del(th_solution_vector)
  call del(sparse_matrix_a)
  call del(cel_comm)
  call del(x_distr_vector)
  call del(r_distr_vector)
  call del(d_distr_vector)
  call del(boundary_distr_vector)
  call del(solution_distr_vector)
  call del(cel_omp_shared_work)
  call del(cel_cg)
  call del(cg_vector_x)
  call del(cg_vector_r)
  call del(cg_vector_d)
  call del(cg_vector_q)
  call del(cg_vector_tmp)
  
  if(allocated(matrix_diagonal)) deallocate(matrix_diagonal)
  if(allocated(matrix_diagonal)) deallocate(matrix_up_triangl)
  if(allocated(matrix_diagonal)) deallocate(column_indices)
  if(allocated(matrix_diagonal)) deallocate(row_ptr)
 ! if(allocated(vtk_nodes)) deallocate(vtk_nodes)
  if(allocated(procs_array)) deallocate(procs_array)
  if(allocated(vtxdist)) deallocate(vtxdist)
  
  deallocate(cel_omps)
  CALL MPI_FINALIZE(ierr);
  
end program cel_cg_program

