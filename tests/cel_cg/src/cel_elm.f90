
program cel_elm
use cel_types_module
use cel_omp_module
use cel_timer_interface_module
use cel_perf_module
use cel_perf_distr_module
use cel_cpu_param_module
use cel_vtkio_module
use cel_omp_shared_work_module
use cel_profile_module
use cel_comm_module
use cel_sp_mat_generate_module
use cel_sp_mat_module
use cel_sp_mat_distr_vec_module
use cel_sp_mat_distr_vec_module
use cel_sp_mat_distr_check_module
use cel_sp_mat_jad_converter_module
use cel_fop_mv_module
use cel_fop_vec_distr_module
use cel_alg_module
use cel_sp_mat_distr_gl_module
use cel_cgalg_module
use cel_gmresalg_module
use cel_sp_mat_distr_format_module
use cel_input_parameters_module
use cel_elmsolver_module
use cel_timer_interface_module
!mpi and cel_omps libraries
use OMP_LIB
use MPI
use, intrinsic :: iso_c_binding
implicit none
  
  integer(kind = mpik) :: mpi_comm_elm
  integer(kind = mpik)                          :: ierr, mpi_rank !mpi parameters
  integer(kind = mpik)                          :: mpi_size, mpi_type_ik, ik_mpi !mpi parameters
  integer(kind = mpik)                          :: mpi_type_rk, rk_mpi, mpi_type_perf_ik, perf_ik_mpi, mpi_provided !mpi parameters
  type(cel_sp_mat_type)                         :: sparse_matrix_coo_global!filled by proc 0 (before the distribution of the matrix)
  type(cel_sp_mat_type)                         :: sparse_matrix_coo_local!filled by proc 0 (before the distribution of the matrix)
  integer(kind=ik)                              :: ii, core_id, omp_indx, num_procs_per_socket
  integer(kind=ik)                              :: err_code, num_threads
  type(cel_omp_type), dimension(:), allocatable :: cel_omps !parameters for threads; thread num, master or worker thread and so on
  type(cel_vec_type)                            :: boundary_vector !boundary vector of system to solve
  type(cel_vec_type)                            :: th_solution_vector !theoretical solution of the system
  integer(kind=ik), dimension(5)                :: buffer
  integer(kind=ik), dimension(:), allocatable   :: vtxdist 
  integer(kind=ik) :: num, cur_row, cur_col
  real(kind=rk) :: cur_value
  integer(kind=ik), dimension(2)                :: nnzero
  real(kind=rk), dimension(:), allocatable      :: bmat
  real(kind=rk), dimension(:), allocatable      :: xmat
  real(kind=perf_rk)                            :: time_start, time_end !perf. evaluation
  integer(kind=ik) :: iierr,niter
  real(kind=rk) :: ftol
  real(kind=rk) :: rand2
  call cel_input_parameters_read(input_parameters)!read user input parameters from command line
  !call cel_input_parameters_print(input_parameters)
  num_threads=1_ik  
  call omp_set_dynamic(.false._lk)
  call omp_set_nested(.false._lk)
  call omp_set_num_threads(num_threads)
  
  mpi_comm_elm = MPI_COMM_WORLD
  ik_mpi = INT(ik, kind = mpik)
  rk_mpi = INT(rk, kind = mpik)
  perf_ik_mpi = INT(perf_ik, kind = mpik)
  mpi_provided=0

  call MPI_INIT_THREAD(MPI_THREAD_FUNNELED,mpi_provided,ierr)
!  call MPI_INIT(ierr)
  call MPI_COMM_RANK(mpi_comm_elm, mpi_rank, ierr)
  call MPI_COMM_SIZE(mpi_comm_elm, mpi_size, ierr)
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

  !allocate memory for thread and processes topology cel_omps
  allocate(cel_omps(num_threads))
  
  num_procs_per_socket = cel_cpu_param_num_sockets
    
  omp_indx = cel_master_out_thread+1

  core_id = omp_get_thread_num()
  omp_indx = omp_get_thread_num() + 1
  ii = int(mpi_rank,kind=ik)/num_procs_per_socket
  !initialize thread and processes topology cel_omps
  call init(cel_omps(omp_indx), num_threads, omp_indx-1,&
      mpi_size, mpi_rank, mpi_comm_elm, mpi_type_ik, mpi_type_rk, mpi_type_perf_ik, input_parameters%verbosity)
  cel_omps(omp_indx)%is_master = .TRUE._lk
  cel_omps(omp_indx)%master_num = mpi_rank
  cel_omps(omp_indx)%num_masters = mpi_size
  if(input_parameters%read_matrix .GT. 0_ik) then
      if(mpi_rank .eq. 0_mpik) then
        if(input_parameters%read_matrix .EQ. 1_ik) then
         !read coo matrix in market format
          call cel_sp_mat_read_matrix_market_file(sparse_matrix_coo_global, &
              input_parameters%matrix_data_dir%directory_name, &
              input_parameters%matrix_data_dir%file_name, &
              input_parameters%block_rows_num,&
              input_parameters%block_cols_num,&
              cel_omps(omp_indx),&
              err_code)
          allocate(th_solution_vector%values(sparse_matrix_coo_global%rows_num))
          th_solution_vector%values = 0.0_rk
          allocate(boundary_vector%values(sparse_matrix_coo_global%rows_num))
          boundary_vector%values = 0.0_rk
        elseif(input_parameters%read_matrix .EQ. 2) then
          call cel_error("Error in main: to read fmps matrix, please link the wp4 libraries!", &
           &1_ik,cel_is_in_debug,.TRUE._lk)
        end if
      endif
      !fill sp_mat_local with the local part of the global matrix sp_mat_global
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
      call cel_sp_mat_distr_local(sparse_matrix_coo_local, vtxdist, sparse_matrix_coo_global,&
        cel_omps(omp_indx), err_code)
      call del(sparse_matrix_coo_global)
      !if(mpi_rank .eq. 0_mpik) then
      ! call print("p-matrix:",sparse_matrix_coo_local,.FALSE._lk)
      !end if
      call cel_error("Error in main: cel_sp_mat_distr_local!", &
        &err_code,cel_is_in_debug, .FALSE._lk)
  end if
 
  deallocate(cel_omps)
  
  num = vtxdist(mpi_rank+2)-vtxdist(mpi_rank+1)
  nnzero(1)=1_ik
  nnzero(2)= int(real((sparse_matrix_coo_local%elements_num / num),kind=rk)*(1.1_rk),kind=ik)
  call cel_elmSolverInitialization(sparse_matrix_coo_global%rows_num,&
    num,vtxdist(mpi_rank+1),nnzero,mpi_comm_elm)
 
  call MPI_Barrier(mpi_comm_elm,ierr)
  time_start = get_time()
  do ii=1_ik,sparse_matrix_coo_local%elements_num
    cur_row = sparse_matrix_coo_local%rows(ii)
    cur_col = sparse_matrix_coo_local%cols(ii)
    cur_value = sparse_matrix_coo_local%values(ii)
    call cel_elmSolverDatarowStart(cur_row)
    call cel_elmSolverDataInsertion(cur_col,cur_value)
  end do
  time_end = get_time()
  call MPI_Barrier(mpi_comm_elm,ierr)
  if(mpi_rank .eq. 0_mpik) then
    write (cel_output_unit,'(I0,A,E15.6)',advance="yes") mpi_rank, &
     ":data insertion time:",time_end-time_start
  end if

  !write(*,'(I0,A,I0)') mpi_rank,": main elm_num_local_elements:",elm_num_local_elements
  !write(*,'(I0,A,I0)') mpi_rank,": main elm_num_local_rows:",elm_num_local_rows
  
  call MPI_Barrier(mpi_comm_elm,ierr)
  time_start = get_time()
  call cel_elmSolverMatrixCompleted()
  time_end = get_time()
  call MPI_Barrier(mpi_comm_elm,ierr)

  if(mpi_rank .eq. 0_mpik) then
    write (cel_output_unit,'(I0,A,E15.6)',advance="yes") mpi_rank, &
     ":matrix completed time:",time_end-time_start
  end if

  allocate(bmat(sparse_matrix_coo_local%rows_num))
          
  call RANDOM_SEED()
  do ii=1_ik, sparse_matrix_coo_local%rows_num
   call RANDOM_NUMBER(bmat(ii))
   call RANDOM_NUMBER(rand2)
   if(rand2 .gt. 0.5_rk) then
     bmat(ii) = - bmat(ii)/rand2
   else
     bmat(ii) =  bmat(ii)/rand2
   endif
  end do
  
  
  call MPI_Barrier(mpi_comm_elm,ierr)
  time_start = get_time()
  call cel_elmSolverDoIt(bmat,iierr,niter,ftol)
  time_end = get_time()
  call MPI_Barrier(mpi_comm_elm,ierr)
  if(mpi_rank .eq. 0_mpik) then
    write (cel_output_unit,'(I0,A,E15.6)',advance="yes") mpi_rank, &
     ":solution time:",time_end-time_start
    write (cel_output_unit,'(I0,A,E15.6)',advance="yes") mpi_rank, &
     ":ftol:",ftol
    write (cel_output_unit,'(I0,A,I0)',advance="yes") mpi_rank, &
     ":iter:",niter
  endif
  
  !call cel_elmSolverPerformance()

   call MPI_Barrier(mpi_comm_elm,ierr)
  time_start = get_time()
  call cel_elmSolverCollectSolution(xmat)
  time_end = get_time()
  !if(mpi_rank .eq. 30_ik) then
  ! call print("bmat",bmat)
  ! call print("xmat",xmat)
  !end if
  call MPI_Barrier(mpi_comm_elm,ierr)
  if(mpi_rank .eq. 0_mpik) then
    write (cel_output_unit,'(I0,A,E15.6)',advance="yes") mpi_rank, &
     "collect solustion time:",time_end-time_start
  end if

  do ii=1_ik, 0_ik
    call cel_elmSolverReDoIt(bmat,iierr,niter,ftol)
    call cel_elmSolverCollectSolution(xmat)
    if(mpi_rank .eq. 0) then
     write(*,'(I0,A,I0)')mpi_rank,"solution number:",ii
    end if

  end do
  call cel_elmSolverDestroy()
  deallocate(bmat)
  CALL MPI_FINALIZE(ierr);
  
end program cel_elm

