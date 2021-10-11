
program cel_mpi_omp
!library fbase
use cel_types_module
use cel_input_parameters_module
use cel_omp_module
use cel_timer_interface_module
!library fspmt
use cel_vec_module
use cel_sp_mat_module
use cel_sp_mat_generate_module
!library fsp_mat_distr
use cel_sp_mat_distr_module
use cel_sp_mat_distr_opt_module
use cel_sp_mat_distr_gl_module
use cel_sp_mat_distr_mv_module
!external libraries
use OMP_LIB
use MPI
use, intrinsic :: iso_c_binding

implicit none
  type(cel_input_parameters_type)               :: input_parameters !user input
  integer(kind = mpik)                          :: ierr, mpi_comm, mpi_rank !mpi parameters
  integer(kind = mpik)                          :: mpi_size, mpi_ik_type, ik_mpi !mpi parameters
  logical(kind=lk)                              :: output_on !output parameters
  logical(kind=lk)                              :: to_stop !debuffer parameters
  integer(kind=ik)                              :: mx_cases, max_cases, case_type, test_idx !perf. evaluation
  integer(kind=ik)                              :: where_cross !perf. evaluation
  integer                                       :: scratch_unit=1001 !perf. evaluation
  character(len=80),allocatable, dimension(:)   :: cases !perf. evaluation
  character(len=256)                            :: string, filename !perf. evaluation
  character(len=64)                             :: number_string !perf. evaluation
  integer (kind=ik)                             :: curr_num_threads !threads number
  type(cel_omp_type), dimension(:), allocatable :: omp !parameters for threads; thread num, master or worker thread and so on
  integer(kind=ik)                              :: num_threads, omp_indx !OpenMP / threads control
  type(cel_sp_mat_type)                       :: sparse_matrix_coo !sparse matrix in coordinate format of system to solve 
  integer(kind=ik)                              :: block_indx !block_id of coo matrix
  type(cel_vec_type)                            :: boundary_vector !boundary vector of system to solve
  type(cel_vec_type)                            :: th_solution_vector !theoretical solution of the system
  real(kind=rk)                                 :: time_start, time_end !perf. evaluation
  integer(kind=ik)                              :: rows_num, num_elements, global_num_elements, global_num_rows !matrix parameters
  integer(kind=ik)                              :: nn, err_code !help variable
  integer(kind=c_int)                           :: num_threads_c !matrix properties
  logical(kind=lk)                              :: is_index_based_one, is_to_convert_coo, flag_error !information aboiut init matrix
  type(cel_sp_mat_distr_type)                         :: sp_mat_distr!set of blocked matrix for iteration
  type(cel_raw_blocks_type)                     :: raw_blocks!raw blocks of matrix distribution
  type(cel_hv_blocks_type)                      :: blocks!blocks of matrix distribution
  integer(kind=ik)                              :: temp_1, temp_2, ii
  integer(kind=ik)                              :: num_chunks, chunk_size, block_size, num_raw_blocks, global_flop
  type(cel_blvector_type)                       :: blvector_x
  type(cel_vec_type)                            :: vector_y
  integer(kind=mpik)                            :: mpi_win
  
  num_chunks = 8
  chunk_size = 32
  block_size = num_chunks*chunk_size
  to_stop = .TRUE._lk
  mpi_comm = MPI_COMM_WORLD
  ik_mpi = INT(ik, kind = mpik)
  call MPI_INIT(ierr)
  call cel_error_proc("MPI_INIT", mpi_rank, ierr, output_on, to_stop)
  call MPI_COMM_RANK(mpi_comm, mpi_rank, ierr)
  call cel_error_proc("MPI_COMM_RANK", mpi_rank, ierr, output_on, to_stop)
  call MPI_COMM_SIZE(mpi_comm, mpi_size, ierr)
  call cel_error_proc("MPI_COMM_SIZE", mpi_rank, ierr, output_on, to_stop)
  call MPI_TYPE_MATCH_SIZE(MPI_TYPECLASS_INTEGER, ik_mpi, mpi_ik_type, ierr)
  call cel_error_proc("MPI_TYPE_MATCH_SIZE", mpi_rank, ierr, output_on, to_stop)
  if (mpi_rank .EQ. 0) then
    cel_is_in_debug = .TRUE._lk
  else
    cel_is_in_debug = .FALSE._lk
  endif
  call cel_input_parameters_read(input_parameters)!read user input parameters from command line
  num_chunks = input_parameters%num_chunks
  chunk_size = input_parameters%chunk_size

  if (mpi_rank .EQ. 0) then
    output_on = input_parameters%verbosity
    open(newunit=scratch_unit,file='scratch.dat')!open sratch file on proc 0 for performance evaluation
    !initiliaze structures for performance evaluation
    mx_cases=input_parameters%omp%num_threads-1
    max_cases = mx_cases
    allocate(cases(mx_cases))
    case_type=0
    do nn=1,mx_cases
      write(string,'(I0,A)') nn+1, trim("_threads_CSR blocked / ")
      where_cross=index(string,'#')
      if(where_cross>0) string=string(1:where_cross-1)
      if(trim(string) == '') cycle
      case_type=case_type+1
      cases(case_type)=trim(string)
      max_cases=case_type
    enddo
    filename='result'
    write(number_string,'(a,i0,a,i0)') 'OMP_FROM_2_TO_',input_parameters%omp%num_threads, &
      "_MPI_", mpi_size
    filename=trim(adjustl(filename))//'_'//trim(number_string)
  else
    mx_cases=input_parameters%omp%num_threads-1
    max_cases = mx_cases
    allocate(cases(0))
    output_on = .FALSE._lk
  endif
  curr_num_threads =  input_parameters%omp%num_threads
  call omp_set_dynamic(.FALSE.)
  do curr_num_threads = input_parameters%omp%num_threads, input_parameters%omp%num_threads

    num_threads = INT(curr_num_threads,kind=c_ik)  
    num_threads_c = INT(num_threads,kind=c_ik)  
    call omp_set_num_threads(num_threads_c)
    if(allocated(omp)) deallocate(omp)
    allocate(omp(num_threads))
    omp_indx = cel_master_out_thread+1
    !$omp parallel private(omp_indx, err_code) shared(mpi_rank,global_num_rows,global_num_elements, rows_num, num_elements) 
    omp_indx = omp_get_thread_num() + 1
    call init(omp(omp_indx), num_threads, omp_indx-1,&
      mpi_size, mpi_rank, mpi_comm, input_parameters%verbosity)
    block_indx = cel_index_undef
    if(.not. (omp(omp_indx)%num_workers .EQ. curr_num_threads-1) ) then
        call cel_error("MAIN check number of threads: could not run the request number", &
          &1_ik,.TRUE._lk,.TRUE._lk)
    end if
    test_idx=1
    !do test_idx=1, input_parameters%tests_num
      !initialize matrix for the test_idx
      rows_num = 0_ik
      num_elements = 0_ik

      !$omp barrier

      !make coo matrix with poisson
      if(input_parameters%generate_matrix .EQ. 1) then
        if(omp(omp_indx)%is_master) then
          is_to_convert_coo = .TRUE._lk
          block_indx = omp(omp_indx)%master_num
          call cel_sp_mat_generate(sparse_matrix_coo, &
            boundary_vector, th_solution_vector, &
            input_parameters, omp(omp_indx), test_idx, err_code, &
            block_indx=block_indx)
          is_index_based_one = .TRUE.
          rows_num = sparse_matrix_coo%rows_num
          num_elements = sparse_matrix_coo%elements_num
          call MPI_ALLREDUCE(rows_num, global_num_rows, 1_mpik, mpi_ik_type, MPI_SUM, mpi_comm, ierr)
          call MPI_ALLREDUCE(num_elements, global_num_elements, 1_mpik, mpi_ik_type, MPI_SUM, mpi_comm, ierr)
          call cel_error_proc("MPI_ALLREDUCE", mpi_rank, ierr, output_on, to_stop)
          !$OMP FLUSH (rows_num, num_elements, global_num_elements, global_num_rows)
        endif
      else if(input_parameters%generate_matrix .EQ. 2) then
         if(omp(omp_indx)%is_master) then
          is_to_convert_coo = .TRUE._lk
          block_indx = omp(omp_indx)%master_num
          call cel_sp_mat_generate_simple_test(sparse_matrix_coo, &
            boundary_vector, th_solution_vector, &
            input_parameters, omp(omp_indx), test_idx, err_code)
           is_index_based_one = .TRUE.
          rows_num = sparse_matrix_coo%rows_num
          num_elements = sparse_matrix_coo%elements_num
          call MPI_ALLREDUCE(rows_num, global_num_rows, 1_mpik, mpi_ik_type, MPI_SUM, mpi_comm, ierr)
          call MPI_ALLREDUCE(num_elements, global_num_elements, 1_mpik, mpi_ik_type, MPI_SUM, mpi_comm, ierr)
          call cel_error_proc("MPI_ALLREDUCE", mpi_rank, ierr, output_on, to_stop)
          !$OMP FLUSH (rows_num, num_elements, global_num_elements, global_num_rows)
        end if
      else
        !$omp barrier
        if(omp(omp_indx)%is_master) then
          call cel_error_proc("Error matrix can not be generated; check command line", mpi_rank, 1_mpik, output_on, .FALSE._lk)
          call MPI_FINALIZE(ierr);
          stop
        endif
      end if

      !coo matrix with poisson is done
      !$omp barrier
      if(omp(omp_indx)%is_master .and. omp(omp_indx)%master_num .eq. 0_ik ) then
          write(*,'(9(A,I0))') "test_idx:", &
          test_idx,"; num_threads:",omp_get_num_threads(), "; g_num_rows:", global_num_rows, &
            "; gl_num_elements:", global_num_elements,"; rows_num:",rows_num,&
             "; num_elements:",num_elements,&
             "; num_procs:",omp(omp_indx)%num_masters,&
             "; num_chunks:",num_chunks,"; chunk_size:", chunk_size
          !call MPI_Barrier(omp(omp_indx)%master_comm, ierr)
          !if(omp(omp_indx)%master_num .EQ. 0)  call print("0 B:",boundary_vector%values)
          !call MPI_Barrier(omp(omp_indx)%master_comm, ierr)
          !if(omp(omp_indx)%master_num .EQ. 1)  call print("1 B:",boundary_vector%values)
          !call MPI_Barrier(omp(omp_indx)%master_comm, ierr)
          !if(omp(omp_indx)%master_num .EQ. 0)  call print("0 X ORIG:",th_solution_vector%values)
          !call MPI_Barrier(omp(omp_indx)%master_comm, ierr)
          !if(omp(omp_indx)%master_num .EQ. 1)  call print("1 X ORIG:",th_solution_vector%values)
          !call MPI_Barrier(omp(omp_indx)%master_comm, ierr)
      end if
      !divide the generated coo matrix in small sub matricies of size (num_chunks * chunk_size)
      !also initialize the h_blocks and v_blocks, these are used by remote access to the vector data
      call cel_sp_mat_distr_csr_from_coo(sp_mat_distr, raw_blocks, blocks,&
        num_chunks,&!number of chunks (eq. to the size of one remote access data) in sub matrix in one direction
        chunk_size,&!number of vector elements in a chunk (eq. to the size of one remote access data)
        rows_num,&!local number of rows of the matrix
        rows_num,&!local number of cols of the matrix
        num_elements,&!local number of elements of the matrix
        global_num_rows,&!global number of rows of the matrix
        sparse_matrix_coo%block_rows_num,&! 1 
        sparse_matrix_coo%block_cols_num,&! 1
        sparse_matrix_coo%rows,&!original matrix in coo format
        sparse_matrix_coo%cols,&!original matrix in coo format
        sparse_matrix_coo%values,&!original matrix in coo format
        omp(omp_indx),&!topology of threads and processes
        3_ik,&!place in l1 for 3 addional vectors per raw block
        output_on,&
        err_code)
        call cel_error_proc("cel_sp_mat_distr_csr_from_coo", mpi_rank, INT(err_code,kind=mpik), output_on, .TRUE._lk)
        num_raw_blocks = raw_blocks%num_raw_blocks
      if(omp(omp_indx)%is_master) then
        call new(vector_y,chunk_size*num_chunks*blocks%num_v_blocks)
      end if
      !$omp barrier

      call cel_blvector_new(blvector_x, blocks , num_chunks, chunk_size, omp(omp_indx))
      if(omp(omp_indx)%is_master) then
        sp_mat_distr%flop=0_ik
        sp_mat_distr%mv_time_sec=0_rk
      end if
      !$omp barrier
      mpi_win = 0_mpik
      time_start = get_time()
      do ii=1, input_parameters%attempts_num 
        call cel_sp_mat_distr_mv(vector_y%values, th_solution_vector%values, blvector_x, sp_mat_distr, num_raw_blocks, &
          blocks,  omp(omp_indx), num_chunks, chunk_size, mpi_win)
        !$omp barrier
  !      if(omp(omp_indx)%is_master) then
  !        call MPI_BARRIER(mpi_comm, ierr); call cel_error_proc("MPI_BARRIER", mpi_rank, ierr, output_on, to_stop)
  !      end if
      end do
      time_end  = get_time()
      if(omp(omp_indx)%is_master) then
        call MPI_ALLREDUCE(sp_mat_distr%flop, global_flop, 1_mpik, mpi_ik_type, MPI_SUM, mpi_comm, ierr)
        call MPI_BARRIER(mpi_comm, ierr); call cel_error_proc("MPI_BARRIER", mpi_rank, ierr, output_on, to_stop)
      end if
      if(omp(omp_indx)%is_master) then
          write (*,'(A,I0,5(A,E15.6))') "; proc:",omp(omp_indx)%master_num, &
          "; time:",time_end-time_start,&
          "; mv time:",sp_mat_distr%mv_time_sec/omp(omp_indx)%num_workers,&
           "; flops:", REAL(global_flop,kind=rk)/(time_end-time_start),&
           "; my flops:", REAL(sp_mat_distr%flop,kind=rk)/(time_end-time_start),&
           "; mv flops:", REAL(sp_mat_distr%flop,kind=rk)/(sp_mat_distr%mv_time_sec/omp(omp_indx)%num_workers)
           
      end if
     
      if(omp(omp_indx)%is_master) then
        flag_error = .false._lk
        do ii=1, rows_num
          if(ABS(vector_y%values(ii)-boundary_vector%values(ii)) .GE. 0.0000001) then
            flag_error = .true._lk
            !write (*,'(I0,A,I0,A,E16.6,A,E16.6)') omp(omp_indx)%master_num,&
            ! ": error in ii:",ii,"; calc:", vector_y%values(ii), "; th:",boundary_vector%values(ii)
          end if
        end do
        if(flag_error) write (*,*) "MV with ERROR"
      end if
      if(omp(omp_indx)%is_master) then
     
!        call MPI_BARRIER(mpi_comm, ierr); call cel_error_proc("MPI_BARRIER", mpi_rank, ierr, output_on, to_stop)

        call del(sparse_matrix_coo)
        call del(boundary_vector)
        call del(th_solution_vector)
        call del(vector_y)
        call del(blvector_x)
        call del(sp_mat_distr)
        call del(raw_blocks)
        call del(blocks)
      endif
     
   ! end do
    !$omp end parallel
    call MPI_Win_Free(mpi_win,ierr)
  end do

  CALL MPI_FINALIZE(ierr);
  call cel_error_proc("MPI_COMM_SIZE", mpi_rank, ierr, output_on, to_stop)
end program cel_mpi_omp

