!interface for elmfire
!Task 1. initializing of cel obejects:
!cel_elmSolverInitialization - 
!  a. here we have to decide about topology of the linear solver
!  b. define the mapping between the linear solver mpi processes (master processes in cel) and 
!     elmfire mpi processes (in comm) 
!  b. initialize omp_threads array(contained mpi ranks as master_num, omp thread pids and 
module cel_elmsolver_module


use cel_types_module
use cel_sp_mat_module
use cel_input_parameters_module
use cel_perf_module
use cel_cpu_param_module
use cel_comm_module
use cel_sp_mat_distr_module
use cel_omp_shared_work_module
use cel_sp_mat_distr_vec_module
use cel_sp_mat_distr_gl_module
use cel_perf_distr_module
use cel_cgalg_module
use cel_gmresalg_module
use cel_profile_module
use cel_sp_mat_distr_format_module
use cel_prec_module
use MPI
use OMP_LIB
implicit none



private

integer :: elm_nzero(2),elm_thisRow,elm_totSize,elm_blckSize,elm_blckLow
integer(kind=mpik) :: elm_comm
integer(kind=mpik) :: elm_reduce_comm
integer(kind=ik), dimension(:), allocatable :: elm_local_rows ! row index array of local matrix part on each mpi
integer(kind=ik), dimension(:), allocatable :: elm_local_cols ! col index array of local matrix part on each mpi
real(kind=rk), dimension(:), allocatable    :: elm_local_values ! value array of local matrix part on each mpi
real(kind=rk), dimension(:), allocatable    :: elm_local_diag_idx ! value array of local matrix diagonal part on each mpi

integer(kind=ik) :: elm_num_local_elements !number of matrix elements on of local matrix part on each mpi
integer(kind=ik) :: elm_num_local_rows!number of matrix rows on of local matrix part on each mpi

integer(kind=mpik) :: solver_mpi_comm!MPI communicator for cel
type(cel_sp_mat_type)  :: sparse_matrix_coo_local!reduced matrix to the master procs
type(cel_sp_mat_type), dimension(:), allocatable  :: solver_sparse_matrix_a !sparse matrix to solve in cel format
type(cel_comm_type)                               :: solver_comm 
type(cel_sp_mat_distr_vec_type)                   :: vector_x
integer(kind=ik), dimension(:), allocatable   :: vtxdist ! array of matrix distribution (i-th el. is the index of the first row on i-th process)

type(cel_profile_type),dimension(:), allocatable               :: perf_profiles
type(cel_perf_counter_type)                                    :: perf_counter_mv
type(cel_perf_counter_type)                                    :: perf_counter_vec
type(cel_perf_counter_type)                                    :: solver_local_perf
type(cel_perf_counter_type)                                    :: solver_total_perf

integer(kind=mpik), allocatable, dimension(:) :: receive_count_rows
integer(kind=mpik), allocatable, dimension(:) :: receive_displ_rows
integer(kind=ik) :: solver_num_local_rows, num_reduce_threads

type(cel_gmresalg_parameter_type)                              :: solver_gmres
type(cel_cgalg_parameter_type)                                 :: solver_cg
type(cel_sp_mat_distr_vec_type)                                :: vector_r
type(cel_sp_mat_distr_vec_type)                                :: vector_d
type(cel_sp_mat_distr_vec_type)                                :: vector_q
type(cel_sp_mat_distr_vec_type)                                :: vector_v
type(cel_sp_mat_distr_vec_type),dimension(:), allocatable      :: gmres_vector_v
type(cel_sp_mat_distr_vec_type),dimension(:), allocatable      :: gmres_vector_z
type(cel_sp_mat_distr_vec_type),dimension(:), allocatable      :: gmres_vector_w
type(cel_sp_mat_distr_vec_type)                                :: vector_tmp
type(cel_sp_mat_distr_vec_type)                                :: boundary_vector
real(kind=rk), allocatable, dimension(:,:)        :: gmres_hh
real(kind=rk), allocatable, dimension(:)          :: gmres_cc
real(kind=rk), allocatable, dimension(:)          :: gmres_ss
real(kind=rk), allocatable, dimension(:)          :: gmres_gamma
real(kind=rk), allocatable, dimension(:)          :: gmres_alpha
integer(kind=mpik) :: mpi_type_rk
integer(kind=mpik) :: mpi_type_ik
integer(kind=mpik) :: mpi_type_perf_ik
type(cel_omp_type), dimension(:), allocatable :: solver_omps !parameters for threads; thread num, master or worker thread and so on
type(cel_omp_shared_work_type)                :: solver_omp_shared_work !manages the number of the worker threads, which are involved in computational tasks

public cel_elmSolverInitialization, cel_elmSolverDatarowStart, &
 cel_elmSolverDataInsertion, cel_elmSolverDataReinsertion, &
 cel_elmSolverDatarowCompleted, cel_elmSolverMatrixCompleted, &
 cel_elmSolverDoIt, cel_elmSolverReDoIt, cel_elmSolverCollectSolution, &
 cel_elmSolverDestroy,cel_elmSolverPerformance



contains
!ToDo: initializing of common objects:
!omp
!input parameters:
!in itotSize - number of global rows
!in ibs - number of local rows
!in ibl - first local index of rows
!in nnzero(1) - number of nonzeros per row in DIAGONAL portion of local submatrix (same value is used for all local rows) 
!in nnzero(2) - number of nonzeros per row in the OFF-DIAGONAL portion of local submatrix (same value is used for all local rows)
!in comm -MPI communicator 
subroutine cel_elmSolverInitialization(itotSize,ibs,ibl,nnzero,comm)
    integer(kind=ik),intent(IN) :: itotSize,ibs,ibl,nnzero(2)
    integer(kind=mpik),intent(IN) :: comm
    integer(kind=ik) :: array_size
    integer(kind=ik) :: sys_stat
    integer(kind = mpik)      :: mpi_rank, ierr

   elm_comm = comm
   elm_blckLow = ibl
   elm_num_local_rows = ibs
   array_size = (nnzero(1)+nnzero(2))*ibs

   call MPI_COMM_RANK(elm_comm, mpi_rank, ierr)
   if(mpi_rank .eq. 0_ik) then
     input_parameters%verbosity = input_parameters%verbosity
   else
     input_parameters%verbosity = input_parameters%verbosity
   endif
   allocate(elm_local_rows(array_size), stat=sys_stat)       
   call cel_error("Error cel_elmSolverInitialization not enought memory", sys_stat,&
                        input_parameters%verbosity,.TRUE._lk)
   allocate(elm_local_cols(array_size), stat=sys_stat)       
   call cel_error("Error cel_elmSolverInitialization not enought memory", sys_stat,&
                        input_parameters%verbosity,.TRUE._lk)
   allocate(elm_local_values(array_size), stat=sys_stat)       
   call cel_error("Error cel_elmSolverInitialization not enought memory", sys_stat,&
                        input_parameters%verbosity,.TRUE._lk)
   allocate(elm_local_diag_idx(ibs), stat=sys_stat)       
   call cel_error("Error cel_elmSolverInitialization not enought memory", sys_stat,&
                        input_parameters%verbosity,.TRUE._lk)               
  elm_num_local_elements = 0_ik
  input_parameters%cel_omp%num_threads = input_parameters%num_reduce_threads
  input_parameters%is_matrix_symmetric = .FALSE._lk
  input_parameters%block_rows_num = 1_ik
  input_parameters%block_cols_num = 1_ik
  input_parameters%write_performance = .FALSE._lk
  input_parameters%write_profile = .FALSE._lk
  input_parameters%iter_max = 5000_ik
  input_parameters%eps = 0.1e-07_rk
  input_parameters%mv_algorithm = MV_CSR
  input_parameters%mv_algorithm_off_diag = MV_COO
  input_parameters%mv_communicator = MV_COMM_ISEND_COMM_BUFFER
  input_parameters%cpu_freq_idx = 0
  input_parameters%benchmark_id  = BENCHMARK_GMRES
  input_parameters%matrix_check_null_diagonal = .FALSE._lk
  input_parameters%preconditioner = PRECONDITIONER_JACOBI
  input_parameters%diagonal_block_size = 1_ik
  input_parameters%matrix_scaler  = MATRIX_SCALER_DIAGONAL
  input_parameters%gmres_max_restart  = 4_ik - 1_ik
  ! input_parameters%num_reduce_threads = cel_cpu_param_cores
end subroutine

subroutine cel_elmSolverDatarowStart(row)
integer(kind=ik),intent(IN) :: row
 elm_thisRow = row
 return
end subroutine

subroutine cel_elmSolverDataInsertion(col,value)
integer,intent(IN) :: col
real(kind=rk),intent(IN) :: value

  elm_num_local_elements = elm_num_local_elements + 1_ik
  elm_local_rows(elm_num_local_elements) = elm_thisRow
  elm_local_cols(elm_num_local_elements) = col
  elm_local_values(elm_num_local_elements) = value
  if(col .eq. elm_thisRow) then
    elm_local_diag_idx(elm_thisRow-elm_blckLow+1_ik) = elm_num_local_elements
  endif  
end subroutine

subroutine cel_elmSolverDataReinsertion(col,value)
integer(kind=ik),intent(IN) :: col
real(kind=rk),intent(IN) :: value
integer(kind=ik) :: idx
  
   if(col .eq. elm_thisRow) then
    idx = elm_local_diag_idx(elm_thisRow-elm_blckLow+1_ik)
    elm_local_values(idx) = value
  else
   call cel_error("Error cel_elmSolverDataReinsertion not implemented for off diagonal elements", 1_ik,&
                 input_parameters%verbosity,.TRUE._lk)
  endif
  return
end subroutine

subroutine cel_elmSolverDatarowCompleted
 return
end subroutine

subroutine cel_elmSolverMatrixCompleted
  integer(kind = mpik)      :: mpi_size, ik_mpi, mpi_new_size !mpi parameters
  integer(kind = mpik)      :: rk_mpi, perf_ik_mpi
  integer(kind = mpik) ::  mpi_rank_red, solver_mpi_rank, num_mpik,mpi_rank
  integer(kind=mpik) :: num_reduce_threads_mpik, num_masters_mpik
  integer(kind=ik) :: rank, master_num, count, ii, sys_stat, num_masters
  integer(kind = mpik), allocatable, dimension(:) :: solver_ranks_incl, reduce_ranks_incl
  integer(kind=mpik) :: elmf_group, solver_group, reduce_group
  integer(kind=mpik)  :: ierr, master_num_mpik, tag, source
  integer(kind=mpik), allocatable, dimension(:) :: receive_count_el
  integer(kind=mpik), allocatable, dimension(:) :: receive_displ_el
  integer(kind=ik) :: solver_num_local_elements,idx1, idx2
  integer(kind=ik) ::  thread_num
  integer(kind=mpik) :: elm_num_local_elements_mpik, elm_num_local_rows_mpik
  integer(kind=ik), dimension(:), allocatable   :: procs_array! array with proc_id for the communication elements in sp_mat_local_array
  integer(kind=ik), dimension(:), allocatable   :: procs_array_first_idx ! array with first idx of data needed from the process in procs_array
  integer(kind=ik), dimension(:), allocatable   :: procs_array_last_idx  ! array with last idx of data needed from the process in procs_array
  integer(kind=ik), dimension(:,:), allocatable :: com_buffer_perm
  integer(kind=ik), dimension(:), allocatable   :: length_comm_buffer
  integer(kind=ik)                              :: omp_indx, err_code
  integer(kind=ik), dimension(0) :: dummy_array_int
  real(kind=rk), dimension(0) :: dummy_array_real
  type(cel_perf_counter_type)                   :: solver_local_perf
  type(cel_omp_shared_work_type)                :: solver_omp_shared_work !manages the number of the worker threads, which are involved in computational tasks
  logical(kind=lk)                              :: to_scale

  ik_mpi = INT(ik, kind = mpik)
  rk_mpi = INT(rk, kind = mpik)
  perf_ik_mpi = INT(perf_ik, kind = mpik)
  
  call MPI_COMM_RANK(elm_comm, mpi_rank, ierr)
  call MPI_COMM_SIZE(elm_comm, mpi_size, ierr)
  call MPI_TYPE_MATCH_SIZE(MPI_TYPECLASS_INTEGER, ik_mpi, mpi_type_ik, ierr)
  call MPI_TYPE_MATCH_SIZE(MPI_TYPECLASS_INTEGER, perf_ik_mpi, mpi_type_perf_ik, ierr)
  call MPI_TYPE_MATCH_SIZE(MPI_TYPECLASS_REAL, rk_mpi, mpi_type_rk, ierr)
 

 
  num_reduce_threads = input_parameters%num_reduce_threads
  num_reduce_threads_mpik = int(num_reduce_threads,kind=mpik)
  num_masters = mpi_size / num_reduce_threads
  num_masters_mpik = int(num_masters,kind=mpik)

  if(num_masters*num_reduce_threads .ne. mpi_size) then
  call cel_error(&
               "Error cel_elmSolverMatrixCompleted: number of mpi procs  must always be divisible"//&
               " by number of reduce_threads",&
                1_ik, input_parameters%verbosity,.TRUE._lk)
  end if
  rank = int(mpi_rank,kind=ik)
  thread_num = modulo(mpi_rank,num_reduce_threads)
  master_num = rank / num_reduce_threads
  master_num_mpik = int(master_num,kind=mpik)
  !create new communicator for the solver
  call MPI_Comm_group(elm_comm,elmf_group,ierr)
  allocate(solver_ranks_incl(num_masters))
  count=0_ik
  do ii=0_ik,num_masters-1_ik
    count=count+1_ik
    solver_ranks_incl(count)=int(ii*num_reduce_threads,kind=mpik)
  end do
  call MPI_Group_incl(elmf_group,num_masters_mpik,solver_ranks_incl(1_ik:num_masters),&
  solver_group,ierr)
  call cel_error(&
             "Error cel_elmSolverMatrixCompleted: could not create solver group!", &
              ierr, input_parameters%verbosity,.TRUE._lk)
  call MPI_Comm_create(elm_comm, solver_group,solver_mpi_comm,ierr)
  call cel_error(&
             "Error cel_elmSolverMatrixCompleted: could not create solver communicator!", &
              ierr, input_parameters%verbosity,.TRUE._lk)
  call MPI_Group_free(solver_group,ierr)
  if(thread_num .eq. 0_ik) then
   call MPI_COMM_RANK(solver_mpi_comm, solver_mpi_rank, ierr)
  end if

  !create new communicator for the data reorganization
  allocate(reduce_ranks_incl(num_reduce_threads))
  count=0
  do ii=0_ik,num_reduce_threads-1_ik
    count=count+1_ik
    reduce_ranks_incl(count)=int(master_num*num_reduce_threads+ii,kind=mpik)
  end do

  call MPI_Group_incl(elmf_group,num_reduce_threads_mpik,reduce_ranks_incl(1_ik:num_reduce_threads),& 
    reduce_group,ierr)
  call cel_error(&
               "Error cel_elmSolverMatrixCompleted: could not create reduce group!", &
                ierr, input_parameters%verbosity,.TRUE._lk)
  call MPI_Comm_create(elm_comm, reduce_group,elm_reduce_comm,ierr)
  call cel_error(&
               "Error cel_elmSolverMatrixCompleted: could not create reduce communicator!", &
                ierr, input_parameters%verbosity,.TRUE._lk)
  call MPI_Group_free(reduce_group,ierr)
  call MPI_COMM_RANK(elm_reduce_comm, mpi_rank_red, ierr)
  call MPI_Group_free(elmf_group,ierr)
  
  elm_num_local_elements_mpik = int(elm_num_local_elements,kind=mpik)
  elm_num_local_rows_mpik = int(elm_num_local_rows,kind=mpik)
  if(mpi_rank_red .eq. 0_ik) then
    allocate(receive_count_el(num_reduce_threads))
    allocate(receive_displ_el(num_reduce_threads))
    allocate(receive_count_rows(num_reduce_threads))
    allocate(receive_displ_rows(num_reduce_threads))
    receive_count_el(1:num_reduce_threads) = 0_mpik
    receive_displ_el(1:num_reduce_threads) = 0_mpik
    receive_count_rows(1:num_reduce_threads) = 0_mpik
    call MPI_Gather(elm_num_local_elements_mpik,1_mpik,MPI_INTEGER,&
     receive_count_el(1_ik:num_reduce_threads), 1_mpik, MPI_INTEGER,0_mpik,elm_reduce_comm,ierr) 
    call MPI_Gather(elm_num_local_rows_mpik,1_mpik,MPI_INTEGER,&
     receive_count_rows(1_ik:num_reduce_threads), 1_mpik, MPI_INTEGER,0_mpik,elm_reduce_comm,ierr) 
  else
    allocate(receive_count_el(0_ik))
    allocate(receive_displ_el(0_ik))
    allocate(receive_count_rows(0_ik))
    allocate(receive_displ_rows(0_ik))
    call MPI_Gather(elm_num_local_elements_mpik,1_mpik,MPI_INTEGER,&
     receive_count_el, 1_mpik, MPI_INTEGER,0_mpik,elm_reduce_comm,ierr) 
    call MPI_Gather(elm_num_local_rows_mpik,1_mpik,MPI_INTEGER,&
     receive_count_rows, 1_mpik, MPI_INTEGER,0_mpik,elm_reduce_comm,ierr) 
  endif
  !Collect the matrix from all mpi processes to the master processes
  if(mpi_rank_red .eq. 0_ik) then
    solver_num_local_elements = sum(receive_count_el(1_ik:num_reduce_threads))
    solver_num_local_rows = sum(receive_count_rows(1_ik:num_reduce_threads))
    receive_displ_el(1_ik)=1_mpik
    receive_displ_rows(1_ik)=1_mpik
    do ii=2_ik,num_reduce_threads
      receive_displ_el(ii)=receive_displ_el(ii-1_ik)+receive_count_el(ii-1_ik)
      receive_displ_rows(ii)=receive_displ_rows(ii-1_ik)+receive_count_rows(ii-1_ik)
    end do
    receive_displ_el(1_ik:num_reduce_threads) = receive_displ_el(1_ik:num_reduce_threads) - 1_mpik
    receive_displ_rows(1_ik:num_reduce_threads) = receive_displ_rows(1_ik:num_reduce_threads) - 1_mpik
    call cel_sp_mat_new(sparse_matrix_coo_local, solver_num_local_rows, solver_num_local_rows,&
     input_parameters%block_rows_num , input_parameters%block_cols_num ,solver_num_local_elements)
    call MPI_Gatherv(elm_local_cols(1_ik:elm_num_local_elements), elm_num_local_elements_mpik, &
      mpi_type_ik, sparse_matrix_coo_local%cols(1_ik:solver_num_local_elements), receive_count_el(1_ik:num_reduce_threads), &
     receive_displ_el(1_ik:num_reduce_threads),mpi_type_ik,0_mpik,elm_reduce_comm,ierr)
    call MPI_Gatherv(elm_local_rows(1_ik:elm_num_local_elements), elm_num_local_elements_mpik, &
      mpi_type_ik, sparse_matrix_coo_local%rows(1_ik:solver_num_local_elements), receive_count_el(1_ik:num_reduce_threads), &
      receive_displ_el(1_ik:num_reduce_threads),mpi_type_ik,0_mpik,elm_reduce_comm,ierr)
    call MPI_Gatherv(elm_local_values(1_ik:elm_num_local_elements), elm_num_local_elements_mpik, &
      mpi_type_rk, sparse_matrix_coo_local%values(1_ik:solver_num_local_elements), receive_count_el(1_ik:num_reduce_threads), &
      receive_displ_el(1_ik:num_reduce_threads),mpi_type_rk,0_mpik,elm_reduce_comm,ierr)
  else
    call MPI_Gatherv(elm_local_cols(1_ik:elm_num_local_elements), elm_num_local_elements_mpik, &
      mpi_type_ik, dummy_array_int, receive_count_el, &
      receive_displ_el,mpi_type_ik,0_mpik,elm_reduce_comm,ierr)
    call MPI_Gatherv(elm_local_rows(1_ik:elm_num_local_elements), elm_num_local_elements_mpik, &
      mpi_type_ik, dummy_array_int, receive_count_el, &
      receive_displ_el,mpi_type_ik,0_mpik,elm_reduce_comm,ierr)
    call MPI_Gatherv(elm_local_values(1_ik:elm_num_local_elements), elm_num_local_elements_mpik, &
      mpi_type_rk, dummy_array_real, receive_count_el, &
      receive_displ_el,mpi_type_rk,0_mpik,elm_reduce_comm,ierr)
  end if

  !only the master procs convert the matrix
  if(mpi_rank_red .eq. 0_ik) then
    allocate(solver_omps(1_ik))
    call init(solver_omps(1_ik), 1_ik, 0_ik,&
      num_masters, master_num, solver_mpi_comm, mpi_type_ik, mpi_type_rk, mpi_type_perf_ik,&
      input_parameters%verbosity)
    solver_omps(1_ik)%is_master = .TRUE._lk
    solver_omps(1_ik)%master_comm = solver_mpi_comm
    solver_omps(1_ik)%worker_num = 1_ik
    solver_omps(1_ik)%num_workers = 1_ik
    call cel_omp_shared_work_new(solver_omp_shared_work,solver_omps(1_ik))
    call reset(solver_local_perf)
    allocate(vtxdist(num_masters+1_ik))
    vtxdist(1)=1_ik
    call MPI_ALLGATHER(solver_num_local_rows,1_mpik,solver_omps(1_ik)%mpi_type_ik,&
      vtxdist(2:size(vtxdist)),1_mpik,solver_omps(1_ik)%mpi_type_ik,&
      solver_omps(1_ik)%master_comm,ierr)
    call cel_error_proc("Error in cel_elmSolverMatrixCompleted MPI_ALLREDUCE", solver_mpi_rank,&
      ierr, input_parameters%verbosity, .TRUE._lk)
    do ii=1,solver_omps(1_ik)%num_masters
      vtxdist(ii+1_ik)=vtxdist(ii+1_ik)+vtxdist(ii)
    end do
   !calculate who is the neighbours processes (involved in the communication)
   !and convert matrix in csr format
   call cel_sp_mat_format_sparse_31_sort(sparse_matrix_coo_local, .FALSE._lk, err_code)
   call cel_error_proc("Error in cel_elmSolverMatrixCompleted MPI_ALLREDUCE", solver_mpi_rank,&
      ierr, input_parameters%verbosity, .TRUE._lk)
           !   allocate(th_solution_vector%values(sparse_matrix_coo_global%rows_num))
    call cel_sp_mat_distr_from_sp_mat_local(solver_sparse_matrix_a, procs_array, &
         procs_array_first_idx, procs_array_last_idx,&
         com_buffer_perm,length_comm_buffer,&
         vtxdist, sparse_matrix_coo_local,&
         solver_omps(1_ik), input_parameters, err_code)
    call cel_error("Error in main: cel_elmSolverMatrixCompleted!", &
         &err_code,input_parameters%verbosity,.TRUE._lk)
    call del(sparse_matrix_coo_local)
    call cel_comm_get_connections(solver_comm, procs_array,&
         procs_array_first_idx, procs_array_last_idx,&
         com_buffer_perm,length_comm_buffer,&
         solver_omps(1_ik), err_code)
    call cel_error("Error in main: cel_comm_get_connections!", &
         &err_code,cel_is_in_debug,.FALSE._lk)
    call cel_comm_init_vtxdist_counts(solver_comm, vtxdist, solver_omps(1_ik), err_code)
    call cel_error("Error in main: cel_comm_init_vtxdist_counts!", &
         &err_code,cel_is_in_debug,.FALSE._lk)
   !create the communicators for the halo data transfer within communicators
    if(input_parameters%mv_communicator .eq. MV_COMM_ISEND_GROUPS .or. &
     input_parameters%mv_communicator .eq. MV_COMM_IBCAST) then
      call cel_comm_get_communicators(solver_comm, solver_omps(1_ik), err_code)
      call cel_error("Error in main: cel_comm_get_communicators!", &
         &err_code,cel_is_in_debug,.FALSE._lk)
    else
        allocate(solver_comm%get_comms(0))
    endif
    
    call cel_sp_mat_distr_set_mv_parameters(solver_sparse_matrix_a,solver_omps(1_ik),&
                                 input_parameters,err_code)
    call MPI_BARRIER(solver_omps(1_ik)%master_comm, ierr)
    if(input_parameters%matrix_scaler .eq. MATRIX_SCALER_DIAGONAL) then
       call cel_sp_mat_distr_scale(solver_sparse_matrix_a,&
          solver_omps(1_ik),solver_omp_shared_work,solver_local_perf, input_parameters%verbosity, err_code)
       call MPI_BARRIER(solver_omps(1_ik)%master_comm, ierr)
    end if
    call cel_sp_mat_distr_format(solver_sparse_matrix_a,vtxdist,solver_omp_shared_work,&
      solver_omps(1_ik),solver_local_perf,input_parameters%verbosity,err_code)
    call MPI_BARRIER(solver_omps(1_ik)%master_comm, ierr)
    if(input_parameters%preconditioner .eq. PRECONDITIONER_JACOBI) then
       to_scale = input_parameters%matrix_scaler .eq. MATRIX_SCALER_DIAGONAL
       call cel_sp_mat_distr_set_jacobi(solver_sparse_matrix_a,to_scale,&
                                solver_omps(1_ik),solver_omp_shared_work,&
                                solver_local_perf, input_parameters%verbosity,err_code)
     end if
     call MPI_BARRIER(solver_omps(1_ik)%master_comm, ierr)
   
     deallocate(solver_omps)
     call del(solver_omp_shared_work)
  endif
 

  call del(sparse_matrix_coo_local)
  deallocate(elm_local_rows) 
  deallocate(elm_local_cols) 
  deallocate(elm_local_values) 
  deallocate(elm_local_diag_idx)
  deallocate(solver_ranks_incl)
  deallocate(reduce_ranks_incl)
  deallocate(receive_count_el)
  deallocate(receive_displ_el)
  

 return
end subroutine cel_elmSolverMatrixCompleted

subroutine cel_elmSolverDoIt(bmat,iierr,niter,ftol)
  real(kind=rk),intent(IN) :: bmat(:)
  integer(kind=ik),intent(OUT) :: iierr,niter
  real(kind=rk),intent(OUT) :: ftol

  integer(kind = mpik)      :: mpi_size,  mpi_new_size !mpi parameters
  integer(kind = mpik) ::  mpi_rank_red, solver_mpi_rank, num_mpik,mpi_rank, ierr
  integer(kind=ik) :: thread_num, master_num, num_masters, ii
  integer(kind=mpik) :: target_mpik
  integer(kind=mpik),dimension(:), allocatable :: send_sync_request
  real(kind=rk), allocatable, dimension(:)          :: reduced_boundary_vect
  integer(kind=mpik) :: get_sync_request
  real(kind=rk),dimension(3_ik) :: sync_out_buff
  logical(kind=lk) :: sync_complited
  real(kind=rk), dimension(0_ik) :: dummy_array_real
  integer(kind=mpik) , dimension(0_ik) :: dummy_array_int
  integer(kind=ik) :: omp_indx,  err_code
  integer(kind=mpik) :: solver_mpi_size, elm_num_local_rows_mpik
  integer(kind=mpik) :: solver_num_local_rows_mpik

  
  call MPI_COMM_RANK(elm_comm, mpi_rank, ierr)
  call MPI_COMM_SIZE(elm_comm, mpi_size, ierr)
  call MPI_COMM_RANK(elm_reduce_comm, mpi_rank_red, ierr)

  
  thread_num = modulo(mpi_rank,num_reduce_threads)
  master_num = mpi_rank / num_reduce_threads  
  num_reduce_threads = input_parameters%num_reduce_threads
  num_masters = mpi_size / num_reduce_threads

 !collect boundary vector
  elm_num_local_rows_mpik = int(elm_num_local_rows,kind=mpik)
  solver_num_local_rows_mpik = int(solver_num_local_rows,kind=mpik)
  if(mpi_rank_red .eq. 0_ik) then
    allocate(reduced_boundary_vect(solver_num_local_rows))
    call MPI_Gatherv(bmat(1_ik:elm_num_local_rows), elm_num_local_rows_mpik, &
      mpi_type_rk, reduced_boundary_vect(1_ik:solver_num_local_rows),&
      receive_count_rows(1_ik:num_reduce_threads), receive_displ_rows(1_ik:num_reduce_threads),&
      mpi_type_rk,0_mpik,elm_reduce_comm,ierr)
  else
     allocate(reduced_boundary_vect(0_ik))
     call MPI_Gatherv(bmat(1_ik:elm_num_local_rows), elm_num_local_rows_mpik, &
      mpi_type_rk, dummy_array_real,&
      dummy_array_int, dummy_array_int,&
      mpi_type_rk,0_mpik,elm_reduce_comm,ierr)
  end if

  if(mpi_rank_red .eq. 0_mpik) then
    iierr = 0_ik
    niter = 0_ik
    ftol = 0_rk
    call MPI_COMM_RANK(solver_mpi_comm, solver_mpi_rank, ierr)
    call MPI_COMM_SIZE(solver_mpi_comm, solver_mpi_size, ierr)
   !start omp threads
   call omp_set_dynamic(.false._lk)
   call omp_set_nested(.false._lk)
   call omp_set_num_threads(num_reduce_threads)
   allocate(solver_omps(num_reduce_threads))
   !$omp parallel default(none) private(ii,err_code,omp_indx)&
    !$omp& shared(solver_sparse_matrix_a, vtxdist,boundary_vector,vector_x,&
      !$omp& vector_r,vector_d,vector_q,vector_v,vector_tmp,&
      !$omp& gmres_hh,gmres_cc,gmres_ss,gmres_gamma,gmres_alpha,perf_counter_mv,&
      !$omp& solver_local_perf,perf_counter_vec,input_parameters,solver_omps,&
      !$omp& solver_mpi_comm,solver_mpi_rank, mpi_type_ik, mpi_type_rk, mpi_type_perf_ik,&
      !$omp& solver_comm,gmres_vector_v,gmres_vector_z,gmres_vector_w,&
      !$omp& num_reduce_threads,solver_mpi_size, solver_cg, solver_gmres,&
      !$omp& solver_omp_shared_work,perf_profiles,reduced_boundary_vect,niter,ftol,iierr,mpi_rank )

     omp_indx = omp_get_thread_num() + 1
     !initialize thread and processes topology solver_omps
     call init(solver_omps(omp_indx), num_reduce_threads, omp_indx-1,&
       solver_mpi_size, solver_mpi_rank, solver_mpi_comm, &
       mpi_type_ik, mpi_type_rk, mpi_type_perf_ik, input_parameters%verbosity)
     call new(solver_omp_shared_work,solver_omps(omp_indx))
     !initialize vectors
     if(input_parameters%benchmark_id .eq. BENCHMARK_CG) then
      call cel_cgalg_parameter_new(solver_cg,vtxdist,&
        vector_x,vector_r,vector_d,vector_q,vector_v,vector_tmp,&
        solver_local_perf,perf_counter_mv, perf_counter_vec,perf_profiles,&
        solver_omps(omp_indx),solver_comm, input_parameters,err_code)
     elseif(input_parameters%benchmark_id .eq. BENCHMARK_GMRES) then
      call cel_gmresalg_parameter_new(solver_gmres,vtxdist,&
        vector_x,vector_r,gmres_vector_v,gmres_vector_z,gmres_vector_w,vector_tmp,&
        gmres_hh,gmres_cc,gmres_ss,gmres_gamma,gmres_alpha,&
        solver_local_perf,perf_counter_mv,perf_counter_vec,perf_profiles,&
        solver_omps(omp_indx),solver_comm, input_parameters,err_code)
     end if
    !!make diagonal -ToDo
    !if(input_parameters%preconditioner .eq. PRECONDITIONER_JACOBI) then
    !  call cel_prec_build_diagonal(solver_sparse_matrix_a, &
    !        input_parameters%diagonal_block_size,solver_omp_shared_work,&
    !        solver_omps(omp_indx), input_parameters,err_code)
    !  !$omp barrier
    !  call cel_prec_build_diagonal_inv(solver_sparse_matrix_a, &
    !        solver_omp_shared_work,&
    !        solver_omps(omp_indx), input_parameters,err_code)
    !  !$omp barrier
    !  if(err_code .gt. 0) then
    !    !$omp critical
    !      if(input_parameters%preconditioner .ne. PRECONDITIONER_JACOBI) then
    !         input_parameters%preconditioner = PRECONDITIONER_NULL
    !         call cel_error("Warning in cel_elmSolverDoIt  block jacobi can not be applied !", &
    !                     &1_ik,input_parameters%verbosity,.FALSE._lk) 
    !      end if
    !    !$omp end critical
    !  end if
    !end if
    !!if(mpi_rank .eq. 0_ik) then
    !!  if(solver_omps(omp_indx)%is_master) then
    !!   call print("diagnonal",solver_sparse_matrix_a(1_ik)%diagonal(1:5_ik))
    !!   do ii=1_ik,5_ik!size(solver_sparse_matrix_a(1_ik)%diagonal_blocks,dim=3)
    !!   call print("diagnonal 3x3",solver_sparse_matrix_a(1_ik)%diagonal_blocks(:,:,ii))
    !!   end do
    !!  end if
    !!end if
    !fill boundary vector
    call new(boundary_vector,solver_comm%get_from,vtxdist, solver_omps(omp_indx))
    call cel_sp_mat_distr_local_vec(boundary_vector,reduced_boundary_vect,vtxdist,solver_omps(omp_indx))
    !$omp barrier
    if(input_parameters%matrix_scaler .eq. MATRIX_SCALER_DIAGONAL) then
     call cel_sp_mat_distr_boundary_diag_scale(boundary_vector,&
      solver_sparse_matrix_a(1_ik)%scaled_diagonal,&
      solver_omps(omp_indx),solver_omp_shared_work,solver_local_perf,&
      input_parameters%verbosity,err_code)
     !$omp barrier
    endif
    !solve the system
    select case (input_parameters%benchmark_id)
    case (BENCHMARK_CG)
      if(input_parameters%preconditioner .eq. PRECONDITIONER_JACOBI) then
       call cel_cgalg_jacobi(solver_cg,solver_sparse_matrix_a,vtxdist, boundary_vector, &
         vector_x,vector_r,vector_d,vector_q,vector_v,vector_tmp,&
         solver_omps(omp_indx),solver_comm,solver_omp_shared_work, &
         solver_local_perf,perf_counter_mv,perf_counter_vec,perf_profiles,&
         input_parameters%verbosity,err_code)
      else
        call cel_cgalg(solver_cg,solver_sparse_matrix_a,vtxdist, boundary_vector, &
         vector_x,vector_r,vector_d,vector_q,&
         vector_tmp,solver_omps(omp_indx),solver_comm,solver_omp_shared_work, &
         solver_local_perf,perf_counter_mv,perf_counter_vec,perf_profiles,&
         input_parameters%verbosity,err_code)
      endif
      if(solver_omps(omp_indx)%is_master) then
        niter=solver_cg%iter_done
        ftol = solver_cg%norm2
        iierr = err_code
      end if
       call cel_error("Error in cel_elmSolverDoIt  SOLVER_CG !", &
                         &err_code,input_parameters%verbosity,.FALSE._lk) 
    case (BENCHMARK_GMRES)
      if(input_parameters%preconditioner .eq. PRECONDITIONER_JACOBI) then
        call cel_gmresalg_2jacobi(solver_gmres,solver_sparse_matrix_a,vtxdist,boundary_vector,&
          vector_x,vector_r,gmres_vector_v,gmres_vector_z,gmres_vector_w,vector_tmp,&
          gmres_hh,gmres_cc,gmres_ss,gmres_gamma,gmres_alpha,&
          solver_omps(omp_indx),solver_comm,solver_omp_shared_work,&
          solver_local_perf,perf_counter_mv,perf_counter_vec,perf_profiles,&
          input_parameters%verbosity,err_code)
      else              
        call cel_gmresalg_2(solver_gmres,solver_sparse_matrix_a,vtxdist,boundary_vector,&
          vector_x,vector_r,gmres_vector_v,gmres_vector_w,vector_tmp,&
          gmres_hh,gmres_cc,gmres_ss,gmres_gamma,gmres_alpha,&
          solver_omps(omp_indx),solver_comm,solver_omp_shared_work,&
          solver_local_perf,perf_counter_mv,perf_counter_vec,perf_profiles,&
          input_parameters%verbosity,err_code)
      endif
      if(solver_omps(omp_indx)%is_master) then
        niter=solver_gmres%iter_done
        ftol = solver_gmres%residuum_norm
        iierr = err_code
      end if
      call cel_error("Error in cel_elmSolverDoIt  SOLVER_GMRES !", &
                     &err_code,input_parameters%verbosity,.FALSE._lk) 
      case default
        err_code = -4523_ik
        call cel_error("Error in main  solver is undefined!", &
                         &err_code,input_parameters%verbosity,.TRUE._lk)
    end select
    
    !$omp end parallel
   
    sync_out_buff(1)=real(iierr,kind=rk)
    sync_out_buff(2)=real(niter,kind=rk)
    sync_out_buff(3)=real(ftol,kind=rk)
    allocate(send_sync_request(num_reduce_threads-1_ik))
    do ii=1,num_reduce_threads-1_ik
      target_mpik = int(ii,kind=mpik)
      call MPI_Isend(sync_out_buff,3_mpik,mpi_type_rk,target_mpik,target_mpik,&
        elm_reduce_comm,send_sync_request(ii),ierr)
    end do
    do ii=1,num_reduce_threads-1_ik
      !call MPI_Request_free(send_sync_request(ii),ierr)
    end do
    deallocate(send_sync_request)

  else
    call MPI_Irecv(sync_out_buff,3_mpik,mpi_type_rk,0_mpik,mpi_rank_red,elm_reduce_comm,&
      get_sync_request,ierr)
    sync_complited = .false.
    do while(.not. sync_complited)
      call MPI_Test(get_sync_request,sync_complited,MPI_STATUSES_IGNORE,ierr)
      call sleep(1_ik)
     !sync_complited = .true.
    end do
    
    iierr=int(sync_out_buff(1),kind=ik)
    niter=int(sync_out_buff(2),kind=ik)
    ftol=int(sync_out_buff(3),kind=rk)
  endif
 return
end subroutine cel_elmSolverDoIt

subroutine cel_elmSolverReDoIt(bmat,iierr,niter,ftol)

  real(kind=rk),intent(IN) :: bmat(:)
  integer(kind=ik),intent(OUT) :: iierr,niter
  real(kind=rk),intent(OUT) :: ftol
  integer(kind=mpik) , dimension(0_ik) :: dummy_array_int
  real(kind=rk), dimension(0_ik) :: dummy_array_real
  real(kind=rk), allocatable, dimension(:)     :: reduced_boundary_vect
  integer(kind=mpik) :: mpi_rank_red, ierr, err_code, target_mpik
  integer(kind=mpik) :: elm_num_local_rows_mpik, solver_num_local_rows_mpik
  integer(kind=ik) :: ii, omp_indx
  real(kind=rk),dimension(3_ik) :: sync_out_buff
  integer(kind=mpik),dimension(:), allocatable :: send_sync_request
  integer(kind=mpik) :: get_sync_request
  logical(kind=lk) :: sync_complited


  
  call MPI_COMM_RANK(elm_reduce_comm, mpi_rank_red, ierr)
  !collect boundary vector
  elm_num_local_rows_mpik = int(elm_num_local_rows,kind=mpik)
  solver_num_local_rows_mpik = int(solver_num_local_rows,kind=mpik)
  if(mpi_rank_red .eq. 0_ik) then
    allocate(reduced_boundary_vect(solver_num_local_rows))
    call MPI_Gatherv(bmat(1_ik:elm_num_local_rows), elm_num_local_rows_mpik, &
      mpi_type_rk, reduced_boundary_vect(1_ik:solver_num_local_rows),&
      receive_count_rows(1_ik:num_reduce_threads), receive_displ_rows(1_ik:num_reduce_threads),&
      mpi_type_rk,0_mpik,elm_reduce_comm,ierr)
  else
    allocate(reduced_boundary_vect(0_ik))
    call MPI_Gatherv(bmat(1_ik:elm_num_local_rows), elm_num_local_rows_mpik, &
      mpi_type_rk, dummy_array_real,&
      dummy_array_int, dummy_array_int,&
      mpi_type_rk,0_mpik,elm_reduce_comm,ierr)
  end if
  if(mpi_rank_red .eq. 0_ik) then
    !start omp threads
    call omp_set_dynamic(.false._lk)
    call omp_set_nested(.false._lk)
    call omp_set_num_threads(num_reduce_threads)
    !$omp parallel default(none) private(ii,err_code,omp_indx)&
    !$omp& shared(solver_sparse_matrix_a, vtxdist,boundary_vector,vector_x,&
      !$omp& vector_r,vector_d,vector_q,vector_v,vector_tmp,&
      !$omp& gmres_hh,gmres_cc,gmres_ss,gmres_gamma,gmres_alpha,perf_counter_mv,&
      !$omp& solver_local_perf,perf_counter_vec,input_parameters,solver_omps,&
      !$omp& solver_comm,gmres_vector_v,gmres_vector_z,gmres_vector_w,&
      !$omp& num_reduce_threads, solver_cg, solver_gmres,&
      !$omp& solver_omp_shared_work,perf_profiles,reduced_boundary_vect,niter,ftol,iierr )
       
       omp_indx = omp_get_thread_num() + 1
       
    !fill boundary vector
    call cel_sp_mat_distr_local_vec(boundary_vector,reduced_boundary_vect,vtxdist,solver_omps(omp_indx))
    !solve the system
    select case (input_parameters%benchmark_id)
    case (BENCHMARK_CG)
      if(input_parameters%preconditioner .eq. PRECONDITIONER_JACOBI) then
       call cel_cgalg_jacobi(solver_cg,solver_sparse_matrix_a,vtxdist, boundary_vector, &
         vector_x,vector_r,vector_d,vector_q,vector_v,vector_tmp,&
         solver_omps(omp_indx),solver_comm,solver_omp_shared_work, &
         solver_local_perf,perf_counter_mv,perf_counter_vec,perf_profiles,&
         input_parameters%verbosity,err_code)
      else
        call cel_cgalg(solver_cg,solver_sparse_matrix_a,vtxdist, boundary_vector, &
         vector_x,vector_r,vector_d,vector_q,&
         vector_tmp,solver_omps(omp_indx),solver_comm,solver_omp_shared_work, &
         solver_local_perf,perf_counter_mv,perf_counter_vec,perf_profiles,&
         input_parameters%verbosity,err_code)
      endif
      if(solver_omps(omp_indx)%is_master) then
        niter=solver_cg%iter_done
        ftol = solver_cg%norm2
        iierr = err_code
      end if
       call cel_error("Error in cel_elmSolverDoIt  SOLVER_CG !", &
                         &err_code,input_parameters%verbosity,.FALSE._lk) 
    case (BENCHMARK_GMRES)
      if(input_parameters%preconditioner .eq. PRECONDITIONER_JACOBI) then
        call cel_gmresalg_2jacobi(solver_gmres,solver_sparse_matrix_a,vtxdist,boundary_vector,&
          vector_x,vector_r,gmres_vector_v,gmres_vector_z,gmres_vector_w,vector_tmp,&
          gmres_hh,gmres_cc,gmres_ss,gmres_gamma,gmres_alpha,&
          solver_omps(omp_indx),solver_comm,solver_omp_shared_work,&
          solver_local_perf,perf_counter_mv,perf_counter_vec,perf_profiles,&
          input_parameters%verbosity,err_code)
      else              
        call cel_gmresalg_2(solver_gmres,solver_sparse_matrix_a,vtxdist,boundary_vector,&
          vector_x,vector_r,gmres_vector_v,gmres_vector_w,vector_tmp,&
          gmres_hh,gmres_cc,gmres_ss,gmres_gamma,gmres_alpha,&
          solver_omps(omp_indx),solver_comm,solver_omp_shared_work,&
          solver_local_perf,perf_counter_mv,perf_counter_vec,perf_profiles,&
          input_parameters%verbosity,err_code)
      endif
      if(solver_omps(omp_indx)%is_master) then
        niter=solver_gmres%iter_done
        ftol = solver_gmres%residuum_norm
        iierr = err_code
      end if
      call cel_error("Error in cel_elmSolverDoIt  SOLVER_GMRES !", &
                     &err_code,input_parameters%verbosity,.FALSE._lk) 
      case default
        err_code = -4523_ik
        call cel_error("Error in main  solver is undefined!", &
                         &err_code,input_parameters%verbosity,.TRUE._lk)
    end select
    
    !$omp end parallel
    sync_out_buff(1)=real(iierr,kind=rk)
    sync_out_buff(2)=real(niter,kind=rk)
    sync_out_buff(3)=real(ftol,kind=rk)
    allocate(send_sync_request(num_reduce_threads-1_ik))
    do ii=1,num_reduce_threads-1_ik
      target_mpik = int(ii,kind=mpik)
      call MPI_Isend(sync_out_buff,3_mpik,mpi_type_rk,target_mpik,target_mpik,&
        elm_reduce_comm,send_sync_request(ii),ierr)
    end do
    do ii=1,num_reduce_threads-1_ik
      !call MPI_Request_free(send_sync_request(ii),ierr)
    end do
    deallocate(send_sync_request)

  else
    call MPI_Irecv(sync_out_buff,3_mpik,mpi_type_rk,0_mpik,mpi_rank_red,elm_reduce_comm,&
      get_sync_request,ierr)
    sync_complited = .false.
    do while(.not. sync_complited)
      call MPI_Test(get_sync_request,sync_complited,MPI_STATUSES_IGNORE,ierr)
      call sleep(1_ik)
     !sync_complited = .true.
    end do
    
    iierr=int(sync_out_buff(1),kind=ik)
    niter=int(sync_out_buff(2),kind=ik)
    ftol=int(sync_out_buff(3),kind=rk)
  endif
  return
end subroutine

subroutine cel_elmSolverCollectSolution(bmat)
  real(kind=rk),allocatable, dimension(:),intent(OUT) :: bmat(:)
  integer(kind=mpik) :: ierr,mpi_rank_red,mpi_size_red
  integer(kind=mpik) :: elm_num_local_rows_mpik
  real(kind=rk), dimension(0) :: dummy_real
  integer(kind=mpik), dimension(0) :: dummy_int
  integer(kind=ik) :: ii

  call MPI_COMM_RANK(elm_reduce_comm, mpi_rank_red, ierr)
  call MPI_COMM_SIZE(elm_reduce_comm, mpi_size_red, ierr)
  elm_num_local_rows_mpik = int(elm_num_local_rows,kind=mpik)
  allocate(bmat(elm_num_local_rows))

  if(mpi_rank_red .eq. 0_mpik) then

    call MPI_Scatterv(vector_x%values(1_ik:solver_num_local_rows,1_ik), receive_count_rows(1_mpik:mpi_size_red), &
      receive_displ_rows(1_ik:mpi_size_red), mpi_type_rk,&
       bmat(1_ik:elm_num_local_rows),elm_num_local_rows_mpik,mpi_type_rk, 0_mpik,elm_reduce_comm,ierr)
  else
    call MPI_Scatterv(dummy_real, dummy_int, &
      dummy_int, mpi_type_rk,&
      bmat(1_ik:elm_num_local_rows),elm_num_local_rows_mpik,mpi_type_rk, 0_mpik,elm_reduce_comm,ierr)
  end if

return
end subroutine cel_elmSolverCollectSolution

subroutine cel_elmSolverDestroy
  integer(kind=ik) :: ii, length
  integer(kind=mpik) :: err_code,mpi_rank_red

  call MPI_COMM_RANK(elm_reduce_comm, mpi_rank_red, err_code)
  if(mpi_rank_red .eq. 0_mpik) then
    call MPI_Comm_free(solver_mpi_comm,err_code)
    call del(solver_omp_shared_work)
  end if
  call MPI_Comm_free(elm_reduce_comm,err_code)  

  call del(solver_sparse_matrix_a)
  call del(solver_comm )
  if(allocated(solver_omps)) deallocate(solver_omps)
  call del(vector_x)
  if(allocated(vtxdist)) deallocate(vtxdist)
  if(allocated(perf_profiles)) deallocate(perf_profiles)
  if(allocated(receive_count_rows)) deallocate(receive_count_rows)
  if(allocated(receive_displ_rows)) deallocate(receive_displ_rows)
  call del(vector_r)
  call del(vector_d)
  call del(vector_q)
  call del(vector_v)
  if(allocated(gmres_vector_v)) then
    length = size(gmres_vector_v)
    do ii=1_ik, length
      if(allocated(gmres_vector_v(ii)%values)) deallocate(gmres_vector_v(ii)%values)
      if(allocated(gmres_vector_v(ii)%send_comm_buffer)) deallocate(gmres_vector_v(ii)%send_comm_buffer)
      if(allocated(gmres_vector_v(ii)%get_comm_buffer)) deallocate(gmres_vector_v(ii)%get_comm_buffer)
    end do
    deallocate(gmres_vector_v)
  end if
  if(allocated(gmres_vector_z)) then
    length = size(gmres_vector_z)
    do ii=1_ik, length
      if(allocated(gmres_vector_z(ii)%values)) deallocate(gmres_vector_z(ii)%values)
      if(allocated(gmres_vector_z(ii)%send_comm_buffer)) deallocate(gmres_vector_z(ii)%send_comm_buffer)
      if(allocated(gmres_vector_z(ii)%get_comm_buffer)) deallocate(gmres_vector_z(ii)%get_comm_buffer)
    end do
    deallocate(gmres_vector_z)
  end if
  if(allocated(gmres_vector_z)) then
    length = size(gmres_vector_z)
    do ii=1_ik, length
      if(allocated(gmres_vector_z(ii)%values)) deallocate(gmres_vector_z(ii)%values)
      if(allocated(gmres_vector_z(ii)%send_comm_buffer)) deallocate(gmres_vector_z(ii)%send_comm_buffer)
      if(allocated(gmres_vector_z(ii)%get_comm_buffer)) deallocate(gmres_vector_z(ii)%get_comm_buffer)
    end do
    deallocate(gmres_vector_v)
  end if

  call del(vector_tmp)
  call del(boundary_vector)
  if(allocated(gmres_hh)) deallocate(gmres_hh)
  if(allocated(gmres_cc)) deallocate(gmres_cc)
  if(allocated(gmres_ss)) deallocate(gmres_ss)
  if(allocated(gmres_gamma)) deallocate(gmres_gamma)
  if(allocated(gmres_alpha)) deallocate(gmres_alpha)

 return
end subroutine cel_elmSolverDestroy

subroutine cel_elmSolverPerformance
  integer(kind=mpik) :: mpi_rank_red,ierr
  integer(kind=ik) :: omp_indx
  call MPI_COMM_RANK(elm_reduce_comm, mpi_rank_red, ierr)
  
  if(mpi_rank_red .EQ. 0_mpik) then
   !collect some performance statistics and check result
    call omp_set_dynamic(.false._lk)
    call omp_set_nested(.false._lk)
    call omp_set_num_threads(num_reduce_threads)
    !$omp parallel default(none) private(omp_indx)&
    !$omp& shared(solver_total_perf,solver_local_perf,solver_omps,num_reduce_threads)
       
    omp_indx = omp_get_thread_num() + 1
    call reset(solver_total_perf)
    call cel_perf_distr_collect(solver_total_perf,solver_local_perf,solver_omps(omp_indx),.true._lk)
    if(solver_omps(omp_indx)%is_master) then
      call print(solver_total_perf)
    end if
     !$omp end parallel
  end if
end subroutine cel_elmSolverPerformance

end module cel_elmsolver_module
