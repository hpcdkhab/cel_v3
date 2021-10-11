!     
! File:   cel_input_parameters_module.f90
! Project CRESTA (see details on https://www.cresta-project.eu) Exascale library
! Send you email about the bugs to
! The High Performance Computing Center Stuttgart (HLRS) of the University of Stuttgart
! Dmitry Khabi,  (email: khabi@hlrs.de)
! Created on Jun 14, 2013
!

module cel_input_parameters_module
use cel_types_module
use cel_cpu_param_module
use cel_error_module
use cel_extract_command_module
use cel_omp_module
implicit none


type cel_input_parameters_type
  logical(kind=lk) :: verbosity
  integer(kind=ik) :: verbosity_level
  type(cel_omp_type) :: cel_omp
  integer(kind=ik) :: read_matrix
  integer(kind=ik) :: generate_matrix
  type(cel_io_path_type) :: matrix_data_dir
  type(cel_io_path_type) :: comm_graph_dir
  integer(kind=ik) ::nx
  integer(kind=ik) ::ny
  integer(kind=ik) ::nz
  integer(kind=ik) ::dx
  integer(kind=ik) ::dy
  integer(kind=ik) ::dz
  integer(kind=ik) ::num_chunks
  integer(kind=ik) ::chunk_size
  real(kind=rk)    :: increase_factor
  integer(kind=ik) :: tests_num
  logical(kind=lk) :: is_matrix_symmetric
  integer(kind=ik) :: matrix_scaler
  logical(kind=lk) :: matrix_check_null_diagonal
  integer(kind=ik) :: block_rows_num
  integer(kind=ik) :: block_cols_num
  integer(kind=ik) :: attempts_num
  integer(kind=ik) :: prefetcher
  integer(kind=ik) :: write_gnuplot
  logical(kind=lk) :: write_matrix
  logical(kind=lk) :: write_performance
  logical(kind=lk) :: write_profile
  logical(kind=lk) :: write_comm_graph
  integer(kind=ik) :: iter_max
  integer(kind=ik) :: mv_algorithm
  integer(kind=ik) :: mv_algorithm_off_diag
  integer(kind=ik) :: mv_communicator
  integer(kind=ik) :: preconditioner
  integer(kind=ik) :: diagonal_block_size
  real(kind=rk)    :: eps
  real(kind=rk)    ::resolution_duration
  character(len=1024), allocatable, dimension(:) :: matrix_data_dir_array
  character(len=1024) :: performance_filename
  character(len=1024) :: profile_filename
  logical(kind=lk) :: check_result
  integer(kind=ik) :: cpu_freq_idx
  integer(kind=ik) :: benchmark_id
  integer(kind=ik) :: memory_level
  integer(kind=ik) :: gmres_max_restart
  integer(kind=ik) :: num_reduce_threads
end type cel_input_parameters_type

type(cel_input_parameters_type)   :: input_parameters !user input such as type of solver, verbosity level and so on
  
interface def
  module procedure cel_input_parameters_def
end interface def

interface check
  module procedure cel_input_parameters_check
end interface check

interface print
  module procedure cel_input_parameters_print
end interface print

contains 

subroutine cel_input_parameters_print(input_parameters)
type(cel_input_parameters_type),intent(in) :: input_parameters
  write(cel_output_unit,*)"cel library parameters:"
  write(cel_output_unit,*)"verbosity:", input_parameters%verbosity
  write(cel_output_unit,*)"read_matrix:", input_parameters%read_matrix
  write(cel_output_unit,*)"generate_matrix:",input_parameters%generate_matrix
  write(cel_output_unit,*)"nx:",input_parameters%nx
  write(cel_output_unit,*)"ny:",input_parameters%ny
  write(cel_output_unit,*)"nz:",input_parameters%nz
  write(cel_output_unit,*)"dx:",input_parameters%dx
  write(cel_output_unit,*)"dy:",input_parameters%dy
  write(cel_output_unit,*)"dz:",input_parameters%dz
  write(cel_output_unit,*)"num_threads:",input_parameters%cel_omp%num_threads
  write(cel_output_unit,*)"is_matrix_symmetric:",input_parameters%is_matrix_symmetric
  write(cel_output_unit,*)"block_rows_num:",input_parameters%block_rows_num
  write(cel_output_unit,*)"block_cols_num:",input_parameters%block_cols_num
  write(cel_output_unit,*)"attempts_num:",input_parameters%attempts_num
  write(cel_output_unit,*)"prefetcher:",input_parameters%prefetcher
  write(cel_output_unit,*)"increase_factor:",input_parameters%increase_factor
  write(cel_output_unit,*)"tests_num:",input_parameters%tests_num
  write(cel_output_unit,*)"write_gnuplot:",input_parameters%write_gnuplot
  write(cel_output_unit,*)"write_performance:",input_parameters%write_performance
  write(cel_output_unit,*)"write_profile:",input_parameters%write_profile
  write(cel_output_unit,*)"num_chunks:",input_parameters%num_chunks
  write(cel_output_unit,*)"chunk_size:",input_parameters%chunk_size
  write(cel_output_unit,*)"iter_max:",input_parameters%iter_max
  write(cel_output_unit,*)"eps:",input_parameters%eps
  write(cel_output_unit,*)"write_matrix:",input_parameters%write_matrix
  write(cel_output_unit,*)"check_result:",input_parameters%check_result
  write(cel_output_unit,*)"mv_algorithm:",input_parameters%mv_algorithm
  write(cel_output_unit,*)"mv_algorithm_off_diag:",input_parameters%mv_algorithm_off_diag
  write(cel_output_unit,*)"mv_communicator:",input_parameters%mv_communicator
  write(cel_output_unit,*)"cpu_freq_idx:",input_parameters%cpu_freq_idx
  write(cel_output_unit,*)"benchmark_id:",input_parameters%benchmark_id
  write(cel_output_unit,*)"memory_level:",input_parameters%memory_level
  write(cel_output_unit,*)"write_comm_graph:",input_parameters%write_comm_graph
  write(cel_output_unit,*)"matrix_check_null_diagonal:",input_parameters%matrix_check_null_diagonal
  write(cel_output_unit,*)"preconditioner:",input_parameters%preconditioner
  write(cel_output_unit,*)"matrix_scaler:",input_parameters%matrix_scaler
  write(cel_output_unit,*)"matrix_scaler:",input_parameters%gmres_max_restart
  write(cel_output_unit,*)"num_reduce_threads:",input_parameters%num_reduce_threads
  write(cel_output_unit,*)"diagonal_block_size:",input_parameters%diagonal_block_size
end subroutine cel_input_parameters_print

subroutine cel_input_parameters_def(input_parameters)
type(cel_input_parameters_type) :: input_parameters

  input_parameters%verbosity = .FALSE._lk
  input_parameters%read_matrix = 0_ik
  input_parameters%generate_matrix = 1_ik
  input_parameters%nx = 0_ik
  input_parameters%ny = 0_ik
  input_parameters%nz = 0_ik
  input_parameters%dx = 0_ik
  input_parameters%dy = 0_ik
  input_parameters%dz = 0_ik
  input_parameters%num_reduce_threads= cel_cpu_param_cores
  input_parameters%cel_omp%num_threads = 2_ik
  input_parameters%is_matrix_symmetric = .FALSE._lk
  input_parameters%block_rows_num = 1_ik
  input_parameters%block_cols_num = 1_ik
  input_parameters%attempts_num = 1_ik
  input_parameters%prefetcher = 0_ik
  input_parameters%increase_factor = 1.0_rk
  input_parameters%tests_num = 1_ik
  input_parameters%write_gnuplot = 0_ik
  input_parameters%write_performance = .FALSE._lk
  input_parameters%write_profile = .FALSE._lk
  input_parameters%num_chunks = 8_ik
  input_parameters%chunk_size = 8_ik
  input_parameters%iter_max = 1000_ik
  input_parameters%eps = 0.1e-06_rk
  input_parameters%write_matrix = .FALSE._lk
  input_parameters%check_result = .FALSE._lk
  input_parameters%mv_algorithm = MV_CSR
  input_parameters%mv_algorithm_off_diag = MV_COO
  input_parameters%mv_communicator = MV_COMM_ISEND_COMM_BUFFER
  input_parameters%cpu_freq_idx = 0
  input_parameters%benchmark_id  = BENCHMARK_MV
  input_parameters%memory_level  = 1
  input_parameters%write_comm_graph = .FALSE._lk
  input_parameters%matrix_check_null_diagonal = .FALSE._lk
  input_parameters%preconditioner = preconditioner_JACOBI
  input_parameters%diagonal_block_size  = 1_ik
  input_parameters%matrix_scaler  = MATRIX_SCALER_NULL
  input_parameters%gmres_max_restart  = GMRES_MAX_RESTART_DEF
end subroutine cel_input_parameters_def


subroutine cel_input_parameters_check(input_parameters)
type(cel_input_parameters_type), intent(in) :: input_parameters
  logical(kind=lk) :: flag
  integer(kind=ik) :: aa, bb

 ! if(input_parameters%read_matrix .EQ. 0 .AND. &
 !    input_parameters%generate_matrix .EQ. 0) then
 !   call cel_error("One of the parameters read_matrix or generate_matrix must be defined",&
 !    1_ik, .TRUE._lk, .TRUE._lk)
 ! end if
  if(input_parameters%cel_omp%num_threads .LE. 0_ik .OR. &
    input_parameters%cel_omp%num_threads .GT. cel_cpu_param_max_num_threads) then
    call cel_error("Parameters omp_threads has wrong value",&
     1_ik, .TRUE._lk, .TRUE._lk)
  end if
  if(input_parameters%attempts_num .LE. 0 ) then
    call cel_error("Parameters attempts_num has wrong value",&
     1_ik, .TRUE._lk, .TRUE._lk)
  end if
  if(input_parameters%generate_matrix .EQ. 1) then
    aa = input_parameters%nx * input_parameters%ny* &
      input_parameters%nz
    bb = input_parameters%nx * input_parameters%ny* &
      input_parameters%nz
    if(aa .LE. 0 ) then
          call cel_error("Parameters nx, ny, nz are wrong",&
       1_ik, .TRUE._lk, .TRUE._lk)
    end if
    if(bb .LE. 0 ) then
          call cel_error("Parameters dx, dy, dz are wrong",&
       1_ik, .TRUE._lk, .TRUE._lk)
    end if
    aa = input_parameters%nx/input_parameters%dx
    bb = aa*input_parameters%dx
    flag = bb .NE. input_parameters%nx
    if(flag) then
      call cel_error("Parameters nx and dx are not compatible",&
       1_ik, .TRUE._lk, .TRUE._lk)
    end if
    aa = input_parameters%ny/input_parameters%dy
    bb = aa*input_parameters%dy
    flag = bb .NE. input_parameters%ny
    if(flag) then
      call cel_error("Parameters ny and dy are not compatible",&
       1_ik, .TRUE._lk, .TRUE._lk)
    end if
    aa = input_parameters%nz/input_parameters%dz
    bb = aa*input_parameters%dz
    flag = bb .NE. input_parameters%nz
    if(flag) then
      call cel_error("Parameters nz and dz are not compatible",&
       1_ik, .TRUE._lk, .TRUE._lk)
    end if
  endif
!  if(input_parameters%read_matrix .EQ. 1 ) then
!     call cel_error("read matrix is not implemented",&
!       1_ik, .TRUE._lk, .TRUE._lk)
!  end if

end subroutine cel_input_parameters_check

subroutine cel_input_parameters_read(cel_input_parameters)
  type(cel_input_parameters_type), intent(out) :: cel_input_parameters
  character(len=2048) :: cmd
  character(len=2048) :: message
  integer(kind=ik) :: help
  integer(kind=ik) :: int_temp
  help = -1_ik
  int_temp = 0_ik
  call def(cel_input_parameters)
  call get_command(cmd)
  !write(*,'(3a)') 'command_line is "',trim(cmd),'"'
  message='parameters of matrix vector multiplication' &
      &//' -threads <int>' &
      &//' -verbosity <int> ' &
      &//' -readmatrix <int>' &
      &//' -generatematrix <int>' &
      &//' -topdir <filenames> ' &
      &//' -matrixdir <filenames> ' &
      &//' -filename <filenames> ' &
      &//' -comm_graph_dir <filenames>' &
      &//' -nx <int>' &
      &//' -ny <int>' &
      &//' -nz <int>' &
      &//' -dx <int>' &
      &//' -dy <int>' &
      &//' -dz <int>' &
      &//' -mv_algorithm  <int>' &
      &//' -mv_algorithm_off_diag <int>' &
      &//' -benchmark_id  <int>' &
      &//' -mv_communicator  <int>' &
      &//' -matrix_scaler <int>' &
      &//' -matrix_check_null_diagonal <int>' &
      &//' -preconditioner  <int>' &
      &//' -cpu_freq_idx  <int>' &
      &//' -num_chunks <int>' &
      &//' -chunk_size <int>' &
      &//' -iter_max <int>' &
      &//' -block_rows_num <int>' &
      &//' -block_cols_num <int>' &
      &//' -eps <real>' &
      &//' -mv <int>' &
      &//' -attempts <int>' &
      &//' -prefetcher <int>' &
      &//' -increase_factor <real>' &
      &//' -tests_num <int>' &
      &//' -write_gnuplot <int>' &
      &//' -write_performance <int>' &
      &//' -write_profile <int>' &
      &//' -check_result <int>' &
      &//' -write_matrix <int>' &
      &//' -write_comm_graph' &
      &//' -performance_filename <filenames> ' &
      &//' -profile_filename <filenames> ' &
      &//' -memory_level <int>' &
      &//' -gmres_max_restart <int>' &
      &//' -num_reduce_threads <int>' &
      &//' -diagonal_block_size <int>' &
      &//' -help <int>'

  if(cel_extract_command_parameter(cmd,'-help',stop_on_error=.false.,&
                value=help, syntax=message) ==0) then
      if(help .GE. 0) then
        write (*,'(A)') trim(message)
        if(help .GT. 0) then
          stop
        endif
      end if
  endif
  
  if(cel_extract_command_parameter(cmd,'-threads',stop_on_error=.false.,&
                value=cel_input_parameters%cel_omp%num_threads, syntax=message) /=0) then
  endif
  if(cel_extract_command_parameter(cmd,'-num_reduce_threads',stop_on_error=.false.,&
                value=cel_input_parameters%num_reduce_threads, syntax=message) /=0) then
  endif
  if(cel_extract_command_parameter(cmd,'-verbosity',stop_on_error=.false.,&
                value=cel_input_parameters%verbosity_level,syntax=message) /=0) then
  endif
  if(cel_extract_command_parameter(cmd,'-readmatrix',stop_on_error=.false.,&
                value=cel_input_parameters%read_matrix,syntax=message) /=0) then
  endif
  if(cel_extract_command_parameter(cmd,'-generatematrix',stop_on_error=.false.,&
                value=cel_input_parameters%generate_matrix,syntax=message) /=0) then
  endif
  if(cel_extract_command_parameter(cmd,'-topdir',stop_on_error=.false.,&
                value=cel_input_parameters%matrix_data_dir%top_directory_name,&
                syntax=message) /=0) then
  endif
  if(cel_extract_command_parameter(cmd,'-comm_graph_dir',stop_on_error=.false.,&
                value=cel_input_parameters%comm_graph_dir%directory_name,&
                syntax=message) /=0) then
  endif
  if(cel_extract_command_parameter(cmd,'-matrixdir',stop_on_error=.false.,&
                value=cel_input_parameters%matrix_data_dir%directory_name,&
                syntax=message) /=0) then
  endif
  if(cel_extract_command_parameter(cmd,'-filename',stop_on_error=.false.,&
                value=cel_input_parameters%matrix_data_dir%file_name,syntax=message) /=0) then
  endif
  if(cel_extract_command_parameter(cmd,'-performance_filename',stop_on_error=.false.,&
                value=cel_input_parameters%performance_filename,syntax=message) /=0) then
  endif
  if(cel_extract_command_parameter(cmd,'-profile_filename',stop_on_error=.false.,&
                value=cel_input_parameters%profile_filename,syntax=message) /=0) then
  endif
  if(cel_extract_command_parameter(cmd,'-nx',stop_on_error=.false.,&
                value=cel_input_parameters%nx,syntax=message) /=0) then
  endif
   if(cel_extract_command_parameter(cmd,'-ny',stop_on_error=.false.,&
                value=cel_input_parameters%ny,syntax=message) /=0) then
  endif
  if(cel_extract_command_parameter(cmd,'-nz',stop_on_error=.false.,&
                value=cel_input_parameters%nz,syntax=message) /=0) then
  endif
  if(cel_extract_command_parameter(cmd,'-dx',stop_on_error=.false.,&
                value=cel_input_parameters%dx,syntax=message) /=0) then
  endif
   if(cel_extract_command_parameter(cmd,'-dy',stop_on_error=.false.,&
                value=cel_input_parameters%dy,syntax=message) /=0) then
  endif
  if(cel_extract_command_parameter(cmd,'-dz',stop_on_error=.false.,&
                value=cel_input_parameters%dz,syntax=message) /=0) then
  endif
  if(cel_extract_command_parameter(cmd,'-mv_algorithm',stop_on_error=.false.,&
                value=cel_input_parameters%mv_algorithm,syntax=message) /=0) then
  endif
  if(cel_extract_command_parameter(cmd,'-mv_algorithm_off_diag',stop_on_error=.false.,&
                value=cel_input_parameters%mv_algorithm_off_diag,syntax=message) /=0) then
  endif  
  if(cel_extract_command_parameter(cmd,'-mv_communicator',stop_on_error=.false.,&
                value=cel_input_parameters%mv_communicator,syntax=message) /=0) then
  endif
  if(cel_extract_command_parameter(cmd,'-preconditioner',stop_on_error=.false.,&
                value=cel_input_parameters%preconditioner,syntax=message) /=0) then
  endif
  if(cel_extract_command_parameter(cmd,'-cpu_freq_idx',stop_on_error=.false.,&
                value=cel_input_parameters%cpu_freq_idx,syntax=message) /=0) then
  endif
  if(cel_extract_command_parameter(cmd,'-num_chunks',stop_on_error=.false.,&
                value=cel_input_parameters%num_chunks,syntax=message) /=0) then
  endif
  if(cel_extract_command_parameter(cmd,'-chunk_size',stop_on_error=.false.,&
                value=cel_input_parameters%chunk_size,syntax=message) /=0) then
  endif
  if(cel_extract_command_parameter(cmd,'-iter_max',stop_on_error=.false.,&
                value=cel_input_parameters%iter_max,syntax=message) /=0) then
  endif
  if(cel_extract_command_parameter(cmd,'-block_rows_num',stop_on_error=.false.,&
                value=cel_input_parameters%block_rows_num,syntax=message) /=0) then
  endif
  if(cel_extract_command_parameter(cmd,'-block_cols_num',stop_on_error=.false.,&
                value=cel_input_parameters%block_cols_num,syntax=message) /=0) then
  endif
  if(cel_extract_command_parameter(cmd,'-eps',stop_on_error=.false.,&
                value=cel_input_parameters%eps,syntax=message) /=0) then
  endif
  if(cel_extract_command_parameter(cmd,'-attempts',stop_on_error=.false.,&
                value=cel_input_parameters%attempts_num,syntax=message) /=0) then
  endif
  if(cel_extract_command_parameter(cmd,'-threads',stop_on_error=.false.,&
                value=cel_input_parameters%cel_omp%num_threads,syntax=message) /=0) then
  endif
  if(cel_extract_command_parameter(cmd,'-prefetcher',stop_on_error=.false.,&
                value=cel_input_parameters%prefetcher,syntax=message) /=0) then
  endif
  if(cel_extract_command_parameter(cmd,'-increase_factor',stop_on_error=.false.,&
                value=cel_input_parameters%increase_factor,syntax=message) /=0) then
  endif
  if(cel_extract_command_parameter(cmd,'-matrix_scaler',stop_on_error=.false.,&
                value=cel_input_parameters%matrix_scaler,syntax=message) /=0) then
  endif
  if(cel_extract_command_parameter(cmd,'-tests_num',stop_on_error=.false.,&
                value=cel_input_parameters%tests_num,syntax=message) /=0) then
  endif
  if(cel_extract_command_parameter(cmd,'-write_gnuplot',stop_on_error=.false.,&
                value=cel_input_parameters%write_gnuplot,syntax=message) /=0) then
  endif
  int_temp = -1_ik
  if(cel_extract_command_parameter(cmd,'-write_performance',stop_on_error=.false.,&
                value=int_temp,syntax=message) /=0) then
  endif
  if(int_temp .GT. 0_ik) cel_input_parameters%write_performance = .TRUE._lk
  if(int_temp .EQ. 0_ik) cel_input_parameters%write_performance = .false._lk
  int_temp = -1_ik
  if(cel_extract_command_parameter(cmd,'-write_profile',stop_on_error=.false.,&
                value=int_temp,syntax=message) /=0) then
  endif
  if(int_temp .GT. 0_ik) cel_input_parameters%write_profile = .TRUE._lk
  if(int_temp .EQ. 0_ik) cel_input_parameters%write_profile = .false._lk
  int_temp = -1_ik
  if(cel_extract_command_parameter(cmd,'-matrix_check_null_diagonal',stop_on_error=.false.,&
                value=int_temp,syntax=message) /=0) then
  endif
  if(int_temp .GT. 0_ik) cel_input_parameters%matrix_check_null_diagonal = .TRUE._lk
  if(int_temp .EQ. 0_ik) cel_input_parameters%matrix_check_null_diagonal = .false._lk
  int_temp = -1_ik
  if(cel_extract_command_parameter(cmd,'-check_result',stop_on_error=.false.,&
                value=int_temp,syntax=message) /=0) then
  endif
  if(int_temp .GT. 0_ik) cel_input_parameters%check_result = .TRUE._lk
  if(int_temp .EQ. 0_ik) cel_input_parameters%check_result = .false._lk
  int_temp = -1_ik
  if(cel_extract_command_parameter(cmd,'-write_matrix',stop_on_error=.false.,&
                value=int_temp,syntax=message) /=0) then
  endif
  if(int_temp .GT. 0_ik) cel_input_parameters%write_matrix = .TRUE._lk
  int_temp = -1_ik
  if(cel_extract_command_parameter(cmd,'-write_comm_graph',stop_on_error=.false.,&
                value=int_temp,syntax=message) /=0) then
  endif
  if(int_temp .GT. 0_ik) cel_input_parameters%write_comm_graph = .TRUE._lk
  if(int_temp .EQ. 0_ik) cel_input_parameters%write_comm_graph = .false._lk
  if(cel_extract_command_parameter(cmd,'-memory_level',stop_on_error=.false.,&
                value=cel_input_parameters%memory_level,syntax=message) /=0) then
  endif
  if(cel_extract_command_parameter(cmd,'-benchmark_id',stop_on_error=.false.,&
                value=cel_input_parameters%benchmark_id,syntax=message) /=0) then
  endif
  if(cel_extract_command_parameter(cmd,'-gmres_max_restart',stop_on_error=.false.,&
                value=cel_input_parameters%gmres_max_restart,syntax=message) /=0) then
  endif
  if(cel_extract_command_parameter(cmd,'-diagonal_block_size',stop_on_error=.false.,&
                value=cel_input_parameters%diagonal_block_size,syntax=message) /=0) then
  endif
  if(cel_input_parameters%verbosity_level .GT. 0_ik) cel_input_parameters%verbosity = .TRUE._lk
  if(cel_input_parameters%read_matrix .GT. 2) then
    call cel_input_parameters_matrix_data_dir(cel_input_parameters%matrix_data_dir_array, &
      &'matrix_data_dir.data', cel_input_parameters)
    cel_input_parameters%tests_num = size(cel_input_parameters%matrix_data_dir_array)
  end if

  call cel_input_parameters_check(cel_input_parameters)
end subroutine cel_input_parameters_read

subroutine cel_input_parameters_matrix_data_dir(matrix_data_dir, parameter_file, input_parameters)
  character(len=1024), allocatable, dimension(:), intent(out) :: matrix_data_dir
  character(len=*), intent(in) :: parameter_file
  type(cel_input_parameters_type), intent(in) :: input_parameters
  integer :: parameter_unit, num_matrix_dirs
  logical :: reading_successfull
  integer(kind=ik) :: alloc_stat, ii
  integer,parameter                        :: mx_cases=20
  integer                                  :: max_cases, case_type, where_cross, iostat
  character(len=1024),dimension(mx_cases)    :: cases
  character(len=1024) :: string

  reading_successfull=.true.
  open(newunit=parameter_unit,iostat=iostat,file=parameter_file,action='read')
  if(iostat /=  0 ) then
    write (*,'(2(a))') 'error opening the file with matrix directories', parameter_file
    print*,'may be that file is not present'
    print*,'must contain paths to the matrix data'
    reading_successfull=.false.
  else
    cases=(/('',ii=1,mx_cases)/)
    max_cases=mx_cases
    case_type=0
    do ii=1,mx_cases
      read(parameter_unit,*,iostat=iostat) string
      if(iostat /= 0 ) exit
      where_cross=index(string,'#')
      if(where_cross>0) string=string(1:where_cross-1)
      if(trim(string) == '') cycle
      case_type=case_type+1
      cases(case_type)=trim(string)
      max_cases=case_type
    enddo
    num_matrix_dirs = max_cases
    if(num_matrix_dirs .GE. 0) then
      allocate(matrix_data_dir(num_matrix_dirs), stat=alloc_stat)
    else
      call cel_error("cel_input_parameters_matrix_data_dir: check options file or delete it", &
        1_ik, .TRUE._lk, .TRUE._lk)
    end if
    matrix_data_dir = cases(1:max_cases)
  end if
    
  if(.not. reading_successfull) then
    if(iostat /=  0 ) then
    print'(a,a)','error reading the file ',trim(parameter_file)
    print'(a)','must contain '
    print'(a)','line 1 : path to 1-th matrix directory'
    print'(a)','line n : path to n-th matrix directory'
    reading_successfull=.false.
  endif 

  if(.not.reading_successfull) then
    close(parameter_unit)
    open(newunit=parameter_unit,file=parameter_file)
    print*
    print'(a,a)','a default parameter will be written to ',trim(parameter_file)
    write(parameter_unit,'(a)') &
      '"/home/hpcdkhab/matrix_bcrs_original/io_bone_matrix_second_order_block3x3_el_1641x1641"'
    write(parameter_unit,'(a)') &
      '"/home/hpcdkhab/matrix_bcrs_original/io_bone_matrix_second_order_block3x3_el_51279x51279"'
    write(parameter_unit,'(a)') &
      '"/home/hpcdkhab/matrix_bcrs_original/io_bone_matrix_second_order_block3x3_el_163611x163611"'
    write(parameter_unit,'(a)') &
      '"/home/hpcdkhab/matrix_bcrs_original/io_bone_matrix_second_order_block3x3_el_439737x439737"'
    write(parameter_unit,'(a)') &
      '"/home/hpcdkhab/matrix_bcrs_original/io_bone_matrix_second_order_block3x3_el_543816x543816"'
    write(parameter_unit,'(a)') &
      '"/home/hpcdkhab/matrix_bcrs_original/io_bone_matrix_second_order_block3x3_el_601185x601185"'
    stop
  endif 
end if


end subroutine cel_input_parameters_matrix_data_dir



end module cel_input_parameters_module
