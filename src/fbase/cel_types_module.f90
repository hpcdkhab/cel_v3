!
! File:   cel_types_module.f90
! Project CRESTA (see details on https://cresta-project.eu) Exascale library
! Send you email about the bugs to
! The High Performance Computing Center Stuttgart (HLRS) of the University of Stuttgart
! Dmitry Khabi,  (email: khabi@hlrs.de)
! Created on May 29, 2013
!

module cel_types_module
use, intrinsic :: ISO_C_BINDING
implicit none
  integer, parameter :: rk = 8 ! type for real
  integer, parameter :: ik = 4 ! type for index format
  integer, parameter :: c_rk = C_DOUBLE !shoul be changed by changing of rk
  integer, parameter :: c_ik = C_INT!C_INT!C_LONG !shoul be changed by changing of ik
  integer, parameter :: mcl = 1024!
  integer, parameter :: itk = 8! timer

  integer, parameter :: c_enum = 4 !enumeration type
  integer, parameter :: ceik = 4 ! type of c error int
  integer, parameter :: crk = 8 ! type for nodes coordinates pckage base (cel c interface)

  integer, parameter :: mpik = 4 !mpi integer format
    
  integer, parameter :: lk = 4 !type for logical
  
  logical(kind=lk) :: cel_is_in_debug = .TRUE._lk !global debug switcher
  
  integer :: cel_output_unit=6
  integer :: cel_error_unit=0

  real(kind=rk), parameter :: cel_real_small = 1.1_rk*EPSILON(cel_real_small)
  real(kind=rk), parameter :: cel_real_undef = huge(1.0_rk) !undefined real value
  integer(kind=ik), parameter :: cel_index_undef = huge(1_ik)  !undefined index value
  integer(kind=ik), parameter :: cel_index_max = huge(1_ik)-1_ik !max index value
  character(len=1),parameter    :: cel_delimiter='/'    ! delimeter for Unix, Linux
  !character(len=1),parameter    :: cel_delimiter='\'    ! delimeter for Windows
  real(kind=rk) :: cel_tsc_cycles  = 1.0_rk !processor's tick frequency  cel_timer_module.clock_tics_per_second()
  
  integer, parameter :: real_kind = rk
  doubleprecision ref_variable
  doubleprecision db_variable
  integer(kind=ik),parameter                         :: r_high=kind(db_variable)
  integer(kind=ik),parameter                         :: complex_kind=kind(ref_variable)
  integer(kind=ik),parameter                         :: str_length=1024
  real(kind=rk),parameter                            :: Pi=3.14159265358979323846264338328_rk
  complex(kind=complex_kind),parameter               :: imag=(0._real_kind,1._rk)
  integer(kind=ik),parameter                         :: integr=1
  integer(kind=ik),parameter                         :: infinity=huge(integr)
  character(len=53),parameter ::  &
&  alpha_chars='ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz_'
!!              123456789 123456789 123456789 123456789 123456789 123456789
  character(len=10),parameter  :: digits='0123456789'
  character(len=1),parameter  :: nonsignificant_chars=' '
  character(len=len(alpha_chars)+len(digits)),parameter  :: variable_chars=alpha_chars//digits

  integer(kind=ik), parameter :: VRL = 64_ik !vector length for jad matrix vector multiplication
  integer(kind=ik), parameter :: MV_CSR=1_ik!matrix vector multiplication with csr algorithm
  integer(kind=ik), parameter :: MV_JAD=2_ik!matrix vector multiplication with jad algorithm
  integer(kind=ik), parameter :: MV_COO=3_ik!matrix vector multiplication with coo algorithm
  integer(kind=ik), parameter :: MV_COMM_ISEND=1_ik!matrix vector with isend
  integer(kind=ik), parameter :: MV_COMM_ISEND_GROUPS=2_ik!matrix vector with isend within smal groups
  integer(kind=ik), parameter :: MV_COMM_IBCAST=3_ik!matrix vector with ibcast within smal groups
  integer(kind=ik), parameter :: MV_COMM_ISEND_IDX=4_ik!matrix vector with isend and cutted buffer (first and last idx of column in matrix sub-block)
  integer(kind=ik), parameter :: MV_COMM_ISEND_COMM_BUFFER=5_ik!matrix vector with isend and additional buffer, only needed elements will be send / get from the neighbours procs
  integer(kind=ik), parameter :: MV_COMM_ISEND_IDX_AND_COMM_BUFFER=6_ik!matrix vector combination of MV_COMM_ISEND_COMM_BUFFER and MV_COMM_ISEND_IDX !not implemented
  integer(kind=ik), parameter :: BENCHMARK_MV=1_ik!benchmark matrix vector multiplication
  integer(kind=ik), parameter :: BENCHMARK_R_D=2_ik!benchmark matrix vector multiplication
  integer(kind=ik), parameter :: BENCHMARK_NORM2=3_ik!benchmark norm 2 of vector (dot product)
  integer(kind=ik), parameter :: BENCHMARK_CG=4_ik!benchmark cg
  integer(kind=ik), parameter :: BENCHMARK_GMRES=5_ik!benchmark gmres
  integer(kind=ik), parameter :: PRECONDITIONER_NULL=0_ik !no preconditioner for preconditioned cg
  integer(kind=ik), parameter :: PRECONDITIONER_JACOBI=1_ik !jacobi preconditioner for preconditioned cg or gmres
  integer(kind=ik), parameter :: MATRIX_SCALER_NULL=0_ik !scaling will be not done
  integer(kind=ik), parameter :: MATRIX_SCALER_DIAGONAL=1_ik !scaling of the matrix with diagonal (1/sqrt(abs(d_i)))

  integer(kind=ik), parameter :: GMRES_MAX_RESTART_DEF=10_ik !diefined max restart for gmres
  integer(kind=ik), parameter :: GMRES_REASON_UNDEX =0_ik !gmres - undef
  integer(kind=ik), parameter :: GMRES_REASON_X0 =1_ik !gmres - the initial x vector is a solution
  integer(kind=ik), parameter :: GMRES_REASON_XN =2_ik !gmres - the initial x vector is a found solution


  integer(kind=ik), parameter :: AW_SMALLEST_OMP_TASK = 4
  integer(kind=ik), parameter :: LOOP_SLEEP_DEFLENGTH=2_ik!about  Xx10 + 4 locks theoretical cycles (nop operations) + 1 write 8 bytes
  integer(kind=ik), parameter :: LOOP_SLEEP_DEFSHORTLENGTH=1_ik!about  Xx10  theoretical cycles (nop operations) + 1 write 8 bytes
  
  integer(kind=ik), parameter :: STRING_MAX_LENGTH=1024_ik!max string length
  !Path to the files
  type cel_io_path_type
    character(len=512)                             :: file_name !file name
    character(len=512)                             :: directory_name !directory name
    character(len=512)                             :: top_directory_name !directory name
  end type cel_io_path_type

  
end module cel_types_module
