!     
! File:   cel_cg_module.f90
! Project CRESTA (see details on https://cresta-project.eu) Exascale library
! Send you email about the bugs to
! The High Performance Computing Center Stuttgart (HLRS) of the University of Stuttgart
! Dmitry Khabi,  (email: khabi@hlrs.de)
! Created on May 7, 2014
!

module cel_alg_module
use cel_types_module
use cel_timer_interface_module
use cel_cpu_param_module
use cel_sp_mat_distr_vec_module
use cel_sp_mat_distr_module
use cel_perf_module
use cel_profile_module
use cel_base_print_module
use cel_omp_module
use cel_omp_shared_work_module
use cel_comm_module
use MPI
use OMP_LIB
implicit none

type cel_cgalg_parameter_type
  real(kind=rk)                                     :: sigma_new
  real(kind=rk)                                     :: sigma_old
  real(kind=rk)                                     :: betta
  real(kind=rk)                                     :: alpha
  real(kind=rk)                                     :: scalar_d_q
  real(kind=rk)                                     :: real_tmp
  real(kind=rk)                                     :: norm2
  real(kind=rk)                                     :: eps
  real(kind=rk)                                     :: residuum
  integer(kind=ik)                                  :: iter_max
  integer(kind=ik)                                  :: iter_done
  logical(kind=lk)                                  :: verbosity
  integer(kind=ik)                                  :: proc_out
end type cel_cgalg_parameter_type
 
type cel_gmresalg_parameter_type
  real(kind=rk)                                     :: residuum_norm
  real(kind=rk)                                     :: betta
  real(kind=rk)                                     :: eps
  real(kind=rk)                                     :: residuum
  integer(kind=ik)                                  :: restart_curr
  integer(kind=ik)                                  :: restart_max
  logical(kind=lk)                                  :: verbosity
  integer(kind=ik)                                  :: proc_out
  integer(kind=ik)                                  :: iter_max
  integer(kind=ik)                                  :: iter_done
  integer(kind=ik)                                  :: conv_reason
end type cel_gmresalg_parameter_type

interface new
  module procedure cel_cgalg_parameter_new
  module procedure cel_gmresalg_parameter_new
end interface new

interface del
  module procedure cel_cgalg_parameter_del
  module procedure cel_gmresalg_parameter_del
end interface del

contains
!initialization of main cg parameters, vectors etc.
subroutine cel_cgalg_parameter_new(cel_cgalg_par,vtxdist,&
  vector_x,vector_r,vector_d,vector_q,vector_v,vector_tmp,&
  cg_perf_counter,cg_perf_counter_mv, cg_perf_counter_vec,cg_cel_profiles,&
  cel_omp,cel_comm,input_parameters,err_code)
  type(cel_cgalg_parameter_type), intent(out)                    :: cel_cgalg_par
  integer(kind=ik),dimension(:), intent(in)                      :: vtxdist
  type(cel_sp_mat_distr_vec_type)  , intent(inout)               :: vector_x
  type(cel_sp_mat_distr_vec_type) , intent(inout)                :: vector_r
  type(cel_sp_mat_distr_vec_type)  , intent(inout)               :: vector_d
  type(cel_sp_mat_distr_vec_type) , intent(inout)                :: vector_q
  type(cel_sp_mat_distr_vec_type) , intent(inout)                :: vector_v
  type(cel_sp_mat_distr_vec_type) , intent(inout)                :: vector_tmp
  type(cel_perf_counter_type)  , intent(inout)                   :: cg_perf_counter
  type(cel_perf_counter_type)  , intent(inout)                   :: cg_perf_counter_mv
  type(cel_perf_counter_type)   , intent(inout)                  :: cg_perf_counter_vec
  type(cel_profile_type),dimension(:),allocatable,intent(inout)  :: cg_cel_profiles
  type(cel_omp_type), intent(in)                                 :: cel_omp
  type(cel_comm_type), intent(in)                                :: cel_comm
  type(cel_input_parameters_type), intent(in)                    :: input_parameters
  integer(kind=ik), intent(out)                                  :: err_code
  integer(kind=ik) :: ii, length_procs, length_vtxdist,sys_stat
  integer(kind=ik) :: my_length, max_length, tmp_int, num_elements
  
  err_code = 0_ik

  call cel_sp_mat_distr_vec_new_from_procs(vector_x,cel_comm%get_from,vtxdist,cel_omp)
  call cel_sp_mat_distr_vec_new_from_procs(vector_r,cel_comm%get_from,vtxdist,cel_omp)
  call cel_sp_mat_distr_vec_new_from_procs(vector_d,cel_comm%get_from,vtxdist,cel_omp)
  call cel_sp_mat_distr_vec_new_from_procs(vector_q,cel_comm%get_from,vtxdist,cel_omp)
  call cel_sp_mat_distr_vec_new_from_procs(vector_v,cel_comm%get_from,vtxdist,cel_omp)
  call cel_sp_mat_distr_vec_new_from_procs(vector_tmp,cel_comm%get_from,vtxdist,cel_omp)

  if(cel_omp%is_master) then
    cel_cgalg_par%sigma_new = 0.0_rk
    cel_cgalg_par%sigma_old = 0.0_rk
    cel_cgalg_par%betta = 0.0_rk
    cel_cgalg_par%alpha = 0.0_rk
    cel_cgalg_par%scalar_d_q = 0.0_rk
    cel_cgalg_par%real_tmp = 0.0_rk
    cel_cgalg_par%norm2 = 0.0_rk
    cel_cgalg_par%eps = input_parameters%eps
    cel_cgalg_par%residuum = 0.0_rk
    cel_cgalg_par%iter_max = input_parameters%iter_max
    cel_cgalg_par%iter_done = 0_ik
    cel_cgalg_par%verbosity = .false._lk
    cel_cgalg_par%proc_out = 0_ik
    call cel_perf_counter_reset(cg_perf_counter)
    call cel_perf_counter_reset(cg_perf_counter_vec)
    call cel_perf_counter_reset(cg_perf_counter_mv)
    my_length = vtxdist(cel_omp%master_num+2_ik)-vtxdist(cel_omp%master_num+1_ik)
    cg_perf_counter%size = my_length
   ! cg_perf_counter%num_elements = num_elements
   if(input_parameters%write_profile) then
    allocate(cg_cel_profiles(input_parameters%iter_max), stat=sys_stat)
        call cel_error("Error cel_cg_parameter_new alloc ", sys_stat,&
                        .TRUE._lk,.TRUE._lk)
   end if
  end if
  !$omp barrier
end subroutine cel_cgalg_parameter_new


!delete of main cg parameters, vectors etc.
subroutine cel_cgalg_parameter_del(cel_cgalg_par)
  type(cel_cgalg_parameter_type)                 :: cel_cgalg_par
  

end subroutine cel_cgalg_parameter_del


!initialization of main gmres parameters, vectors etc.
subroutine cel_gmresalg_parameter_new(cel_gmres,vtxdist,&
  vector_x,vector_r,vector_v,vector_z,vector_w,vector_tmp,&
  cel_gmres_hh,cel_gmres_cc,cel_gmres_ss,cel_gmres_gamma,cel_gmres_alpha,&
  gmres_perf_counter,gmres_perf_counter_mv, gmres_perf_counter_vec,gmres_cel_profiles,&
  cel_omp,cel_comm,input_parameters,err_code)
  type(cel_gmresalg_parameter_type), intent(out)                         :: cel_gmres
  integer(kind=ik),dimension(:), intent(in)                              :: vtxdist
  type(cel_sp_mat_distr_vec_type)  , intent(inout)                       :: vector_x
  type(cel_sp_mat_distr_vec_type) , intent(inout)                        :: vector_r
  type(cel_sp_mat_distr_vec_type),dimension(:),allocatable,intent(inout) :: vector_v
  type(cel_sp_mat_distr_vec_type),dimension(:),allocatable,intent(inout) :: vector_z
  type(cel_sp_mat_distr_vec_type) ,dimension(:),allocatable,intent(inout):: vector_w
  type(cel_sp_mat_distr_vec_type) , intent(inout)                        :: vector_tmp
  real(kind=rk), allocatable, dimension(:,:), intent(inout)              :: cel_gmres_hh
  real(kind=rk), allocatable, dimension(:)  , intent(inout)              :: cel_gmres_cc
  real(kind=rk), allocatable, dimension(:)  , intent(inout)              :: cel_gmres_ss
  real(kind=rk), allocatable, dimension(:)  , intent(inout)              :: cel_gmres_gamma
  real(kind=rk), allocatable, dimension(:)  , intent(inout)              :: cel_gmres_alpha
  type(cel_perf_counter_type)  , intent(inout)                           :: gmres_perf_counter
  type(cel_perf_counter_type)  , intent(inout)                           :: gmres_perf_counter_mv
  type(cel_perf_counter_type)   , intent(inout)                          :: gmres_perf_counter_vec
  type(cel_profile_type),dimension(:),allocatable,intent(inout)          :: gmres_cel_profiles
  type(cel_omp_type), intent(in)                                         :: cel_omp
  type(cel_comm_type), intent(in)                                        :: cel_comm
  type(cel_input_parameters_type), intent(in)                            :: input_parameters
  integer(kind=ik), intent(out)                                          :: err_code
  integer(kind=ik) :: ii, max_j, sys_stat, vec_length
  integer(kind=ik) :: hh_size_dim1, hh_size_dim2, v_size_dim1
  integer(kind=ik) :: cc_size_dim1, ss_size_dim1, gamma_size_dim1, alpha_size_dim1
  err_code = 0_ik
  !length of the vectors local part (a number of rows in the local part of the matrix A)
  vec_length =  vtxdist(cel_comm%get_from(1)+2_ik)- vtxdist(cel_comm%get_from(1)+1_ik)
  !restart size
   max_j = input_parameters%gmres_max_restart
  !hessenbergmatrix size
  hh_size_dim1 = max_j + 1_ik
  hh_size_dim2 = max_j + 1_ik
  !rotation matrix size
  cc_size_dim1 = max_j + 1_ik
  ss_size_dim1 = max_j + 1_ik
  !gamma size
  gamma_size_dim1 = max_j + 1_ik
  !aplha size
  alpha_size_dim1 = max_j + 1_ik
  !num of v vectors
  v_size_dim1 = max_j + 1_ik

  !check restart
  if(max_j .le. 0_ik) then
    call cel_error("Error cel_gmresalg_parameter_new restart size is equal or less to zero ", 1_ik,&
                        .TRUE._lk,.TRUE._lk)
  end if
  if(cel_omp%is_master) then
    !write(*,*) "max_j:",max_j
    !write(*,*) "hh_size_dim1:",hh_size_dim1
    !write(*,*) "hh_size_dim2:",hh_size_dim2
    !write(*,*) "cc_size_dim1:",cc_size_dim1
    !write(*,*) "ss_size_dim1:",ss_size_dim1
    !write(*,*) "gamma_size_dim1:",gamma_size_dim1
    !write(*,*) "alpha_size_dim1:",alpha_size_dim1
    !write(*,*) "v_size_dim1:",v_size_dim1
    !allocate memory
    if(allocated(cel_gmres_hh)) deallocate(cel_gmres_hh)
    allocate(cel_gmres_hh(hh_size_dim1,hh_size_dim2), stat=sys_stat)  
    call cel_error("Error cel_gmresalg_parameter_new alloc hessenbergmatrix", sys_stat,&
                          cel_is_in_debug,.TRUE._lk)
    if(allocated(cel_gmres_cc)) deallocate(cel_gmres_cc)
    allocate(cel_gmres_cc(cc_size_dim1), stat=sys_stat)  
    call cel_error("Error cel_gmresalg_parameter_new alloc rotation matrix the first component", sys_stat,&
                          cel_is_in_debug,.TRUE._lk)
    if(allocated(cel_gmres_ss)) deallocate(cel_gmres_ss)
    allocate(cel_gmres_ss(ss_size_dim1), stat=sys_stat)  
    call cel_error("Error cel_gmresalg_parameter_new alloc rotation matrix the second component", sys_stat,&
                          cel_is_in_debug,.TRUE._lk)
    if(allocated(cel_gmres_gamma)) deallocate(cel_gmres_gamma)
    allocate(cel_gmres_gamma(gamma_size_dim1), stat=sys_stat)  
    call cel_error("Error cel_gmresalg_parameter_new alloc gamma", sys_stat,&
                          cel_is_in_debug,.TRUE._lk)
    if(allocated(cel_gmres_alpha)) deallocate(cel_gmres_alpha)
    allocate(cel_gmres_alpha(alpha_size_dim1), stat=sys_stat)  
    call cel_error("Error cel_gmresalg_parameter_new alloc alpha", sys_stat,&
                          cel_is_in_debug,.TRUE._lk)
    if(allocated(vector_v)) deallocate(vector_v)
    allocate(vector_v(v_size_dim1), stat=sys_stat)  
    call cel_error("Error cel_gmresalg_parameter_new alloc vector_v", sys_stat,&
                          cel_is_in_debug,.TRUE._lk)
    if(allocated(vector_z)) deallocate(vector_z)
    allocate(vector_z(v_size_dim1), stat=sys_stat)  
    call cel_error("Error cel_gmresalg_parameter_new alloc vector_z", sys_stat,&
                          cel_is_in_debug,.TRUE._lk)
    if(allocated(vector_w)) deallocate(vector_w)
      allocate(vector_w(v_size_dim1), stat=sys_stat)  
      call cel_error("Error cel_gmresalg_parameter_new alloc vector_w", sys_stat,&
                            cel_is_in_debug,.TRUE._lk)
  end if
  !$omp barrier                      
  !initialized all involved vectors for the gmres start
  call cel_sp_mat_distr_vec_new_from_procs(vector_x,cel_comm%get_from,vtxdist,cel_omp)
  call cel_sp_mat_distr_vec_new_from_procs(vector_r,cel_comm%get_from,vtxdist,cel_omp)
  
  do ii=1, v_size_dim1
    call cel_sp_mat_distr_vec_new_from_procs(vector_w(ii),cel_comm%get_from,vtxdist,cel_omp)
    call cel_sp_mat_distr_vec_new_from_procs(vector_v(ii),cel_comm%get_from,vtxdist,cel_omp)
    call cel_sp_mat_distr_vec_new_from_procs(vector_z(ii),cel_comm%get_from,vtxdist,cel_omp)
  end do
  
  call cel_sp_mat_distr_vec_new_from_procs(vector_tmp,cel_comm%get_from,vtxdist,cel_omp)
  !set rest of the gmres parameters
  if(cel_omp%is_master) then
    cel_gmres%residuum_norm = 0_rk
    cel_gmres%betta = 0_rk
    cel_gmres%eps = input_parameters%eps
    cel_gmres%restart_curr = 0_ik
    cel_gmres%restart_max = input_parameters%gmres_max_restart
    cel_gmres%verbosity = input_parameters%verbosity
    cel_gmres%proc_out = 0_ik
    cel_gmres%conv_reason = 0_ik
    cel_gmres%iter_max = input_parameters%iter_max
    cel_gmres%iter_done = 0_ik
    call cel_perf_counter_reset(gmres_perf_counter)
    call cel_perf_counter_reset(gmres_perf_counter_mv)
    call cel_perf_counter_reset(gmres_perf_counter_vec)
    vec_length = vtxdist(cel_omp%master_num+2_ik)-vtxdist(cel_omp%master_num+1_ik)
    gmres_perf_counter%size = vec_length
    if(input_parameters%write_profile) then
      if(allocated(gmres_cel_profiles)) deallocate(gmres_cel_profiles)
      allocate(gmres_cel_profiles(max_j), stat=sys_stat)
      call cel_error("Error cel_gmresalg_parameter_new profile allocation ", sys_stat,&
                        cel_is_in_debug,.TRUE._lk)
    end if
  end if
  !$omp barrier
end subroutine cel_gmresalg_parameter_new



subroutine cel_gmresalg_parameter_del(cel_gmres)
  type(cel_gmresalg_parameter_type)                 :: cel_gmres
  
end subroutine cel_gmresalg_parameter_del


end module cel_alg_module
