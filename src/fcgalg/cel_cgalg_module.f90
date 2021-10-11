!     
! File:   cel_cgalg_module.f90
! Project CRESTA (see details on https://www.cresta-project.eu) Exascale library
! Send your email about the bugs to
! The High Performance Computing Center Stuttgart (HLRS) of the University of Stuttgart
! Dmitry Khabi,  (email: khabi@hlrs.de)
! Created on Apr 8, 2014
!

module cel_cgalg_module
use cel_fop_mv_module
use cel_fop_vec_distr_module
use cel_fop_block_jacobi_module
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
use cel_comm_module
use cel_alg_module
use MPI
use OMP_LIB
implicit none
contains

!CG according to Andreas Meister Numerik linearer Gleichungssysteme
subroutine cel_cgalg(cel_cgalg_par,a_sp_mats,vtxdist,boundary_vector,&
   vector_x,vector_r,vector_d,vector_q,vector_tmp,&
   cel_omp,cel_comm,cel_omp_shared_work,&
   cel_perf_counter,cg_perf_counter_mv,cg_perf_counter_vec,cg_cel_profiles,&
   output_on,err_code)
  type(cel_cgalg_parameter_type)                                    :: cel_cgalg_par
  type(cel_sp_mat_type), dimension(:), allocatable,intent(inout)    :: a_sp_mats
  integer(kind=ik),dimension(:), intent(in)                         :: vtxdist
  type(cel_sp_mat_distr_vec_type), intent(inout)                    :: boundary_vector
  type(cel_sp_mat_distr_vec_type)  , intent(inout)                  :: vector_x
  type(cel_sp_mat_distr_vec_type) , intent(inout)                   :: vector_r
  type(cel_sp_mat_distr_vec_type)  , intent(inout)                  :: vector_d
  type(cel_sp_mat_distr_vec_type) , intent(inout)                   :: vector_q
  type(cel_sp_mat_distr_vec_type) , intent(inout)                   :: vector_tmp
  type(cel_omp_type), intent(in)                                    :: cel_omp
  type(cel_comm_type), intent(inout)                                :: cel_comm
  type(cel_omp_shared_work_type), intent(inout)                     :: cel_omp_shared_work
  type(cel_perf_counter_type), intent(inout)                        :: cel_perf_counter
  type(cel_perf_counter_type)  , intent(inout)                      :: cg_perf_counter_mv
  type(cel_perf_counter_type)   , intent(inout)                     :: cg_perf_counter_vec
  type(cel_profile_type),dimension(:),allocatable,intent(inout)     :: cg_cel_profiles
  logical(kind=lk)  , intent(in)                                    :: output_on
  integer(kind=ik), intent(out)                                     :: err_code
  integer(kind=ik) :: ii
  integer(kind=mpik) :: ierr
  logical(kind=lk) :: flag_error
  err_code = 0_ik
  !calculate r0=b-Ax0 in two steps r0=Ax and r0=b-r0 and d0=r0
  call cel_fop_vec_distr_to_zero(vector_r,cel_omp,cel_omp_shared_work,&
                                cel_perf_counter, output_on,err_code)
  call cel_fop_vec_distr_to_constant(vector_x,1.0_rk,cel_omp,cel_omp_shared_work,&
                                cel_perf_counter, output_on,err_code)
 !$omp barrier
  call cel_fop_mv_axy_distr(a_sp_mats,vector_x,vector_r,vector_tmp,&
                             cel_comm,vtxdist,cel_omp,cel_omp_shared_work,&
                             cel_perf_counter, output_on, err_code)
  ! if(cel_omp%is_master .and. cel_omp%master_num .eq. cel_cgalg_par%proc_out ) then
      !write(*,'(A,I0)') "proc:",cel_omp%master_num
      !call print("array:",vector_r%values(:,1))
  !end if

  
 !norm2=|r|
  call cel_fop_vec_distr_norm2_allreduce(cel_cgalg_par%norm2,vector_r, &
    cel_omp,cel_omp_shared_work,&
    cel_perf_counter, output_on,err_code)
 

  call cel_fop_vec_distr_calc_r_and_d(vector_r,vector_d,&
                                boundary_vector,&
                                cel_omp,cel_omp_shared_work,&
                                cel_perf_counter, output_on,err_code)

  !$omp barrier
  !sigma_new=r0^T * r0 and sigma_old = sigma_new and alpha = 1.0
  call cel_fop_vec_distr_sigma_new_1_allreduce(cel_cgalg_par%sigma_new,cel_cgalg_par%sigma_old, &
                                cel_cgalg_par%alpha, vector_r, &
                                cel_omp,cel_omp_shared_work,&
                                cel_perf_counter, output_on,err_code)

   if(output_on) then
     if(cel_omp%is_master .and. cel_omp%master_num .eq. cel_cgalg_par%proc_out ) then
      write(cel_output_unit,'(5(A,E13.6),2(A,I0))') "sigma_new:",cel_cgalg_par%sigma_new,&
        "; sigma_old:",cel_cgalg_par%sigma_old,&
        "; alpha:",cel_cgalg_par%alpha,&
        "; eps:",cel_cgalg_par%eps,&
        "; d^T*q:", cel_cgalg_par%scalar_d_q,&
        "; iter:", 0_ik,&
        "; iter_max:", cel_cgalg_par%iter_max
     end if
   end if
   ii = 0_ik
   do while(ii<cel_cgalg_par%iter_max .and. cel_cgalg_par%sigma_new .GT. &
    cel_cgalg_par%sigma_old*cel_cgalg_par%eps*cel_cgalg_par%eps .and. &
    abs(cel_cgalg_par%sigma_new) .GT. cel_real_small)
    
     !q_i=A*d_i
     call cel_fop_vec_distr_to_zero(vector_q,cel_omp,cel_omp_shared_work,&
                                cel_perf_counter, output_on,err_code)
     !$omp barrier
     call cel_fop_mv_axy_distr(a_sp_mats,vector_d,vector_q,vector_tmp,&
                             cel_comm,vtxdist,cel_omp,cel_omp_shared_work,&
                             cel_perf_counter, output_on, err_code)
     ii=ii+1_ik
     !alpha=sigma_new/d_i^T*q_i
     call cel_fop_vec_distr_alpha_allreduce(cel_cgalg_par%alpha,cel_cgalg_par%sigma_new,vector_d, &
                                vector_q,cel_omp,cel_omp_shared_work,&
                                cel_perf_counter, output_on,err_code)
     !x_(i+1)<-x_(i)+alpha*d_(i)
     call cel_fop_vec_distr_scale(vector_x,cel_cgalg_par%alpha, vector_d,&
                                cel_omp,cel_omp_shared_work,&
                                cel_perf_counter, output_on,err_code)
     !r_(i+1)=r_(i)-alpha*q_i or r_(i+1)=b-A*x_(i)
     if (.true.) then !MOD(ii+1,10) .GT. 0) then
       !r_(i+1)=r_(i)-alpha*q_i 
       call cel_fop_vec_distr_scale(vector_r,-cel_cgalg_par%alpha, vector_q,&
                                cel_omp,cel_omp_shared_work,&
                                cel_perf_counter, output_on,err_code)
     else
        call cel_fop_vec_distr_to_zero(vector_r,cel_omp,cel_omp_shared_work,&
                                cel_perf_counter, output_on,err_code)
        !r_(i+1)=A*x_(i) and r_(i+1) = b-r_(i+1)
        call cel_fop_mv_axy_distr(a_sp_mats,vector_x,vector_r,vector_tmp,&
                             cel_comm,vtxdist,cel_omp,cel_omp_shared_work,&
                             cel_perf_counter, output_on, err_code)

        call cel_fop_vec_distr_calc_r(vector_r,boundary_vector,&
                                cel_omp,cel_omp_shared_work,&
                                cel_perf_counter, output_on,err_code)
        !$omp barrier
     end if
     !sigma_old = sigma_new and sigma_new=r0^T * r0 
     call cel_fop_vec_distr_sigma_new_2_allreduce(cel_cgalg_par%sigma_new,cel_cgalg_par%sigma_old,vector_r, &
                                cel_omp,cel_omp_shared_work,&
                                cel_perf_counter, output_on,err_code)

     !d_(i+1)=betta*d_(i)+r_(i+1)
    call cel_fop_vec_distr_scale_1(vector_d,cel_cgalg_par%alpha, vector_r,&
                                cel_omp,cel_omp_shared_work,&
                                cel_perf_counter, output_on,err_code)
     !norm2=|r|
      call cel_fop_vec_distr_norm2_allreduce(cel_cgalg_par%norm2,vector_r, &
        cel_omp,cel_omp_shared_work,&
        cel_perf_counter, output_on,err_code)
 
      if(output_on) then
       if(cel_omp%is_master .and. cel_omp%master_num .eq. cel_cgalg_par%proc_out ) then
        write(cel_output_unit,'(6(A,E12.6),2(A,I0))') "sigma_new:",cel_cgalg_par%sigma_new,&
          "; sigma_old:",cel_cgalg_par%sigma_old,&
          "; alpha:",cel_cgalg_par%alpha,&
          "; eps:",cel_cgalg_par%eps,&
          "; d^T*q:", cel_cgalg_par%scalar_d_q,&
          "; |r|:", cel_cgalg_par%norm2,&
          "; iter:", ii,&
          "; iter_max:", cel_cgalg_par%iter_max
       end if
      end if
    end do
    if(cel_omp%is_master) then
      cel_cgalg_par%iter_done = ii
    end if
end subroutine cel_cgalg

!CG according to Andreas Meister Numerik linearer Gleichungssysteme
subroutine cel_cgalg_jacobi(cel_cgalg_par,a_sp_mats,vtxdist,boundary_vector,&
   vector_x,vector_r,vector_z,vector_p,vector_v,vector_tmp,&
   cel_omp,cel_comm,cel_omp_shared_work,&
   cel_perf_counter,cg_perf_counter_mv,cg_perf_counter_vec,cg_cel_profiles,&
   output_on,err_code)
  type(cel_cgalg_parameter_type)                                       :: cel_cgalg_par
  type(cel_sp_mat_type), dimension(:), allocatable,intent(inout)    :: a_sp_mats
  integer(kind=ik),dimension(:), intent(in)                         :: vtxdist
  type(cel_sp_mat_distr_vec_type), intent(inout)                    :: boundary_vector
  type(cel_sp_mat_distr_vec_type)  , intent(inout)                  :: vector_x
  type(cel_sp_mat_distr_vec_type) , intent(inout)                   :: vector_r
  type(cel_sp_mat_distr_vec_type) , intent(inout)                   :: vector_z
  type(cel_sp_mat_distr_vec_type)  , intent(inout)                  :: vector_p
  type(cel_sp_mat_distr_vec_type) , intent(inout)                   :: vector_v
  type(cel_sp_mat_distr_vec_type) , intent(inout)                   :: vector_tmp
  type(cel_omp_type), intent(in)                                    :: cel_omp
  type(cel_comm_type), intent(inout)                                :: cel_comm
  type(cel_omp_shared_work_type), intent(inout)                     :: cel_omp_shared_work
  type(cel_perf_counter_type), intent(inout)                        :: cel_perf_counter
  type(cel_perf_counter_type)  , intent(inout)                      :: cg_perf_counter_mv
  type(cel_perf_counter_type)   , intent(inout)                     :: cg_perf_counter_vec
  type(cel_profile_type),dimension(:),allocatable,intent(inout)     :: cg_cel_profiles
  logical(kind=lk)  , intent(in)                                    :: output_on
  integer(kind=ik), intent(out)                                     :: err_code
  integer(kind=ik) :: ii
  integer(kind=mpik) :: ierr
  logical(kind=lk) :: flag_error
  err_code = 0_ik
  !calculate r0=b-Ax0 in two steps r0=Ax and r0=b-r0 and p0=Pr0
  call cel_fop_vec_distr_to_zero(vector_r,cel_omp,cel_omp_shared_work,&
                                cel_perf_counter, output_on,err_code)
  call cel_fop_vec_distr_to_constant(vector_x,1.0_rk,cel_omp,cel_omp_shared_work,&
                                cel_perf_counter, output_on,err_code)
 !$omp barrier
  call cel_fop_mv_axy_distr(a_sp_mats,vector_x,vector_r,vector_tmp,&
                             cel_comm,vtxdist,cel_omp,cel_omp_shared_work,&
                             cel_perf_counter, output_on, err_code)

  call cel_fop_vec_distr_calc_r(vector_r,boundary_vector,&
                                cel_omp,cel_omp_shared_work,&
                                cel_perf_counter, output_on,err_code)
  !$omp barrier
  call cel_fop_vec_mult_distr_jacobi(vector_p,vector_r,&
                      a_sp_mats(1_ik)%jacobi,cel_omp,&
                      cel_omp_shared_work,&
                      cel_perf_counter, output_on,err_code)
  !call cel_fop_block_jacobi_vec_mult(vector_p,vector_r,a_sp_mats(1_ik)%diagonal_inv_blocks,&
  !                              a_sp_mats(1_ik)%num_last_diag_block_rows,&
  !                              cel_omp,cel_omp_shared_work,&
  !                              cel_perf_counter, output_on,err_code)
  !$omp barrier
  !alpha=(r*p)_2
  call cel_fop_vec_distr_dot_allreduce(cel_cgalg_par%alpha,vector_r, &
                                vector_p,cel_omp,cel_omp_shared_work,&
                                cel_perf_counter, output_on,err_code)
  !norm2=||r||
  call cel_fop_vec_distr_norm2_allreduce(cel_cgalg_par%norm2,vector_r, &
    cel_omp,cel_omp_shared_work,&
    cel_perf_counter, output_on,err_code)

   if(output_on) then
     if(cel_omp%is_master .and. cel_omp%master_num .eq. cel_cgalg_par%proc_out ) then
      write(*,'(2(A,E13.6),2(A,I0))') "alpha:",cel_cgalg_par%alpha,&
        "; |r|:",cel_cgalg_par%norm2,&
        "; iter:", 0_ik,&
        "; iter_max:", cel_cgalg_par%iter_max
     end if
   end if
   ii = 0_ik
   do while(ii<cel_cgalg_par%iter_max .and. &
      abs(cel_cgalg_par%alpha) .gt. cel_real_small .and. &
      cel_cgalg_par%norm2 .gt. cel_cgalg_par%eps)
    !v_i=A*p_i
     call cel_fop_vec_distr_to_zero(vector_v,cel_omp,cel_omp_shared_work,&
                                cel_perf_counter, output_on,err_code)
     !$omp barrier
     call cel_fop_mv_axy_distr(a_sp_mats,vector_p,vector_v,vector_tmp,&
                             cel_comm,vtxdist,cel_omp,cel_omp_shared_work,&
                             cel_perf_counter, output_on, err_code)
     ii=ii+1_ik
     !sigma_new=alpha/(v_i*p_i)_2
     call cel_fop_vec_distr_alpha_allreduce(cel_cgalg_par%sigma_new,cel_cgalg_par%alpha,vector_v, &
                                vector_p,cel_omp,cel_omp_shared_work,&
                                cel_perf_counter, output_on,err_code)
     !x(i+1)<-x_(i)+sigma_new*p_(i)
     call cel_fop_vec_distr_scale(vector_x,cel_cgalg_par%sigma_new, vector_p,&
                                cel_omp,cel_omp_shared_work,&
                                cel_perf_counter, output_on,err_code)
     !r(i+1)=r_(i)-sigma_new*v(i)
     call cel_fop_vec_distr_scale(vector_r,-cel_cgalg_par%sigma_new, vector_v,&
                                cel_omp,cel_omp_shared_work,&
                                cel_perf_counter, output_on,err_code)
     !z(i+1)=P*r(i+1)
     call cel_fop_vec_mult_distr_jacobi(vector_z,vector_r,&
                                a_sp_mats(1_ik)%jacobi,&
                                cel_omp,cel_omp_shared_work,&
                                cel_perf_counter, output_on,err_code)

     !sigma=(r*z)_2/alpha; alpha=(r*z)_2
     call cel_fop_vec_distr_alpha_sigma_r_z_alpha_allreduce(cel_cgalg_par%alpha,&
                                cel_cgalg_par%sigma_new,vector_r, &
                                vector_z,cel_omp,cel_omp_shared_work,&
                                cel_perf_counter, output_on,err_code)
     !p(i+1)=sigma*p(i)+z
     call cel_fop_vec_distr_scale_1(vector_p,cel_cgalg_par%sigma_new, vector_z,&
                                cel_omp,cel_omp_shared_work,&
                                cel_perf_counter, output_on,err_code)
   
     !norm2=|r|
      call cel_fop_vec_distr_norm2_allreduce(cel_cgalg_par%norm2,vector_r, &
        cel_omp,cel_omp_shared_work,&
        cel_perf_counter, output_on,err_code)
 
      if(output_on) then
       if(cel_omp%is_master .and. cel_omp%master_num .eq. cel_cgalg_par%proc_out ) then
        write(cel_output_unit,'(4(A,E12.6),2(A,I0))') "sigma_new:",cel_cgalg_par%sigma_new,&
          "; alpha:",cel_cgalg_par%alpha,&
          "; |r|:", cel_cgalg_par%norm2,&
          "; eps:",cel_cgalg_par%eps,&
          "; iter:", ii,&
          "; iter_max:", cel_cgalg_par%iter_max
       end if
      end if
    end do
    if(cel_omp%is_master) then
      cel_cgalg_par%iter_done = ii
    end if
end subroutine cel_cgalg_jacobi

end module cel_cgalg_module
