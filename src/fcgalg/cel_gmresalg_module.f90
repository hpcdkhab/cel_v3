!     
! File:   cel_cgalg_module.f90
! Project CRESTA (see details on https://www.cresta-project.eu) Exascale library
! Send your email about the bugs to
! The High Performance Computing Center Stuttgart (HLRS) of the University of Stuttgart
! Dmitry Khabi,  (email: khabi@hlrs.de)
! Created on Apr 8, 2014
!

module cel_gmresalg_module
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

!GMRES according to Andreas Meister Numerik linearer Gleichungssysteme
!Use of globale variables cel_gmres, the other version _2 see below
subroutine cel_gmresalg(cel_gmres,a_sp_mats,vtxdist,boundary_vector,&
   vector_x,vector_r,vector_v,vector_w,vector_tmp,&
   cel_gmres_hh,cel_gmres_cc,cel_gmres_ss,cel_gmres_gamma,cel_gmres_alpha,&
   cel_omp,cel_comm,cel_omp_shared_work,&
   cel_perf_counter,cg_perf_counter_mv,cg_perf_counter_vec,cg_cel_profiles,&
   output_on,err_code)
  type(cel_gmresalg_parameter_type), intent(out)                           :: cel_gmres
  type(cel_sp_mat_type), dimension(:), allocatable,intent(inout)           :: a_sp_mats
  integer(kind=ik),dimension(:), intent(in)                                :: vtxdist
  type(cel_sp_mat_distr_vec_type), intent(inout)                           :: boundary_vector
  type(cel_sp_mat_distr_vec_type)  , intent(inout)                         :: vector_x
  type(cel_sp_mat_distr_vec_type) , intent(inout)                          :: vector_r
  type(cel_sp_mat_distr_vec_type), dimension(:), allocatable,intent(inout) :: vector_v
  type(cel_sp_mat_distr_vec_type), dimension(:), allocatable,intent(inout) :: vector_w
  type(cel_sp_mat_distr_vec_type) , intent(inout)                          :: vector_tmp
  real(kind=rk), allocatable, dimension(:,:), intent(inout)                :: cel_gmres_hh
  real(kind=rk), allocatable, dimension(:)  , intent(inout)                :: cel_gmres_cc
  real(kind=rk), allocatable, dimension(:)  , intent(inout)                :: cel_gmres_ss
  real(kind=rk), allocatable, dimension(:)  , intent(inout)                :: cel_gmres_gamma
  real(kind=rk), allocatable, dimension(:)  , intent(inout)                :: cel_gmres_alpha
  type(cel_omp_type), intent(in)                                           :: cel_omp
  type(cel_comm_type), intent(inout)                                       :: cel_comm
  type(cel_omp_shared_work_type), intent(inout)                            :: cel_omp_shared_work
  type(cel_perf_counter_type), intent(inout)                               :: cel_perf_counter
  type(cel_perf_counter_type)  , intent(inout)                             :: cg_perf_counter_mv
  type(cel_perf_counter_type)   , intent(inout)                            :: cg_perf_counter_vec
  type(cel_profile_type),dimension(:),allocatable,intent(inout)            :: cg_cel_profiles
  logical(kind=lk)  , intent(in)                                           :: output_on
  integer(kind=ik), intent(out)                                            :: err_code
  integer(kind=ik) :: ii,jj,kk, vec_length, restart_max, restart_count, conv_reason
  integer(kind=mpik) :: ierr, iter_max, iter_count_local
  logical(kind=lk) :: flag_stop
  integer(kind=ik) :: start_idx,end_idx,omp_start_idx,omp_end_idx,num_threads
  real(kind=rk) :: temp_double
  err_code = 0_ik
  conv_reason = GMRES_REASON_UNDEX
  restart_max = cel_gmres%restart_max
  iter_max = cel_gmres%iter_max
  vec_length = size(vector_x%values,dim=1)
  iter_count_local = 0_ik
  if(.true.) then
    if(cel_omp%is_master .and. cel_omp%master_num .eq. 0_ik) then
      write(cel_output_unit,'(2(A,I0),1(A,E12.6))') &
       "(G var) GMRES parameters: iter_max:",iter_max,&
       "; vec_length:", vec_length,&
       "; eps:",cel_gmres%eps
    end if
  end if
  call cel_fop_vec_distr_to_constant(vector_x,1.0_rk,cel_omp,cel_omp_shared_work,&
                                 cel_perf_counter, output_on,err_code)  
  do while(iter_count_local .lt. iter_max .and. conv_reason .eq. GMRES_REASON_UNDEX)
    iter_count_local = iter_count_local + 1_ik 
    !calculate r0=b-Ax0 in two steps r0=Ax and r0=b-r0 and d0=r0
    call cel_fop_vec_distr_to_zero(vector_r,cel_omp,cel_omp_shared_work,&
                                 cel_perf_counter, output_on,err_code)

    !$omp barrier

    call cel_fop_mv_axy_distr(a_sp_mats,vector_x,vector_r,vector_tmp,&
                                 cel_comm,vtxdist,cel_omp,cel_omp_shared_work,&
                                 cel_perf_counter, output_on, err_code)
    
    !vector_r = boundary_vector - vector_r
    call cel_fop_vec_distr_calc_r(vector_r,boundary_vector,&
                                  cel_omp,cel_omp_shared_work,&
                                  cel_perf_counter, output_on,err_code)
    !norm2 = |r|
    call cel_fop_vec_distr_norm2_start_gmres_allreduce(cel_gmres%residuum_norm,&
                                  cel_gmres_gamma(1_ik), cel_gmres%conv_reason,&
                                  cel_real_small, vector_r, &
                                  cel_omp,cel_omp_shared_work,&
                                  cel_perf_counter, output_on,err_code)

    restart_count = 0_ik
  
    do while(restart_count .le. restart_max .and. conv_reason .eq. GMRES_REASON_UNDEX)
        !======v_1=r*1/norm2======
        
        call cel_fop_vec_distr_scale_inv(vector_v(1_ik),vector_r,cel_gmres%residuum_norm,&
                                    cel_omp,cel_omp_shared_work,&
                                    cel_perf_counter, output_on,err_code)
        !==================
        if(cel_omp%is_master) then
          cel_gmres_gamma(1_ik) = cel_gmres%residuum_norm
        end if
        !$omp barrier
        do jj=1,restart_max
          !======w_j=Av_i <- temporaly result in w======
          call cel_fop_vec_distr_to_zero(vector_w(jj),cel_omp,cel_omp_shared_work,&
                                   cel_perf_counter, output_on,err_code)
          !$omp barrier
          call cel_fop_mv_axy_distr(a_sp_mats,vector_v(jj),&
                                  vector_w(jj),vector_tmp,&
                                   cel_comm,vtxdist,cel_omp,cel_omp_shared_work,&
                                   cel_perf_counter, output_on, err_code)
          !==================
          do ii=1,jj
            !======h_i,j=(v_i,w_j)_2======
            call cel_fop_vec_distr_dot_allreduce(cel_gmres_hh(ii,jj),vector_v(ii), &
                                     vector_w(jj),&
                                     cel_omp,cel_omp_shared_work,&
                                     cel_perf_counter,output_on,err_code)
            !==================
            !======w=w-h_i,j*v_i further calculation of w======
            call cel_fop_vec_distr_scale_minus(vector_w(jj),cel_gmres_hh(ii,jj), vector_v(ii),&
                                    cel_omp,cel_omp_shared_work,&
                                    cel_perf_counter, output_on,err_code)
            !==================
          end do
          !===h_j+1,j = ||w_j||_2===
          call cel_fop_vec_distr_norm2_allreduce(cel_gmres_hh(jj+1_ik,jj),vector_w(jj), &
                                  cel_omp,cel_omp_shared_work,&
                                  cel_perf_counter,output_on,err_code)
          !$omp barrier
          !==============================================================
          
          !======transformiere j-te Spalte der Hessenbergmatrix mit======
          !======akkumulierten Givensrotationen==========================
          call cel_omp_shared_work_start(cel_omp_shared_work,cel_omp%num_threads)
          if(cel_omp%worker_num .eq. 1) then
            do ii=1_ik, jj-1_ik
              cel_gmres_hh(ii,jj) = cel_gmres_hh(ii,jj)*cel_gmres_cc(ii)+&
              cel_gmres_hh(ii+1,jj)*cel_gmres_ss(ii)
              cel_gmres_hh(ii,jj) = -cel_gmres_hh(ii,jj)*cel_gmres_ss(ii)+&
              cel_gmres_hh(ii+1,jj)*cel_gmres_cc(ii)
            end do
            !==================        
            !======calculate further coeff.======
            cel_gmres%betta = sqrt(cel_gmres_hh(jj,jj)**2_ik+cel_gmres_hh(jj+1_ik,jj)**2_ik)
            cel_gmres_ss(jj) = cel_gmres_hh(jj+1,jj) / cel_gmres%betta
            cel_gmres_cc(jj) = cel_gmres_hh(jj,jj) / cel_gmres%betta
            cel_gmres_hh(jj,jj) = cel_gmres%betta
            cel_gmres_gamma(jj+1_ik) = -cel_gmres_ss(jj)*cel_gmres_gamma(jj)
            cel_gmres_gamma(jj) = cel_gmres_cc(jj)*cel_gmres_gamma(jj)
            !==================
          endif
          call cel_omp_shared_work_end_and_wait(cel_omp_shared_work,cel_omp,LOOP_SLEEP_DEFLENGTH)
        

          if(abs(cel_gmres_gamma(jj+1_ik)) .le. cel_gmres%eps) then
            !======calculate solution - vector x======
            call cel_omp_shared_work_start(cel_omp_shared_work,cel_omp%num_threads) 
            call cel_omp_shared_work_distr(cel_omp_shared_work,&
                             omp_start_idx,omp_end_idx,num_threads,&
                             jj,jj,&
                             cel_omp,output_on,err_code)
            if(cel_omp%worker_num .le. num_threads) then
              do ii = omp_end_idx,omp_start_idx,-1
                cel_gmres_alpha(ii) = 0.0_rk
                do kk = ii+1_ik,jj
                  cel_gmres_alpha(ii) = cel_gmres_alpha(ii) + &
                   cel_gmres_hh(ii,kk)*cel_gmres_alpha(kk)
                end do
                cel_gmres_alpha(ii) = (1.0_rk/cel_gmres_hh(ii,ii))*&
                (cel_gmres_gamma(ii) - cel_gmres_alpha(ii))
              end do
            end if
            call cel_omp_shared_work_end_and_wait(cel_omp_shared_work,cel_omp,LOOP_SLEEP_DEFLENGTH)
            
            call cel_omp_shared_work_start(cel_omp_shared_work,cel_omp%num_threads)
            call cel_omp_shared_work_distr(cel_omp_shared_work,&
                             omp_start_idx,omp_end_idx,num_threads,&
                             vec_length,vec_length,&
                             cel_omp,output_on,err_code)
            if(cel_omp%worker_num .le. num_threads) then
              do ii = 1_ik,jj
                temp_double = cel_gmres_alpha(ii)
                vector_x%values(omp_start_idx:omp_end_idx,1_ik) =  &
                 vector_x%values(omp_start_idx:omp_end_idx,1_ik) + temp_double*&
                 vector_v(ii)%values(omp_start_idx:omp_end_idx,1_ik)
              end do
              if(cel_omp%worker_num .eq. 1) cel_gmres%conv_reason = GMRES_REASON_XN
              !==================  
            end if
            call cel_omp_shared_work_end_and_wait(cel_omp_shared_work,cel_omp,LOOP_SLEEP_DEFLENGTH)
            conv_reason = GMRES_REASON_XN
            exit
          else !if(abs(cel_gmres_gamma(jj+1_ik)) .le. cel_gmres%eps)
            call cel_omp_shared_work_start(cel_omp_shared_work,cel_omp%num_threads)
            call cel_omp_shared_work_distr(cel_omp_shared_work,&
                             omp_start_idx,omp_end_idx,num_threads,&
                             vec_length,vec_length,&
                             cel_omp,output_on,err_code)
             if(cel_omp%worker_num .le. num_threads) then
               temp_double = cel_gmres_hh(jj+1,jj)
               vector_v(jj+1_ik)%values(omp_start_idx:omp_end_idx,1_ik) = &
                 vector_w(jj)%values(omp_start_idx:omp_end_idx,1_ik) / &
                 temp_double
             end if
             call cel_omp_shared_work_end_and_wait(cel_omp_shared_work,cel_omp,LOOP_SLEEP_DEFLENGTH)
          end if
        end do
        !$omp barrier
        if(conv_reason .eq. GMRES_REASON_UNDEX) then
          call cel_omp_shared_work_start(cel_omp_shared_work,cel_omp%num_threads)
          if(cel_omp%worker_num .eq. 1) then
            do ii = restart_max,1_ik,-1
              cel_gmres_alpha(ii) = 0.0_rk
              do kk = ii+1_ik,restart_max
                cel_gmres_alpha(ii) = cel_gmres_alpha(ii) + &
                   cel_gmres_hh(ii,kk)*cel_gmres_alpha(kk)
              end do
              cel_gmres_alpha(ii) = (1.0_rk/cel_gmres_hh(ii,ii))*&
                (cel_gmres_gamma(ii) - cel_gmres_alpha(ii))
            end do
            do ii = 1_ik,restart_max
              vector_x%values(1_ik:vec_length,1_ik) =  &
               vector_x%values(1_ik:vec_length,1_ik) + cel_gmres_alpha(ii)*&
               vector_v(ii)%values(1_ik:vec_length,1_ik)
            end do
          end if
        end if
        call cel_omp_shared_work_end_and_wait(cel_omp_shared_work,cel_omp,LOOP_SLEEP_DEFLENGTH)
        !======calculate r0=b-Ax0 in two steps r0=Ax and r0=b-r0 and d0=r0======
        call cel_fop_vec_distr_to_zero(vector_r,cel_omp,cel_omp_shared_work,&
                                     cel_perf_counter, output_on,err_code)
        !$omp barrier
        call cel_fop_mv_axy_distr(a_sp_mats,vector_x,vector_r,vector_tmp,&
                                     cel_comm,vtxdist,cel_omp,cel_omp_shared_work,&
                                     cel_perf_counter, output_on, err_code)
        !==================
        !======vector_r = boundary_vector - vector_r======
        call cel_fop_vec_distr_calc_r(vector_r,boundary_vector,&
                                      cel_omp,cel_omp_shared_work,&
                                      cel_perf_counter, output_on,err_code)
        !==================
        !======norm2 = |r|======
        call cel_fop_vec_distr_norm2_allreduce(cel_gmres%residuum_norm,vector_r, &
                                      cel_omp,cel_omp_shared_work,&
                                      cel_perf_counter, output_on,err_code)
        !$omp barrier
        !==================
        if(output_on) then
          if(cel_omp%is_master .and. cel_omp%master_num .eq. 0_ik) then
            write(cel_output_unit,'(2(A,I0),2(A,E12.6))') &
             "; iter_count_local:",iter_count_local,&
             "; restart_count:", restart_count,&
             "; |r|:", cel_gmres%residuum_norm,&
             "; gamma:", cel_gmres_gamma(jj)
          end if
        end if
        restart_count = restart_count+ 1_ik
    end do
    
  end do
  if(output_on) then
    if(cel_omp%is_master .and. cel_omp%master_num .eq. 0_ik) then
      write(cel_output_unit,'(3(A,I0),1(A,E12.6))') &
            "cel_gmres%conv_reason:", cel_gmres%conv_reason,&
            "; iter_count_local:",iter_count_local,&
            "; iter_max:",iter_max,&
            "; |r|:", cel_gmres%residuum_norm
    end if
  end if

end subroutine cel_gmresalg



!GMRES according to Andreas Meister Numerik linearer Gleichungssysteme
!Use of local variables local_gmres insteed of gloabla like in the first version (see above)
subroutine cel_gmresalg_2(cel_gmres,a_sp_mats,vtxdist,boundary_vector,&
   vector_x,vector_r,vector_v,vector_w,vector_tmp,&
   cel_gmres_hh,cel_gmres_cc,cel_gmres_ss,cel_gmres_gamma,cel_gmres_alpha,&
   cel_omp,cel_comm,cel_omp_shared_work,&
   cel_perf_counter,cg_perf_counter_mv,cg_perf_counter_vec,cg_cel_profiles,&
   output_on,err_code)
  type(cel_gmresalg_parameter_type), intent(out)                           :: cel_gmres
  type(cel_sp_mat_type), dimension(:), allocatable,intent(inout)           :: a_sp_mats
  integer(kind=ik),dimension(:), intent(in)                                :: vtxdist
  type(cel_sp_mat_distr_vec_type), intent(inout)                           :: boundary_vector
  type(cel_sp_mat_distr_vec_type)  , intent(inout)                         :: vector_x
  type(cel_sp_mat_distr_vec_type) , intent(inout)                          :: vector_r
  type(cel_sp_mat_distr_vec_type), dimension(:), allocatable,intent(inout) :: vector_v
  type(cel_sp_mat_distr_vec_type), dimension(:), allocatable,intent(inout) :: vector_w
  type(cel_sp_mat_distr_vec_type) , intent(inout)                          :: vector_tmp
  real(kind=rk), allocatable, dimension(:,:), intent(inout)                :: cel_gmres_hh
  real(kind=rk), allocatable, dimension(:)  , intent(inout)                :: cel_gmres_cc
  real(kind=rk), allocatable, dimension(:)  , intent(inout)                :: cel_gmres_ss
  real(kind=rk), allocatable, dimension(:)  , intent(inout)                :: cel_gmres_gamma
  real(kind=rk), allocatable, dimension(:)  , intent(inout)                :: cel_gmres_alpha
  type(cel_omp_type), intent(in)                                           :: cel_omp
  type(cel_comm_type), intent(inout)                                       :: cel_comm
  type(cel_omp_shared_work_type), intent(inout)                            :: cel_omp_shared_work
  type(cel_perf_counter_type), intent(inout)                               :: cel_perf_counter
  type(cel_perf_counter_type)  , intent(inout)                             :: cg_perf_counter_mv
  type(cel_perf_counter_type)   , intent(inout)                            :: cg_perf_counter_vec
  type(cel_profile_type),dimension(:),allocatable,intent(inout)            :: cg_cel_profiles
  logical(kind=lk)  , intent(in)                                           :: output_on
  integer(kind=ik), intent(out)                                            :: err_code
  integer(kind=ik) :: ii,jj,kk, vec_length, restart_max, restart_count, conv_reason
  integer(kind=mpik) :: ierr, iter_max, iter_count_local
  logical(kind=lk) :: flag_stop
  integer(kind=ik) :: start_idx,end_idx,omp_start_idx,omp_end_idx,num_threads, dim1, dim2
  real(kind=rk) :: temp_double, local_eps, local_residuum
  real(kind=rk), allocatable, dimension(:,:)                :: local_gmres_hh
  real(kind=rk), allocatable, dimension(:)                  :: local_gmres_cc
  real(kind=rk), allocatable, dimension(:)                  :: local_gmres_ss
  real(kind=rk), allocatable, dimension(:)                  :: local_gmres_gamma
  real(kind=rk), allocatable, dimension(:)                  :: local_gmres_alpha
  real(kind=rk)                                             :: local_gmres_betta
  
  dim1 = size(cel_gmres_hh,dim=1)
  dim2 = size(cel_gmres_hh,dim=2)
  allocate(local_gmres_hh(dim1,dim2))
  dim1 = size(cel_gmres_cc,dim=1)
  allocate(local_gmres_cc(dim1))
  dim1 = size(cel_gmres_ss,dim=1)
  allocate(local_gmres_ss(dim1))
  dim1 = size(cel_gmres_gamma,dim=1)
  allocate(local_gmres_gamma(dim1))
  dim1 = size(cel_gmres_alpha,dim=1)
  allocate(local_gmres_alpha(dim1))
  local_eps=cel_gmres%eps
  err_code = 0_ik
  conv_reason = GMRES_REASON_UNDEX
  restart_max = cel_gmres%restart_max
  iter_max = cel_gmres%iter_max
  vec_length = size(vector_x%values,dim=1)
  iter_count_local = 0_ik
  if(output_on) then
    if(cel_omp%is_master .and. cel_omp%master_num .eq. 0_ik) then
      write(cel_output_unit,'(2(A,I0),1(A,E12.6))') &
       "GMRES parameters: iter_max:",iter_max,&
       "; vec_length:", vec_length,&
       "; eps:",cel_gmres%eps
    end if
  end if
  call cel_fop_vec_distr_to_constant(vector_x,1.0_rk,cel_omp,cel_omp_shared_work,&
                                 cel_perf_counter, output_on,err_code)  
  do while(iter_count_local .lt. iter_max .and. conv_reason .eq. GMRES_REASON_UNDEX)
    iter_count_local = iter_count_local + 1_ik 
    !calculate r0=b-Ax0 in two steps r0=Ax and r0=b-r0 and d0=r0
    call cel_fop_vec_distr_to_zero(vector_r,cel_omp,cel_omp_shared_work,&
                                 cel_perf_counter, output_on,err_code)

    !$omp barrier

    call cel_fop_mv_axy_distr(a_sp_mats,vector_x,vector_r,vector_tmp,&
                                 cel_comm,vtxdist,cel_omp,cel_omp_shared_work,&
                                 cel_perf_counter, output_on, err_code)
    
    !vector_r = boundary_vector - vector_r
    call cel_fop_vec_distr_calc_r(vector_r,boundary_vector,&
                                  cel_omp,cel_omp_shared_work,&
                                  cel_perf_counter, output_on,err_code)
    !norm2 = |r|
    call cel_fop_vec_distr_norm2_start_gmres_allreduce(cel_gmres%residuum_norm,&
                                  cel_gmres_gamma(1_ik), cel_gmres%conv_reason,&
                                  cel_real_small, vector_r, &
                                  cel_omp,cel_omp_shared_work,&
                                  cel_perf_counter, output_on,err_code)
    local_gmres_gamma(1_ik) = cel_gmres_gamma(1_ik)
    local_residuum = cel_gmres%residuum_norm
    restart_count = 0_ik
  
    do while(restart_count .le. restart_max .and. conv_reason .eq. GMRES_REASON_UNDEX)
        !======v_1=r*1/norm2======
        
        call cel_fop_vec_distr_scale_inv(vector_v(1_ik),vector_r,local_residuum,&
                                    cel_omp,cel_omp_shared_work,&
                                    cel_perf_counter, output_on,err_code)
        !==================
        if(cel_omp%is_master) then
          cel_gmres_gamma(1_ik) = cel_gmres%residuum_norm
        end if
        local_gmres_gamma(1_ik) = cel_gmres%residuum_norm

        do jj=1,restart_max
          !======w_j=Av_i <- temporaly result in w======
          call cel_fop_vec_distr_to_zero(vector_w(jj),cel_omp,cel_omp_shared_work,&
                                   cel_perf_counter, output_on,err_code)
          !$omp barrier
          call cel_fop_mv_axy_distr(a_sp_mats,vector_v(jj),&
                                  vector_w(jj),vector_tmp,&
                                   cel_comm,vtxdist,cel_omp,cel_omp_shared_work,&
                                   cel_perf_counter, output_on, err_code)
          !==================
          do ii=1,jj
            !======h_i,j=(v_i,w_j)_2======
            call cel_fop_vec_distr_dot_allreduce(cel_gmres_hh(ii,jj),vector_v(ii), &
                                     vector_w(jj),&
                                     cel_omp,cel_omp_shared_work,&
                                     cel_perf_counter,output_on,err_code)
            local_gmres_hh(ii,jj) = cel_gmres_hh(ii,jj)
            !==================
            !======w=w-h_i,j*v_i further calculation of w======
            call cel_fop_vec_distr_scale_minus(vector_w(jj),local_gmres_hh(ii,jj), vector_v(ii),&
                                    cel_omp,cel_omp_shared_work,&
                                    cel_perf_counter, output_on,err_code)
            !==================
          end do
          !===h_j+1,j = ||w_j||_2===
          call cel_fop_vec_distr_norm2_allreduce(cel_gmres_hh(jj+1_ik,jj),vector_w(jj), &
                                  cel_omp,cel_omp_shared_work,&
                                  cel_perf_counter,output_on,err_code)
          local_gmres_hh(ii,jj) = cel_gmres_hh(jj+1_ik,jj)
          !==================
          !======transformiere j-te Spalte der Hessenbergmatrix mit======
          !======akkumulierten Givensrotationen==========================
          
          do ii=1_ik, jj-1_ik
            local_gmres_hh(ii,jj) = local_gmres_hh(ii,jj)*local_gmres_cc(ii)+&
            local_gmres_hh(ii+1,jj)*local_gmres_ss(ii)
            local_gmres_hh(ii,jj) = -local_gmres_hh(ii,jj)*local_gmres_ss(ii)+&
            local_gmres_hh(ii+1,jj)*local_gmres_cc(ii)
          end do
          local_gmres_betta = sqrt(local_gmres_hh(jj,jj)**2_ik+local_gmres_hh(jj+1_ik,jj)**2_ik)
          local_gmres_ss(jj) = local_gmres_hh(jj+1,jj) / local_gmres_betta
          local_gmres_cc(jj) = local_gmres_hh(jj,jj) / local_gmres_betta
          local_gmres_hh(jj,jj) = local_gmres_betta
          local_gmres_gamma(jj+1_ik) = -local_gmres_ss(jj)*local_gmres_gamma(jj)
          local_gmres_gamma(jj) = local_gmres_cc(jj)*local_gmres_gamma(jj)
          !==============================================================
          

          if(abs(local_gmres_gamma(jj+1_ik)) .le. cel_gmres%eps) then
            !======calculate solution - vector x======
            do ii = jj,1,-1
                local_gmres_alpha(ii) = 0.0_rk
                do kk = ii+1_ik,jj
                  local_gmres_alpha(ii) = local_gmres_alpha(ii) + &
                   local_gmres_hh(ii,kk)*local_gmres_alpha(kk)
                end do
                local_gmres_alpha(ii) = (1.0_rk/local_gmres_hh(ii,ii))*&
                (local_gmres_gamma(ii) - local_gmres_alpha(ii))
              end do
            
            call cel_omp_shared_work_start(cel_omp_shared_work,cel_omp%num_threads)
            call cel_omp_shared_work_distr(cel_omp_shared_work,&
                             omp_start_idx,omp_end_idx,num_threads,&
                             vec_length,vec_length,&
                             cel_omp,output_on,err_code)
            if(cel_omp%worker_num .le. num_threads) then
              do ii = 1_ik,jj
                temp_double = local_gmres_alpha(ii)
                vector_x%values(omp_start_idx:omp_end_idx,1_ik) =  &
                 vector_x%values(omp_start_idx:omp_end_idx,1_ik) + temp_double*&
                 vector_v(ii)%values(omp_start_idx:omp_end_idx,1_ik)
              end do
              if(cel_omp%worker_num .eq. 1) cel_gmres%conv_reason = GMRES_REASON_XN
              !==================  
            end if
            call cel_omp_shared_work_end_and_wait(cel_omp_shared_work,cel_omp,LOOP_SLEEP_DEFLENGTH)
            conv_reason = GMRES_REASON_XN
            exit
          else !if(abs(cel_gmres_gamma(jj+1_ik)) .le. cel_gmres%eps)
            call cel_omp_shared_work_start(cel_omp_shared_work,cel_omp%num_threads)
            call cel_omp_shared_work_distr(cel_omp_shared_work,&
                             omp_start_idx,omp_end_idx,num_threads,&
                             vec_length,vec_length,&
                             cel_omp,output_on,err_code)
             if(cel_omp%worker_num .le. num_threads) then
               temp_double = local_gmres_hh(jj+1,jj)
               vector_v(jj+1_ik)%values(omp_start_idx:omp_end_idx,1_ik) = &
                 vector_w(jj)%values(omp_start_idx:omp_end_idx,1_ik) / &
                 temp_double
             end if
             call cel_omp_shared_work_end_and_wait(cel_omp_shared_work,cel_omp,LOOP_SLEEP_DEFLENGTH)
          end if
        end do

        if(conv_reason .eq. GMRES_REASON_UNDEX) then
          do ii = restart_max,1_ik,-1
              local_gmres_alpha(ii) = 0.0_rk
              do kk = ii+1_ik,restart_max
                local_gmres_alpha(ii) = local_gmres_alpha(ii) + &
                   local_gmres_hh(ii,kk)*local_gmres_alpha(kk)
              end do
              local_gmres_alpha(ii) = (1.0_rk/local_gmres_hh(ii,ii))*&
                (local_gmres_gamma(ii) - local_gmres_alpha(ii))
          end do
          call cel_omp_shared_work_start(cel_omp_shared_work,cel_omp%num_threads)
          call cel_omp_shared_work_distr(cel_omp_shared_work,&
                 omp_start_idx,omp_end_idx,num_threads,&
                 vec_length,vec_length,&
                 cel_omp,output_on,err_code)
          if(cel_omp%worker_num .le. num_threads) then
            do ii = 1_ik,restart_max
              vector_x%values(omp_start_idx:omp_end_idx,1_ik) =  &
               vector_x%values(omp_start_idx:omp_end_idx,1_ik) + local_gmres_alpha(ii)*&
               vector_v(ii)%values(omp_start_idx:omp_end_idx,1_ik)
            end do
          end if
        end if
        call cel_omp_shared_work_end_and_wait(cel_omp_shared_work,cel_omp,LOOP_SLEEP_DEFLENGTH)
        !======calculate r0=b-Ax0 in two steps r0=Ax and r0=b-r0 and d0=r0======
        call cel_fop_vec_distr_to_zero(vector_r,cel_omp,cel_omp_shared_work,&
                                     cel_perf_counter, output_on,err_code)
        !$omp barrier
        call cel_fop_mv_axy_distr(a_sp_mats,vector_x,vector_r,vector_tmp,&
                                     cel_comm,vtxdist,cel_omp,cel_omp_shared_work,&
                                     cel_perf_counter, output_on, err_code)
        !==================
        !======vector_r = boundary_vector - vector_r======
        call cel_fop_vec_distr_calc_r(vector_r,boundary_vector,&
                                      cel_omp,cel_omp_shared_work,&
                                      cel_perf_counter, output_on,err_code)
        !==================
        !======norm2 = |r|======
        call cel_fop_vec_distr_norm2_allreduce(cel_gmres%residuum_norm,vector_r, &
                                      cel_omp,cel_omp_shared_work,&
                                      cel_perf_counter, output_on,err_code)
        local_residuum = cel_gmres%residuum_norm
        !==================
        if(output_on) then
          if(cel_omp%is_master .and. cel_omp%master_num .eq. 0_ik) then
            write(cel_output_unit,'(2(A,I0),2(A,E12.6))') &
             "; iter_count_local:",iter_count_local,&
             "; restart_count:", restart_count,&
             "; |r|:", local_residuum,&
             "; gamma:", local_gmres_gamma(jj)
          end if
        end if
        restart_count = restart_count+ 1_ik
    end do
    
  end do
  if(output_on) then
    if(cel_omp%is_master .and. cel_omp%master_num .eq. 0_ik) then
      write(cel_output_unit,'(3(A,I0),1(A,E12.6))') &
            "cel_gmres%conv_reason:", cel_gmres%conv_reason,&
            "; iter_count_local:",iter_count_local,&
            "; iter_max:",iter_max,&
            "; |r|:", cel_gmres%residuum_norm
    end if
  end if
  if(cel_omp%is_master) then
      cel_gmres%iter_done = iter_count_local
  end if
end subroutine cel_gmresalg_2


!GMRES according to Andreas Meister Numerik linearer Gleichungssysteme
!PGMRES
subroutine cel_gmresalg_2jacobi(cel_gmres,a_sp_mats,vtxdist,boundary_vector,&
   vector_x,vector_r,vector_v,vector_z,vector_w,vector_tmp,&
   cel_gmres_hh,cel_gmres_cc,cel_gmres_ss,cel_gmres_gamma,cel_gmres_alpha,&
   cel_omp,cel_comm,cel_omp_shared_work,&
   cel_perf_counter,cg_perf_counter_mv,cg_perf_counter_vec,cg_cel_profiles,&
   output_on,err_code)
  type(cel_gmresalg_parameter_type), intent(out)                           :: cel_gmres
  type(cel_sp_mat_type), dimension(:), allocatable,intent(inout)           :: a_sp_mats
  integer(kind=ik),dimension(:), intent(in)                                :: vtxdist
  type(cel_sp_mat_distr_vec_type), intent(inout)                           :: boundary_vector
  type(cel_sp_mat_distr_vec_type)  , intent(inout)                         :: vector_x
  type(cel_sp_mat_distr_vec_type) , intent(inout)                          :: vector_r
  type(cel_sp_mat_distr_vec_type), dimension(:), allocatable,intent(inout) :: vector_v
  type(cel_sp_mat_distr_vec_type), dimension(:), allocatable,intent(inout) :: vector_z
  type(cel_sp_mat_distr_vec_type), dimension(:), allocatable,intent(inout) :: vector_w
  type(cel_sp_mat_distr_vec_type) , intent(inout)                          :: vector_tmp
  real(kind=rk), allocatable, dimension(:,:), intent(inout)                :: cel_gmres_hh
  real(kind=rk), allocatable, dimension(:)  , intent(inout)                :: cel_gmres_cc
  real(kind=rk), allocatable, dimension(:)  , intent(inout)                :: cel_gmres_ss
  real(kind=rk), allocatable, dimension(:)  , intent(inout)                :: cel_gmres_gamma
  real(kind=rk), allocatable, dimension(:)  , intent(inout)                :: cel_gmres_alpha
  type(cel_omp_type), intent(in)                                           :: cel_omp
  type(cel_comm_type), intent(inout)                                       :: cel_comm
  type(cel_omp_shared_work_type), intent(inout)                            :: cel_omp_shared_work
  type(cel_perf_counter_type), intent(inout)                               :: cel_perf_counter
  type(cel_perf_counter_type)  , intent(inout)                             :: cg_perf_counter_mv
  type(cel_perf_counter_type)   , intent(inout)                            :: cg_perf_counter_vec
  type(cel_profile_type),dimension(:),allocatable,intent(inout)            :: cg_cel_profiles
  logical(kind=lk)  , intent(in)                                           :: output_on
  integer(kind=ik), intent(out)                                            :: err_code
  integer(kind=ik) :: ii,jj,kk, vec_length, restart_max, restart_count, conv_reason
  integer(kind=mpik) :: ierr, iter_max, iter_count_local
  logical(kind=lk) :: flag_stop
  integer(kind=ik) :: start_idx,end_idx,omp_start_idx,omp_end_idx,num_threads, dim1, dim2
  real(kind=rk) :: temp_double, local_eps, local_residuum
  real(kind=rk), allocatable, dimension(:,:)                :: local_gmres_hh
  real(kind=rk), allocatable, dimension(:)                  :: local_gmres_cc
  real(kind=rk), allocatable, dimension(:)                  :: local_gmres_ss
  real(kind=rk), allocatable, dimension(:)                  :: local_gmres_gamma
  real(kind=rk), allocatable, dimension(:)                  :: local_gmres_alpha
  real(kind=rk)                                             :: local_gmres_betta
  
  dim1 = size(cel_gmres_hh,dim=1)
  dim2 = size(cel_gmres_hh,dim=2)
  allocate(local_gmres_hh(dim1,dim2))
  dim1 = size(cel_gmres_cc,dim=1)
  allocate(local_gmres_cc(dim1))
  dim1 = size(cel_gmres_ss,dim=1)
  allocate(local_gmres_ss(dim1))
  dim1 = size(cel_gmres_gamma,dim=1)
  allocate(local_gmres_gamma(dim1))
  dim1 = size(cel_gmres_alpha,dim=1)
  allocate(local_gmres_alpha(dim1))
  local_eps=cel_gmres%eps
  err_code = 0_ik
  conv_reason = GMRES_REASON_UNDEX
  restart_max = cel_gmres%restart_max
  iter_max = cel_gmres%iter_max
  vec_length = size(vector_x%values,dim=1)
  iter_count_local = 0_ik
  if(output_on) then
    if(cel_omp%is_master .and. cel_omp%master_num .eq. 0_ik) then
      write(cel_output_unit,'(2(A,I0),1(A,E12.6))') &
       "(P)GMRES parameters: iter_max:",iter_max,&
       "; vec_length:", vec_length,&
       "; eps:",cel_gmres%eps
    end if
  end if
  call cel_fop_vec_distr_to_constant(vector_x,1.0_rk,cel_omp,cel_omp_shared_work,&
                                 cel_perf_counter, output_on,err_code)  
  do while(iter_count_local .le. iter_max .and. conv_reason .eq. GMRES_REASON_UNDEX)
    iter_count_local = iter_count_local + 1_ik 
    !calculate r0=b-Ax0 in two steps r0=Ax and r0=b-r0 and d0=r0
    call cel_fop_vec_distr_to_zero(vector_r,cel_omp,cel_omp_shared_work,&
                                 cel_perf_counter, output_on,err_code)

    !$omp barrier

    call cel_fop_mv_axy_distr(a_sp_mats,vector_x,vector_r,vector_tmp,&
                                 cel_comm,vtxdist,cel_omp,cel_omp_shared_work,&
                                 cel_perf_counter, output_on, err_code)
    
    !vector_r = boundary_vector - vector_r
    call cel_fop_vec_distr_calc_r(vector_r,boundary_vector,&
                                  cel_omp,cel_omp_shared_work,&
                                  cel_perf_counter, output_on,err_code)
    !norm2 = |r|
    call cel_fop_vec_distr_norm2_start_gmres_allreduce(cel_gmres%residuum_norm,&
                                  cel_gmres_gamma(1_ik), cel_gmres%conv_reason,&
                                  cel_real_small, vector_r, &
                                  cel_omp,cel_omp_shared_work,&
                                  cel_perf_counter, output_on,err_code)
    local_gmres_gamma(1_ik) = cel_gmres_gamma(1_ik)
    local_residuum = cel_gmres%residuum_norm
    restart_count = 0_ik
  
    do while(restart_count .le. restart_max .and. conv_reason .eq. GMRES_REASON_UNDEX)
        !======v_1=r*1/norm2======
        
        call cel_fop_vec_distr_scale_inv(vector_v(1_ik),vector_r,local_residuum,&
                                    cel_omp,cel_omp_shared_work,&
                                    cel_perf_counter, output_on,err_code)
        !==================
        if(cel_omp%is_master) then
          cel_gmres_gamma(1_ik) = cel_gmres%residuum_norm
        end if
        local_gmres_gamma(1_ik) = cel_gmres%residuum_norm

        do jj=1,restart_max
          !======w_j=A(M^-1)*v_i <- temporaly result in w======
          call cel_fop_vec_distr_to_zero(vector_w(jj),cel_omp,cel_omp_shared_work,&
                                   cel_perf_counter, output_on,err_code)
          !$omp barrier
          call cel_fop_vec_mult_distr_jacobi(vector_z(jj),vector_v(jj),&
                                a_sp_mats(1_ik)%jacobi,cel_omp,&
                                cel_omp_shared_work,&
                                cel_perf_counter, output_on,err_code)
          !call cel_fop_block_jacobi_vec_mult(vector_z(jj),vector_v(jj),&
          !                      a_sp_mats(1_ik)%diagonal_inv_blocks,&
          !                      a_sp_mats(1_ik)%num_last_diag_block_rows,&
          !                      cel_omp,&
          !                      cel_omp_shared_work,&
          !                      cel_perf_counter, output_on,err_code)
          !$omp barrier
          call cel_fop_mv_axy_distr(a_sp_mats,vector_z(jj),&
                                  vector_w(jj),vector_tmp,&
                                   cel_comm,vtxdist,cel_omp,cel_omp_shared_work,&
                                   cel_perf_counter, output_on, err_code)
          !==================
          do ii=1,jj
            !======h_i,j=(v_i,w_j)_2======
            call cel_fop_vec_distr_dot_allreduce(cel_gmres_hh(ii,jj),vector_v(ii), &
                                     vector_w(jj),&
                                     cel_omp,cel_omp_shared_work,&
                                     cel_perf_counter,output_on,err_code)
            local_gmres_hh(ii,jj) = cel_gmres_hh(ii,jj)
            !==================
            !======w=w-h_i,j*v_i further calculation of w======
            call cel_fop_vec_distr_scale_minus(vector_w(jj),local_gmres_hh(ii,jj), vector_v(ii),&
                                    cel_omp,cel_omp_shared_work,&
                                    cel_perf_counter, output_on,err_code)
            !==================
          end do
          !===h_j+1,j = ||w_j||_2===
          call cel_fop_vec_distr_norm2_allreduce(cel_gmres_hh(jj+1_ik,jj),vector_w(jj), &
                                  cel_omp,cel_omp_shared_work,&
                                  cel_perf_counter,output_on,err_code)
          local_gmres_hh(ii,jj) = cel_gmres_hh(jj+1_ik,jj)
          !==================
          !======transformiere j-te Spalte der Hessenbergmatrix mit======
          !======akkumulierten Givensrotationen==========================
          
          do ii=1_ik, jj-1_ik
            local_gmres_hh(ii,jj) = local_gmres_hh(ii,jj)*local_gmres_cc(ii)+&
            local_gmres_hh(ii+1,jj)*local_gmres_ss(ii)
            local_gmres_hh(ii,jj) = -local_gmres_hh(ii,jj)*local_gmres_ss(ii)+&
            local_gmres_hh(ii+1,jj)*local_gmres_cc(ii)
          end do
          local_gmres_betta = sqrt(local_gmres_hh(jj,jj)**2_ik+local_gmres_hh(jj+1_ik,jj)**2_ik)
          local_gmres_ss(jj) = local_gmres_hh(jj+1,jj) / local_gmres_betta
          local_gmres_cc(jj) = local_gmres_hh(jj,jj) / local_gmres_betta
          local_gmres_hh(jj,jj) = local_gmres_betta
          local_gmres_gamma(jj+1_ik) = -local_gmres_ss(jj)*local_gmres_gamma(jj)
          local_gmres_gamma(jj) = local_gmres_cc(jj)*local_gmres_gamma(jj)
          !==============================================================
          

          if(abs(local_gmres_gamma(jj+1_ik)) .le. cel_gmres%eps) then
            !======calculate solution - vector x======
            do ii = jj,1,-1
                local_gmres_alpha(ii) = 0.0_rk
                do kk = ii+1_ik,jj
                  local_gmres_alpha(ii) = local_gmres_alpha(ii) + &
                   local_gmres_hh(ii,kk)*local_gmres_alpha(kk)
                end do
                local_gmres_alpha(ii) = (1.0_rk/local_gmres_hh(ii,ii))*&
                (local_gmres_gamma(ii) - local_gmres_alpha(ii))
              end do
            
            call cel_omp_shared_work_start(cel_omp_shared_work,cel_omp%num_threads)
            call cel_omp_shared_work_distr(cel_omp_shared_work,&
                             omp_start_idx,omp_end_idx,num_threads,&
                             vec_length,vec_length,&
                             cel_omp,output_on,err_code)
            if(cel_omp%worker_num .le. num_threads) then
              do ii = 1_ik,jj
                temp_double = local_gmres_alpha(ii)
                vector_x%values(omp_start_idx:omp_end_idx,1_ik) =  &
                 vector_x%values(omp_start_idx:omp_end_idx,1_ik) + temp_double*&
                 vector_z(ii)%values(omp_start_idx:omp_end_idx,1_ik)
              end do
              if(cel_omp%worker_num .eq. 1) cel_gmres%conv_reason = GMRES_REASON_XN
              !==================  
            end if
            call cel_omp_shared_work_end_and_wait(cel_omp_shared_work,cel_omp,LOOP_SLEEP_DEFLENGTH)
            conv_reason = GMRES_REASON_XN
            exit
          else !if(abs(cel_gmres_gamma(jj+1_ik)) .le. cel_gmres%eps)
            call cel_omp_shared_work_start(cel_omp_shared_work,cel_omp%num_threads)
            call cel_omp_shared_work_distr(cel_omp_shared_work,&
                             omp_start_idx,omp_end_idx,num_threads,&
                             vec_length,vec_length,&
                             cel_omp,output_on,err_code)
             if(cel_omp%worker_num .le. num_threads) then
               temp_double = local_gmres_hh(jj+1,jj)
               vector_v(jj+1_ik)%values(omp_start_idx:omp_end_idx,1_ik) = &
                 vector_w(jj)%values(omp_start_idx:omp_end_idx,1_ik) / &
                 temp_double
             end if
             call cel_omp_shared_work_end_and_wait(cel_omp_shared_work,cel_omp,LOOP_SLEEP_DEFLENGTH)
          end if
        end do

        if(conv_reason .eq. GMRES_REASON_UNDEX) then
          do ii = restart_max,1_ik,-1
              local_gmres_alpha(ii) = 0.0_rk
              do kk = ii+1_ik,restart_max
                local_gmres_alpha(ii) = local_gmres_alpha(ii) + &
                   local_gmres_hh(ii,kk)*local_gmres_alpha(kk)
              end do
              local_gmres_alpha(ii) = (1.0_rk/local_gmres_hh(ii,ii))*&
                (local_gmres_gamma(ii) - local_gmres_alpha(ii))
          end do
          call cel_omp_shared_work_start(cel_omp_shared_work,cel_omp%num_threads)
          call cel_omp_shared_work_distr(cel_omp_shared_work,&
                 omp_start_idx,omp_end_idx,num_threads,&
                 vec_length,vec_length,&
                 cel_omp,output_on,err_code)
          if(cel_omp%worker_num .le. num_threads) then
            do ii = 1_ik,restart_max
              vector_x%values(omp_start_idx:omp_end_idx,1_ik) =  &
               vector_x%values(omp_start_idx:omp_end_idx,1_ik) + local_gmres_alpha(ii)*&
               vector_z(ii)%values(omp_start_idx:omp_end_idx,1_ik)
            end do
          end if
        end if
        call cel_omp_shared_work_end_and_wait(cel_omp_shared_work,cel_omp,LOOP_SLEEP_DEFLENGTH)
        !======calculate r0=b-Ax0 in two steps r0=Ax and r0=b-r0 and d0=r0======
        call cel_fop_vec_distr_to_zero(vector_r,cel_omp,cel_omp_shared_work,&
                                     cel_perf_counter, output_on,err_code)
        !$omp barrier
        call cel_fop_mv_axy_distr(a_sp_mats,vector_x,vector_r,vector_tmp,&
                                     cel_comm,vtxdist,cel_omp,cel_omp_shared_work,&
                                     cel_perf_counter, output_on, err_code)
        !==================
        !======vector_r = boundary_vector - vector_r======
        call cel_fop_vec_distr_calc_r(vector_r,boundary_vector,&
                                      cel_omp,cel_omp_shared_work,&
                                      cel_perf_counter, output_on,err_code)
        !==================
        !======norm2 = |r|======
        call cel_fop_vec_distr_norm2_allreduce(cel_gmres%residuum_norm,vector_r, &
                                      cel_omp,cel_omp_shared_work,&
                                      cel_perf_counter, output_on,err_code)
        local_residuum = cel_gmres%residuum_norm
        !==================
        if(output_on) then
          if(cel_omp%is_master .and. cel_omp%master_num .eq. 0_ik) then
            write(cel_output_unit,'(2(A,I0),2(A,E12.6))') &
             "; iter_count_local:",iter_count_local,&
             "; restart_count:", restart_count,&
             "; |r|:", local_residuum,&
             "; gamma:", local_gmres_gamma(jj)
          end if
        end if
        restart_count = restart_count+ 1_ik
    end do
    
  end do
  if(output_on) then
    if(cel_omp%is_master .and. cel_omp%master_num .eq. 0_ik) then
      write(cel_output_unit,'(3(A,I0),1(A,E12.6))') &
            "cel_gmres%conv_reason:", cel_gmres%conv_reason,&
            "; iter_count_local:",iter_count_local,&
            "; iter_max:",iter_max,&
            "; |r|:", cel_gmres%residuum_norm
    end if
  end if
  if(cel_omp%is_master) then
      cel_gmres%iter_done = iter_count_local
  end if
end subroutine cel_gmresalg_2jacobi


end module cel_gmresalg_module
