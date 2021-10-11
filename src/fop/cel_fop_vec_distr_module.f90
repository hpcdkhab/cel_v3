! Definition of the distributed matrix vector multiplication
! File:   cel_fop_vec_distr_module.f90
! Project CRESTA Exascale library (see details on https://www.cresta-project.eu)
! Send you email about the bugs to
! The High Performance Computing Center Stuttgart (HLRS) of the University of Stuttgart
! Dmitry Khabi,  (email: khabi@hlrs.de)
! Created on Jun 14, 2013
!

module cel_fop_vec_distr_module
use cel_types_module
use cel_timer_interface_module
use cel_cpu_param_module
use cel_sp_mat_distr_vec_module
use cel_sp_mat_distr_module
use cel_sp_mat_distr_gl_module
use cel_perf_module
use cel_base_print_module
use cel_omp_module
use cel_omp_shared_work_module
use MPI
use OMP_LIB
implicit none
contains


!vector_r = boundary_vector - vector_r
!vector_d = vector_r
subroutine cel_fop_vec_distr_calc_r_and_d(vector_r,vector_d,boundary_vector,&
                                cel_omp,cel_omp_shared_work,&
                                cel_perf_counter, output_on,err_code)
  type(cel_sp_mat_distr_vec_type), intent(inout)                    :: vector_r
  type(cel_sp_mat_distr_vec_type), intent(inout)                    :: vector_d
  type(cel_sp_mat_distr_vec_type), intent(in)                       :: boundary_vector
  type(cel_omp_type), intent(in)                                    :: cel_omp
  type(cel_omp_shared_work_type), intent(inout)                     :: cel_omp_shared_work
  type(cel_perf_counter_type), intent(inout)                        :: cel_perf_counter
  logical(kind=lk)  , intent(in)                                    :: output_on
  integer(kind=ik), intent(out)                                     :: err_code
  integer(kind=ik) :: ii, vec_length
  real(kind=trk) :: time
  integer(kind=ik) :: omp_start_idx,omp_end_idx,num_threads
 
  err_code = 0_ik
  if(cel_omp%worker_num .eq. 1_ik) then
   time = get_time()
  end if
  vec_length = size(vector_r%values,dim=1)
  call cel_omp_shared_work_distr(cel_omp_shared_work,&
                                 omp_start_idx,omp_end_idx,num_threads,&
                                 vec_length,vec_length,&
                                 cel_omp,output_on,err_code)
  if(cel_omp%worker_num .le. num_threads) then
      vector_r%values(omp_start_idx:omp_end_idx,1_ik) = &
          boundary_vector%values(omp_start_idx:omp_end_idx,1_ik)-vector_r%values(omp_start_idx:omp_end_idx,1_ik)
      vector_d%values(omp_start_idx:omp_end_idx,1_ik) = vector_r%values(omp_start_idx:omp_end_idx,1_ik)
  end if

  if(cel_omp%worker_num .eq. 1) then
    cel_perf_counter%time(1) = cel_perf_counter%time(1) + (get_time() - time)
    cel_perf_counter%add= cel_perf_counter%add + vec_length
    cel_perf_counter%load= cel_perf_counter%load + 2_ik*vec_length*rk
    cel_perf_counter%store= cel_perf_counter%store + 2_ik*vec_length*rk
  end if
  
end subroutine cel_fop_vec_distr_calc_r_and_d

!vector_v = vector_r / alpha
subroutine cel_fop_vec_distr_scale_inv(vector_v,vector_r,alpha,&
                                cel_omp,cel_omp_shared_work,&
                                cel_perf_counter, output_on,err_code)
  type(cel_sp_mat_distr_vec_type), intent(inout)                    :: vector_v
  type(cel_sp_mat_distr_vec_type), intent(in)                       :: vector_r
  real(kind=rk), intent(in)                                         :: alpha
  type(cel_omp_type), intent(in)                                    :: cel_omp
  type(cel_omp_shared_work_type), intent(inout)                     :: cel_omp_shared_work
  type(cel_perf_counter_type), intent(inout)                        :: cel_perf_counter
  logical(kind=lk)  , intent(in)                                    :: output_on
  integer(kind=ik), intent(out)                                     :: err_code
  integer(kind=ik) :: ii, vec_length
  real(kind=trk) :: time, tmp_real
  integer(kind=ik) :: omp_start_idx,omp_end_idx,num_threads
 
  err_code = 0_ik
  if(cel_omp%worker_num .eq. 1) then
   time = get_time()
  end if
  vec_length = vector_v%num_elements
  call cel_omp_shared_work_distr(cel_omp_shared_work,&
                                 omp_start_idx,omp_end_idx,num_threads,&
                                 vec_length,vec_length,&
                                 cel_omp,output_on,err_code)

  tmp_real = 1.0_rk/alpha
  if(cel_omp%worker_num .le. num_threads) then
      vector_v%values(omp_start_idx:omp_end_idx,1_ik) = &
          vector_r%values(omp_start_idx:omp_end_idx,1_ik)*tmp_real
  end if

  if(cel_omp%worker_num .eq. 1) then
    cel_perf_counter%time(1) = cel_perf_counter%time(1) + (get_time() - time)
    cel_perf_counter%mult= cel_perf_counter%mult + vec_length
    cel_perf_counter%load= cel_perf_counter%load + 1_ik*vec_length*rk
    cel_perf_counter%store= cel_perf_counter%store + 1_ik*vec_length*rk
  end if
  
end subroutine cel_fop_vec_distr_scale_inv

!vector_r = boundary_vector - vector_r
!vector_d = jacobi*vector_r
subroutine cel_fop_vec_distr_calc_r_and_d_jacobi(vector_r,vector_d,boundary_vector,&
                                jacobi, cel_omp,cel_omp_shared_work,&
                                cel_perf_counter, output_on,err_code)
  type(cel_sp_mat_distr_vec_type), intent(inout)                    :: vector_r
  type(cel_sp_mat_distr_vec_type), intent(inout)                    :: vector_d
  type(cel_sp_mat_distr_vec_type), intent(in)                       :: boundary_vector
  real(kind=rk), dimension(:), allocatable, intent(in)              :: jacobi
  type(cel_omp_type), intent(in)                                    :: cel_omp
  type(cel_omp_shared_work_type), intent(inout)                     :: cel_omp_shared_work
  type(cel_perf_counter_type), intent(inout)                        :: cel_perf_counter
  logical(kind=lk)  , intent(in)                                    :: output_on
  integer(kind=ik), intent(out)                                     :: err_code
  integer(kind=ik) :: ii, vec_length
  real(kind=trk) :: time
  integer(kind=ik) :: omp_start_idx,omp_end_idx,num_threads
 
  err_code = 0_ik
  if(cel_omp%worker_num .eq. 1_ik) then
   time = get_time()
  end if
  vec_length = vector_r%num_elements
  call cel_omp_shared_work_distr(cel_omp_shared_work,&
                                 omp_start_idx,omp_end_idx,num_threads,&
                                 vec_length,vec_length,&
                                 cel_omp,output_on,err_code)

  if(cel_omp%worker_num .le. num_threads) then
      vector_r%values(omp_start_idx:omp_end_idx,1_ik) = &
         boundary_vector%values(omp_start_idx:omp_end_idx,1_ik)-vector_r%values(omp_start_idx:omp_end_idx,1_ik)
      vector_d%values(omp_start_idx:omp_end_idx,1_ik) = &
         jacobi(omp_start_idx:omp_end_idx)*vector_r%values(omp_start_idx:omp_end_idx,1_ik)
  end if

  !the last worker captures time(the last worker has usually more to do)
  if(cel_omp%worker_num .eq. 1_ik) then
    cel_perf_counter%time(1) = cel_perf_counter%time(1) + (get_time() - time)
    cel_perf_counter%add= cel_perf_counter%add + vec_length
    cel_perf_counter%mult= cel_perf_counter%add + vec_length
    cel_perf_counter%load= cel_perf_counter%load + 4_ik*vec_length*rk
    cel_perf_counter%store= cel_perf_counter%store + 2_ik*vec_length*rk
  end if
  
end subroutine cel_fop_vec_distr_calc_r_and_d_jacobi

!vector_r = boundary_vector - vector_r
subroutine cel_fop_vec_distr_calc_r(vector_r,boundary_vector,&
                                cel_omp,cel_omp_shared_work,&
                                cel_perf_counter, output_on,err_code)
  type(cel_sp_mat_distr_vec_type), intent(inout)                    :: vector_r
  type(cel_sp_mat_distr_vec_type), intent(in)                       :: boundary_vector
  type(cel_omp_type), intent(in)                                    :: cel_omp
  type(cel_omp_shared_work_type), intent(inout)                     :: cel_omp_shared_work
  type(cel_perf_counter_type), intent(inout)                        :: cel_perf_counter
  logical(kind=lk)  , intent(in)                                    :: output_on
  integer(kind=ik), intent(out)                                     :: err_code
  integer(kind=ik) :: ii, vec_length
  real(kind=trk) :: time
  integer(kind=ik) :: omp_start_idx,omp_end_idx,num_threads

  err_code = 0_ik
  if(cel_omp%worker_num .eq. 1_ik) then
   time = get_time()
  end if
  vec_length = vector_r%num_elements
  call cel_omp_shared_work_distr(cel_omp_shared_work,&
                                 omp_start_idx,omp_end_idx,num_threads,&
                                 vec_length,vec_length,&
                                 cel_omp,output_on,err_code)

  if(cel_omp%worker_num .le. num_threads) then
    do ii=omp_start_idx, omp_end_idx
      vector_r%values(ii,1_ik) = boundary_vector%values(ii,1_ik)-vector_r%values(ii,1_ik)
    end do
  end if
  
  !the last worker captures time(the last worker has usually more to do)
  if(cel_omp%worker_num .eq. 1_ik) then
    cel_perf_counter%time(1) = cel_perf_counter%time(1) + (get_time() - time)
    cel_perf_counter%add= cel_perf_counter%add + vec_length
    cel_perf_counter%load= cel_perf_counter%load + vec_length*rk
    cel_perf_counter%store= cel_perf_counter%store + vec_length*rk
  end if
  
end subroutine cel_fop_vec_distr_calc_r

!vector = 0.0
subroutine cel_fop_vec_distr_to_zero(vector,cel_omp,cel_omp_shared_work,&
                                cel_perf_counter, output_on,err_code)
  type(cel_sp_mat_distr_vec_type), intent(inout)                    :: vector
  type(cel_omp_type), intent(in)                                    :: cel_omp
  type(cel_omp_shared_work_type), intent(inout)                     :: cel_omp_shared_work
  type(cel_perf_counter_type), intent(inout)                        :: cel_perf_counter
  logical(kind=lk)  , intent(in)                                    :: output_on
  integer(kind=ik), intent(out)                                     :: err_code
  integer(kind=ik) :: vec_length
  real(kind=trk) :: time
  integer(kind=ik) :: omp_start_idx,omp_end_idx,num_threads
  
  err_code = 0_ik
  if(cel_omp%worker_num .eq. 1_ik) then
   time = get_time()
  end if
  vec_length = vector%num_elements
  call cel_omp_shared_work_distr(cel_omp_shared_work,&
                                 omp_start_idx,omp_end_idx,num_threads,&
                                 vec_length,vec_length,&
                                 cel_omp,output_on,err_code)

  if(cel_omp%worker_num .le. num_threads) then
    vector%values(omp_start_idx:omp_end_idx,1_ik) = 0.0_rk
  end if
  if(cel_omp%worker_num .eq. 1_ik) then
    cel_perf_counter%time(1) = cel_perf_counter%time(1) + (get_time() - time)
    cel_perf_counter%store= cel_perf_counter%store + vec_length*rk
  end if
  
end subroutine cel_fop_vec_distr_to_zero

!vector1 = 0.0; vector2 = 0.0
subroutine cel_fop_vecs_distr_to_zero(vector1,vector2,cel_omp,cel_omp_shared_work,&
                                cel_perf_counter, output_on,err_code)
  type(cel_sp_mat_distr_vec_type), intent(inout)                    :: vector1
  type(cel_sp_mat_distr_vec_type), intent(inout)                    :: vector2
  type(cel_omp_type), intent(in)                                    :: cel_omp
  type(cel_omp_shared_work_type), intent(inout)                     :: cel_omp_shared_work
  type(cel_perf_counter_type), intent(inout)                        :: cel_perf_counter
  logical(kind=lk)  , intent(in)                                    :: output_on
  integer(kind=ik), intent(out)                                     :: err_code
  integer(kind=ik) :: vec_length
  real(kind=trk) :: time
  integer(kind=ik) :: omp_start_idx,omp_end_idx,num_threads
  
  err_code = 0_ik
  if(cel_omp%worker_num .eq. 1_ik) then
   time = get_time()
  end if
  vec_length = vector1%num_elements
  call cel_omp_shared_work_distr(cel_omp_shared_work,&
                                 omp_start_idx,omp_end_idx,num_threads,&
                                 vec_length,vec_length,&
                                 cel_omp,output_on,err_code)

  if(cel_omp%worker_num .le. num_threads) then
    vector1%values(omp_start_idx:omp_end_idx,1_ik) = 0.0_rk
    vector2%values(omp_start_idx:omp_end_idx,1_ik) = 0.0_rk
  end if
  if(cel_omp%worker_num .eq. 1_ik) then
    cel_perf_counter%time(1) = cel_perf_counter%time(1) + (get_time() - time)
    cel_perf_counter%store= cel_perf_counter%store + 2_ik*vec_length*rk
  end if
  
end subroutine cel_fop_vecs_distr_to_zero
!vector_r = constant
subroutine cel_fop_vec_distr_to_constant(vector,constant,cel_omp,cel_omp_shared_work,&
                                cel_perf_counter, output_on,err_code)
  type(cel_sp_mat_distr_vec_type), intent(inout)                    :: vector
  real(kind=rk), intent(in)                                         :: constant
  type(cel_omp_type), intent(in)                                    :: cel_omp
  type(cel_omp_shared_work_type), intent(inout)                     :: cel_omp_shared_work
  type(cel_perf_counter_type), intent(inout)                        :: cel_perf_counter
  logical(kind=lk)  , intent(in)                                    :: output_on
  integer(kind=ik), intent(out)                                     :: err_code
  integer(kind=ik) :: vec_length
  real(kind=trk) :: time
  integer(kind=ik) :: omp_start_idx,omp_end_idx,num_threads
  
  err_code = 0_ik
  if(cel_omp%worker_num .eq. 1_ik) then
   time = get_time()
  end if
  vec_length = vector%num_elements
  call cel_omp_shared_work_distr(cel_omp_shared_work,&
                                 omp_start_idx,omp_end_idx,num_threads,&
                                 vec_length,vec_length,&
                                 cel_omp,output_on,err_code)

  if(cel_omp%worker_num .le. num_threads) then
     vector%values(omp_start_idx:omp_end_idx,1_ik) = constant
  end if
  if(cel_omp%worker_num .eq. 1_ik) then
    cel_perf_counter%time(1) = cel_perf_counter%time(1) + (get_time() - time)
    cel_perf_counter%store= cel_perf_counter%store + vec_length*rk
  end if
  
end subroutine cel_fop_vec_distr_to_constant

!vector_p(i) = vector_r(i)*vector_j(i)
subroutine cel_fop_vec_mult_distr(vector_p,vector_r,vector_j,cel_omp,&
                                cel_omp_shared_work,&
                                cel_perf_counter, output_on,err_code)
  type(cel_sp_mat_distr_vec_type), intent(inout)                    :: vector_p
  type(cel_sp_mat_distr_vec_type), intent(inout)                    :: vector_r
  type(cel_sp_mat_distr_vec_type), intent(inout)                    :: vector_j
  type(cel_omp_type), intent(in)                                    :: cel_omp
  type(cel_omp_shared_work_type), intent(inout)                     :: cel_omp_shared_work
  type(cel_perf_counter_type), intent(inout)                        :: cel_perf_counter
  logical(kind=lk)  , intent(in)                                    :: output_on
  integer(kind=ik), intent(out)                                     :: err_code
  integer(kind=ik) :: vec_length
  real(kind=trk) :: time
  integer(kind=ik) :: omp_start_idx,omp_end_idx,num_threads
  
  err_code = 0_ik
  if(cel_omp%worker_num .eq. 1_ik) then
   time = get_time()
  end if
  vec_length = vector_p%num_elements
  call cel_omp_shared_work_distr(cel_omp_shared_work,&
                                 omp_start_idx,omp_end_idx,num_threads,&
                                 vec_length,vec_length,&
                                 cel_omp,output_on,err_code)

  if(cel_omp%worker_num .le. num_threads) then
     vector_p%values(omp_start_idx:omp_end_idx,1_ik) = &
       vector_r%values(omp_start_idx:omp_end_idx,1_ik)*&
       vector_j%values(omp_start_idx:omp_end_idx,1_ik)
  end if
  if(cel_omp%worker_num .eq. 1_ik) then
    cel_perf_counter%time(1) = cel_perf_counter%time(1) + (get_time() - time)
    cel_perf_counter%store= cel_perf_counter%store + vec_length*rk
    cel_perf_counter%load= cel_perf_counter%store + 2_ik*vec_length*rk
    cel_perf_counter%mult= cel_perf_counter%mult + vec_length
  end if
  
end subroutine cel_fop_vec_mult_distr

!vector_p(i) = vector_r(i)*jacobi(i)
subroutine cel_fop_vec_mult_distr_jacobi(vector_p,vector_r,jacobi,cel_omp,&
                                cel_omp_shared_work,&
                                cel_perf_counter, output_on,err_code)
  type(cel_sp_mat_distr_vec_type), intent(inout)                    :: vector_p
  type(cel_sp_mat_distr_vec_type), intent(in)                       :: vector_r
  real(kind=rk), dimension(:), allocatable, intent(in)              :: jacobi
  type(cel_omp_type), intent(in)                                    :: cel_omp
  type(cel_omp_shared_work_type), intent(inout)                     :: cel_omp_shared_work
  type(cel_perf_counter_type), intent(inout)                        :: cel_perf_counter
  logical(kind=lk)  , intent(in)                                    :: output_on
  integer(kind=ik), intent(out)                                     :: err_code
  integer(kind=ik) :: vec_length
  real(kind=trk) :: time
  integer(kind=ik) :: omp_start_idx,omp_end_idx,num_threads
  
  err_code = 0_ik
  if(cel_omp%worker_num .eq. 1_ik) then
   time = get_time()
  end if
  vec_length = vector_p%num_elements
  call cel_omp_shared_work_distr(cel_omp_shared_work,&
                                 omp_start_idx,omp_end_idx,num_threads,&
                                 vec_length,vec_length,&
                                 cel_omp,output_on,err_code)

  if(cel_omp%worker_num .le. num_threads) then
     vector_p%values(omp_start_idx:omp_end_idx,1_ik) = &
       vector_r%values(omp_start_idx:omp_end_idx,1_ik)*&
       jacobi(omp_start_idx:omp_end_idx)
  end if
  if(cel_omp%worker_num .eq. 1_ik) then
    cel_perf_counter%time(1) = cel_perf_counter%time(1) + (get_time() - time)
    cel_perf_counter%store= cel_perf_counter%store + vec_length*rk
    cel_perf_counter%load= cel_perf_counter%store + 2_ik*vec_length*rk
    cel_perf_counter%mult= cel_perf_counter%mult + vec_length
  end if
  
end subroutine cel_fop_vec_mult_distr_jacobi

!norm2=||vec_1*vec_2||_2
subroutine cel_fop_vecs_distr_norm2_allreduce(norm2,distr_vector_1, &
                                distr_vector_2,&
                                cel_omp,cel_omp_shared_work,&
                                cel_perf_counter,output_on,err_code)
  real(kind=rk),intent(out)                                         :: norm2
  type(cel_sp_mat_distr_vec_type), intent(in)                       :: distr_vector_1
  type(cel_sp_mat_distr_vec_type), intent(in)                       :: distr_vector_2
  type(cel_omp_type), intent(in)                                    :: cel_omp
  type(cel_omp_shared_work_type), intent(inout)                     :: cel_omp_shared_work
  type(cel_perf_counter_type), intent(inout)                        :: cel_perf_counter
  logical(kind=lk)  , intent(in)                                    :: output_on
  integer(kind=ik), intent(out)                                     :: err_code
  integer(kind=ik)   :: vec_length, ii
  real(kind=rk)     :: time, tmp_scalar, norm2_local
  integer(kind=ik) :: omp_start_idx,omp_end_idx,num_threads
  integer(kind=mpik) :: ierr
    
  err_code = 0_ik
  if(cel_omp%worker_num .eq. 1_ik) then
   time = get_time()
   norm2 = 0.0_rk
  end if
  !$omp barrier
  vec_length = distr_vector_1%num_elements
  call cel_omp_shared_work_distr(cel_omp_shared_work,&
                                 omp_start_idx,omp_end_idx,num_threads,&
                                 vec_length,vec_length,&
                                 cel_omp,output_on,err_code)

  if(cel_omp%worker_num .le. num_threads) then
    norm2_local = 0.0_rk
    do ii=omp_start_idx, omp_end_idx
      norm2_local = norm2_local + distr_vector_1%values(ii,1_ik)*&
                    distr_vector_2%values(ii,1_ik)
    end do
    !$omp atomic update
    norm2 = norm2 + norm2_local
  end if
  !$omp barrier
  if(cel_omp%is_master) then
    tmp_scalar = 0.0_rk
    call MPI_ALLREDUCE(norm2,tmp_scalar, 1_mpik, cel_omp%mpi_type_rk,&
        MPI_SUM, cel_omp%master_comm, ierr)
    norm2 = SQRT(tmp_scalar)
  end if
  !$omp barrier
  if(cel_omp%worker_num .eq. 1_ik) then
    cel_perf_counter%time(1) = cel_perf_counter%time(1) + (get_time() - time)
    cel_perf_counter%mult= cel_perf_counter%mult + vec_length
    cel_perf_counter%add= cel_perf_counter%add + vec_length
    cel_perf_counter%load= cel_perf_counter%load + 1_ik*vec_length*rk
    cel_perf_counter%num_allreduces = cel_perf_counter%num_allreduces + 1_ik
  end if
end subroutine cel_fop_vecs_distr_norm2_allreduce
!norm2=||vec_1*vec_1||_2
subroutine cel_fop_vec_distr_norm2_allreduce(norm2,distr_vector, &
                                cel_omp,cel_omp_shared_work,&
                                cel_perf_counter, output_on,err_code)
  real(kind=rk),intent(out)                                         :: norm2
  type(cel_sp_mat_distr_vec_type), intent(in)                       :: distr_vector
  type(cel_omp_type), intent(in)                                    :: cel_omp
  type(cel_omp_shared_work_type), intent(inout)                     :: cel_omp_shared_work
  type(cel_perf_counter_type), intent(inout)                        :: cel_perf_counter
  logical(kind=lk)  , intent(in)                                    :: output_on
  integer(kind=ik), intent(out)                                     :: err_code
  integer(kind=ik)   :: vec_length, ii
  real(kind=rk)     :: time, tmp_scalar, norm2_local
  integer(kind=ik) :: omp_start_idx,omp_end_idx,num_threads
  integer(kind=mpik) :: ierr
    
  err_code = 0_ik
  if(cel_omp%worker_num .eq. 1_ik) then
   time = get_time()
   norm2 = 0.0_rk
  end if
  !$omp barrier
  vec_length = distr_vector%num_elements
  call cel_omp_shared_work_distr(cel_omp_shared_work,&
                                 omp_start_idx,omp_end_idx,num_threads,&
                                 vec_length,vec_length,&
                                 cel_omp,output_on,err_code)

  if(cel_omp%worker_num .le. num_threads) then
    norm2_local = 0.0_rk
    do ii=omp_start_idx, omp_end_idx
      norm2_local = norm2_local + distr_vector%values(ii,1_ik)*distr_vector%values(ii,1_ik)
    end do
    !$omp atomic update
    norm2 = norm2 + norm2_local
  end if
  !$omp barrier
  if(cel_omp%is_master) then
    tmp_scalar = 0.0_rk
    call MPI_ALLREDUCE(norm2,tmp_scalar, 1_mpik, cel_omp%mpi_type_rk,&
      MPI_SUM, cel_omp%master_comm, ierr)
    norm2 = SQRT(tmp_scalar)
  end if
  !$omp barrier
  if(cel_omp%worker_num .eq. 1_ik) then
    cel_perf_counter%time(1) = cel_perf_counter%time(1) + (get_time() - time)
    cel_perf_counter%mult= cel_perf_counter%mult + vec_length
    cel_perf_counter%add= cel_perf_counter%add + vec_length
    cel_perf_counter%load= cel_perf_counter%load + 1_ik*vec_length*rk
    cel_perf_counter%num_allreduces = cel_perf_counter%num_allreduces + 1_ik
  end if
end subroutine cel_fop_vec_distr_norm2_allreduce

!Calculate norm of the distr_vector and check wether to start gmres 
!norm2=||vec_1*vec_1||_2
subroutine cel_fop_vec_distr_norm2_start_gmres_allreduce(norm2,&
                                gmres_gamma, gmres_conv_reason,&
                                conv_epsilon,distr_vector, &
                                cel_omp,cel_omp_shared_work,&
                                cel_perf_counter, output_on,err_code)
  real(kind=rk),intent(out)                                         :: norm2
  real(kind=rk),intent(inout)                                       :: gmres_gamma
  integer(kind=ik),intent(out)                                      :: gmres_conv_reason
  real(kind=rk),intent(in)                                          :: conv_epsilon
  type(cel_sp_mat_distr_vec_type), intent(in)                       :: distr_vector
  type(cel_omp_type), intent(in)                                    :: cel_omp
  type(cel_omp_shared_work_type), intent(inout)                     :: cel_omp_shared_work
  type(cel_perf_counter_type), intent(inout)                        :: cel_perf_counter
  logical(kind=lk)  , intent(in)                                    :: output_on
  integer(kind=ik), intent(out)                                     :: err_code
  integer(kind=ik)   :: vec_length, ii
  real(kind=rk)     :: time, tmp_scalar, norm2_local
  integer(kind=ik) :: omp_start_idx,omp_end_idx,num_threads
  integer(kind=mpik) :: ierr
    
  err_code = 0_ik
  if(cel_omp%worker_num .eq. 1_ik) then
   time = get_time()
   norm2 = 0.0_rk
  end if
  !$omp barrier
  vec_length = distr_vector%num_elements
  call cel_omp_shared_work_distr(cel_omp_shared_work,&
                                 omp_start_idx,omp_end_idx,num_threads,&
                                 vec_length,vec_length,&
                                 cel_omp,output_on,err_code)

  if(cel_omp%worker_num .le. num_threads) then
    norm2_local = 0.0_rk
    do ii=omp_start_idx, omp_end_idx
      norm2_local = norm2_local + distr_vector%values(ii,1_ik)*distr_vector%values(ii,1_ik)
    end do
    !$omp atomic update
    norm2 = norm2 + norm2_local
  end if
  !$omp barrier
  if(cel_omp%is_master) then
    tmp_scalar = 0.0_rk
    call MPI_ALLREDUCE(norm2,tmp_scalar, 1_mpik, cel_omp%mpi_type_rk,&
      MPI_SUM, cel_omp%master_comm, ierr)
    norm2 = SQRT(tmp_scalar)
    gmres_gamma = norm2
    if(abs(norm2) .le. conv_epsilon ) then
      gmres_conv_reason = GMRES_REASON_X0
    else
      gmres_conv_reason = GMRES_REASON_UNDEX
    end if
  end if
  !$omp barrier
  if(cel_omp%worker_num .eq. 1_ik) then
    cel_perf_counter%time(1) = cel_perf_counter%time(1) + (get_time() - time)
    cel_perf_counter%mult= cel_perf_counter%mult + vec_length
    cel_perf_counter%add= cel_perf_counter%add + vec_length
    cel_perf_counter%load= cel_perf_counter%load + 1_ik*vec_length*rk
    cel_perf_counter%num_allreduces = cel_perf_counter%num_allreduces + 1_ik
  end if
end subroutine cel_fop_vec_distr_norm2_start_gmres_allreduce


!sigma_new=distr_vector^T * distr_vector and sigma_old = sigma_new and alpha = 1.0
subroutine cel_fop_vec_distr_sigma_new_1_allreduce(sigma_new,sigma_old,alpha,distr_vector, &
                                cel_omp,cel_omp_shared_work,&
                                cel_perf_counter, output_on,err_code)
  real(kind=rk),intent(inout)                                       :: sigma_new
  real(kind=rk),intent(inout)                                       :: sigma_old
  real(kind=rk),intent(inout)                                       :: alpha
  type(cel_sp_mat_distr_vec_type), intent(in)                       :: distr_vector
  type(cel_omp_type), intent(in)                                    :: cel_omp
  type(cel_omp_shared_work_type), intent(inout)                     :: cel_omp_shared_work
  type(cel_perf_counter_type), intent(inout)                        :: cel_perf_counter
  logical(kind=lk)  , intent(in)                                    :: output_on
  integer(kind=ik), intent(out)                                     :: err_code
  integer(kind=ik)   :: vec_length, ii
  real(kind=rk)     :: time, tmp_scalar, sigma_local
  integer(kind=ik) :: omp_start_idx,omp_end_idx,num_threads
  integer(kind=mpik) :: ierr
  
  err_code = 0_ik
  if(cel_omp%worker_num .eq. 1_ik) then
   time = get_time()
   sigma_new = 0.0_rk
  end if
  !$omp barrier
  vec_length = distr_vector%num_elements
  call cel_omp_shared_work_distr(cel_omp_shared_work,&
                                 omp_start_idx,omp_end_idx,num_threads,&
                                 vec_length,vec_length,&
                                 cel_omp,output_on,err_code)

  if(cel_omp%worker_num .le. num_threads) then
    sigma_local = 0.0_rk
    do ii=omp_start_idx, omp_end_idx
      sigma_local = sigma_local + distr_vector%values(ii,1_ik)*distr_vector%values(ii,1_ik)
    end do
    !$omp atomic update
    sigma_new = sigma_new + sigma_local
  end if
  !$omp barrier
  if(cel_omp%is_master) then
    tmp_scalar = 0.0_rk
    call MPI_ALLREDUCE(sigma_new,tmp_scalar, 1_mpik, cel_omp%mpi_type_rk,&
      MPI_SUM, cel_omp%master_comm, ierr)
    sigma_new = tmp_scalar
    sigma_old = tmp_scalar
    alpha = 1.0_rk
  end if
  !$omp barrier
  if(cel_omp%worker_num .eq. 1_ik) then
    cel_perf_counter%time(1) = cel_perf_counter%time(1) + (get_time() - time)
    cel_perf_counter%mult= cel_perf_counter%mult + vec_length
    cel_perf_counter%add= cel_perf_counter%add + vec_length
    cel_perf_counter%load= cel_perf_counter%load + 1_ik*vec_length*rk
    cel_perf_counter%num_allreduces = cel_perf_counter%num_allreduces + 1_ik
  end if
end subroutine cel_fop_vec_distr_sigma_new_1_allreduce

!sigma_old = sigma_new and sigma_new=distr_vector^T * distr_vector 
subroutine cel_fop_vec_distr_sigma_new_2_allreduce(sigma_new,sigma_old, distr_vector, &
                                cel_omp,cel_omp_shared_work,&
                                cel_perf_counter, output_on,err_code)
  real(kind=rk),intent(inout)                                       :: sigma_new
  real(kind=rk),intent(inout)                                       :: sigma_old
  type(cel_sp_mat_distr_vec_type), intent(in)                       :: distr_vector
  type(cel_omp_type), intent(in)                                    :: cel_omp
  type(cel_omp_shared_work_type), intent(inout)                     :: cel_omp_shared_work
  type(cel_perf_counter_type), intent(inout)                        :: cel_perf_counter
  logical(kind=lk)  , intent(in)                                    :: output_on
  integer(kind=ik), intent(out)                                     :: err_code
  integer(kind=ik)   :: vec_length, ii
  real(kind=rk)     :: time, tmp_scalar, sigma_local
  integer(kind=ik) :: omp_start_idx,omp_end_idx,num_threads
  integer(kind=mpik) :: ierr
  
  err_code = 0_ik
  if(cel_omp%worker_num .eq. 1_ik) then
   time = get_time()
   sigma_old = sigma_new
   sigma_new = 0.0_rk
  end if
  !$omp barrier
  vec_length = distr_vector%num_elements
  call cel_omp_shared_work_distr(cel_omp_shared_work,&
                                 omp_start_idx,omp_end_idx,num_threads,&
                                 vec_length,vec_length,&
                                 cel_omp,output_on,err_code)

  if(cel_omp%worker_num .le. num_threads) then
    sigma_local = 0.0_rk
    do ii=omp_start_idx, omp_end_idx
      sigma_local = sigma_local + distr_vector%values(ii,1_ik)*distr_vector%values(ii,1_ik)
    end do
    !$omp atomic update
    sigma_new = sigma_new + sigma_local
  end if
  !$omp barrier
  if(cel_omp%is_master) then
    tmp_scalar = 0.0_rk
    call MPI_ALLREDUCE(sigma_new,tmp_scalar, 1_mpik, cel_omp%mpi_type_rk,&
      MPI_SUM, cel_omp%master_comm, ierr)
    sigma_new = tmp_scalar
  end if
  !$omp barrier
  if(cel_omp%worker_num .eq. 1_ik) then
    cel_perf_counter%time(1) = cel_perf_counter%time(1) + (get_time() - time)
    cel_perf_counter%mult= cel_perf_counter%mult + vec_length
    cel_perf_counter%add= cel_perf_counter%add + vec_length
    cel_perf_counter%load= cel_perf_counter%load + 1_ik*vec_length*rk
    cel_perf_counter%num_allreduces = cel_perf_counter%num_allreduces + 1_ik
  end if
end subroutine cel_fop_vec_distr_sigma_new_2_allreduce


!alpha=sigma_new/distr_vector_d^T*distr_vector_q
subroutine cel_fop_vec_distr_alpha_allreduce(alpha,sigma_new,distr_vector_d, &
                                distr_vector_q,cel_omp,cel_omp_shared_work,&
                                cel_perf_counter, output_on,err_code)
  real(kind=rk),intent(out)                                         :: alpha
  real(kind=rk),intent(in)                                          :: sigma_new
  type(cel_sp_mat_distr_vec_type), intent(in)                       :: distr_vector_d
  type(cel_sp_mat_distr_vec_type), intent(in)                       :: distr_vector_q
  type(cel_omp_type), intent(in)                                    :: cel_omp
  type(cel_omp_shared_work_type), intent(inout)                     :: cel_omp_shared_work
  type(cel_perf_counter_type), intent(inout)                        :: cel_perf_counter
  logical(kind=lk)  , intent(in)                                    :: output_on
  integer(kind=ik), intent(out)                                     :: err_code
  integer(kind=ik)   :: vec_length, ii
  real(kind=rk)     :: time, tmp_scalar, local_scalar
  integer(kind=ik) :: omp_start_idx,omp_end_idx,num_threads
  integer(kind=mpik) :: ierr

  err_code = 0_ik
  if(cel_omp%worker_num .eq. 1_ik) then
   time = get_time()
   alpha = 0.0_rk
  end if
  !$omp barrier
  vec_length = distr_vector_d%num_elements
  call cel_omp_shared_work_distr(cel_omp_shared_work,&
                                 omp_start_idx,omp_end_idx,num_threads,&
                                 vec_length,vec_length,&
                                 cel_omp,output_on,err_code)

  if(cel_omp%worker_num .le. num_threads) then
    local_scalar = 0.0_rk
    do ii=omp_start_idx, omp_end_idx
      local_scalar = local_scalar + distr_vector_d%values(ii,1_ik)*distr_vector_q%values(ii,1_ik)
    end do
    !$omp atomic update
    alpha = alpha + local_scalar
  end if
  !$omp barrier
  if(cel_omp%is_master) then
    if(abs(alpha) .le. cel_real_small) then
      call cel_error("cel_sp_mat_new cel_fop_vec_distr_alpha_allreduce div. by zero (alpha=sigma_new/d^T*q)",&
        err_code, cel_is_in_debug, .TRUE._lk)
    end if
    tmp_scalar = 0.0_rk
    call MPI_ALLREDUCE(alpha,tmp_scalar, 1_mpik, cel_omp%mpi_type_rk,&
      MPI_SUM, cel_omp%master_comm, ierr)
    alpha = sigma_new/tmp_scalar
  end if
  !$omp barrier
  if(cel_omp%worker_num .eq. 1_ik) then
    cel_perf_counter%time(1) = cel_perf_counter%time(1) + (get_time() - time)
    cel_perf_counter%mult= cel_perf_counter%mult + vec_length
    cel_perf_counter%add= cel_perf_counter%add + vec_length
    cel_perf_counter%load= cel_perf_counter%load + 2_ik*vec_length*rk
    cel_perf_counter%num_allreduces = cel_perf_counter%num_allreduces + 1_ik
  end if
end subroutine cel_fop_vec_distr_alpha_allreduce

!alpha=distr_vector_r^T*distr_vector_d
subroutine cel_fop_vec_distr_dot_allreduce(alpha,distr_vector_r, &
                                distr_vector_d,cel_omp,cel_omp_shared_work,&
                                cel_perf_counter, output_on,err_code)
  real(kind=rk),intent(out)                                         :: alpha
  type(cel_sp_mat_distr_vec_type), intent(in)                       :: distr_vector_r
  type(cel_sp_mat_distr_vec_type), intent(in)                       :: distr_vector_d
  type(cel_omp_type), intent(in)                                    :: cel_omp
  type(cel_omp_shared_work_type), intent(inout)                     :: cel_omp_shared_work
  type(cel_perf_counter_type), intent(inout)                        :: cel_perf_counter
  logical(kind=lk)  , intent(in)                                    :: output_on
  integer(kind=ik), intent(out)                                     :: err_code
  integer(kind=ik)   :: vec_length, ii
  real(kind=rk)     :: time, tmp_scalar, local_scalar
  integer(kind=ik) :: omp_start_idx,omp_end_idx,num_threads
  integer(kind=mpik) :: ierr

  err_code = 0_ik
  if(cel_omp%worker_num .eq. 1_ik) then
   time = get_time()
   alpha = 0.0_rk
  end if
  !$omp barrier
  vec_length = distr_vector_r%num_elements
  call cel_omp_shared_work_distr(cel_omp_shared_work,&
                                 omp_start_idx,omp_end_idx,num_threads,&
                                 vec_length,vec_length,&
                                 cel_omp,output_on,err_code)

  if(cel_omp%worker_num .le. num_threads) then
    local_scalar = 0.0_rk
    do ii=omp_start_idx, omp_end_idx
      local_scalar = local_scalar + distr_vector_r%values(ii,1_ik)*distr_vector_d%values(ii,1_ik)
    end do
    !$omp atomic update
    alpha = alpha + local_scalar
  end if
  !$omp barrier
  if(cel_omp%is_master) then
    tmp_scalar = 0.0_rk
    call MPI_ALLREDUCE(alpha,tmp_scalar, 1_mpik, cel_omp%mpi_type_rk,&
      MPI_SUM, cel_omp%master_comm, ierr)
    alpha = tmp_scalar
  end if
  !$omp barrier
  if(cel_omp%worker_num .eq. 1_ik) then
    cel_perf_counter%time(1) = cel_perf_counter%time(1) + (get_time() - time)
    cel_perf_counter%mult= cel_perf_counter%mult + vec_length
    cel_perf_counter%add= cel_perf_counter%add + vec_length
    cel_perf_counter%load= cel_perf_counter%load + 2_ik*vec_length*rk
    cel_perf_counter%num_allreduces = cel_perf_counter%num_allreduces + 1_ik
  end if
end subroutine cel_fop_vec_distr_dot_allreduce


!sigma=(r*z)_2/alpha; alpha=(r*z)_2
subroutine cel_fop_vec_distr_alpha_sigma_r_z_alpha_allreduce(alpha,sigma,distr_vector_r, &
                                distr_vector_z,cel_omp,cel_omp_shared_work,&
                                cel_perf_counter, output_on,err_code)
  real(kind=rk),intent(inout)                                       :: alpha
  real(kind=rk),intent(out)                                         :: sigma
  type(cel_sp_mat_distr_vec_type), intent(in)                       :: distr_vector_r
  type(cel_sp_mat_distr_vec_type), intent(in)                       :: distr_vector_z
  type(cel_omp_type), intent(in)                                    :: cel_omp
  type(cel_omp_shared_work_type), intent(inout)                     :: cel_omp_shared_work
  type(cel_perf_counter_type), intent(inout)                        :: cel_perf_counter
  logical(kind=lk)  , intent(in)                                    :: output_on
  integer(kind=ik), intent(out)                                     :: err_code
  integer(kind=ik)   :: vec_length, ii
  real(kind=rk)     :: time, tmp_scalar, local_scalar
  integer(kind=ik) :: omp_start_idx,omp_end_idx,num_threads
  integer(kind=mpik) :: ierr

  err_code = 0_ik
  if(cel_omp%worker_num .eq. 1_ik) then
   time = get_time()
   sigma = 0.0_rk
  end if
  !$omp barrier
  vec_length = distr_vector_r%num_elements
  call cel_omp_shared_work_distr(cel_omp_shared_work,&
                                 omp_start_idx,omp_end_idx,num_threads,&
                                 vec_length,vec_length,&
                                 cel_omp,output_on,err_code)

  if(cel_omp%worker_num .le. num_threads) then
    local_scalar = 0.0_rk
    do ii=omp_start_idx, omp_end_idx
      local_scalar = local_scalar + distr_vector_r%values(ii,1_ik)*distr_vector_z%values(ii,1_ik)
    end do
    !$omp atomic update
    sigma = sigma + local_scalar
  end if
  !$omp barrier
  if(cel_omp%is_master) then
    tmp_scalar = 0.0_rk
    call MPI_ALLREDUCE(sigma,tmp_scalar, 1_mpik, cel_omp%mpi_type_rk,&
      MPI_SUM, cel_omp%master_comm, ierr)
    sigma = tmp_scalar/alpha
    alpha = tmp_scalar
  end if
  !$omp barrier
  if(cel_omp%worker_num .eq. 1_ik) then
    cel_perf_counter%time(1) = cel_perf_counter%time(1) + (get_time() - time)
    cel_perf_counter%mult= cel_perf_counter%mult + vec_length
    cel_perf_counter%add= cel_perf_counter%add + vec_length
    cel_perf_counter%load= cel_perf_counter%load + 2_ik*vec_length*rk
    cel_perf_counter%num_allreduces = cel_perf_counter%num_allreduces + 1_ik
  end if
end subroutine cel_fop_vec_distr_alpha_sigma_r_z_alpha_allreduce


!x_(i)<-x_(i)+alpha*d_(i)
subroutine cel_fop_vec_distr_scale(distr_vector_x,alpha, distr_vector_d,&
                                cel_omp,cel_omp_shared_work,&
                                cel_perf_counter, output_on,err_code)
  type(cel_sp_mat_distr_vec_type), intent(inout)                    :: distr_vector_x
  real(kind=rk),intent(in)                                          :: alpha
  type(cel_sp_mat_distr_vec_type), intent(in)                       :: distr_vector_d
  type(cel_omp_type), intent(in)                                    :: cel_omp
  type(cel_omp_shared_work_type), intent(inout)                     :: cel_omp_shared_work
  type(cel_perf_counter_type), intent(inout)                        :: cel_perf_counter
  logical(kind=lk)  , intent(in)                                    :: output_on
  integer(kind=ik), intent(out)                                     :: err_code
  integer(kind=ik)   :: vec_length
  integer(kind=ik) :: omp_start_idx,omp_end_idx,num_threads
  real(kind=rk)     :: time

  err_code = 0_ik
  if(cel_omp%worker_num .eq. 1_ik) then
   time = get_time()
  end if
  vec_length = distr_vector_x%num_elements
  call cel_omp_shared_work_distr(cel_omp_shared_work,&
                                 omp_start_idx,omp_end_idx,num_threads,&
                                 vec_length,vec_length,&
                                 cel_omp,output_on,err_code)

  if(cel_omp%worker_num .le. num_threads) then
      distr_vector_x%values(omp_start_idx:omp_end_idx,1_ik) = &
        distr_vector_x%values(omp_start_idx:omp_end_idx,1_ik) + &
        alpha*distr_vector_d%values(omp_start_idx:omp_end_idx,1_ik)
  end if
  !$omp barrier
  if(cel_omp%worker_num .eq. 1_ik) then
    cel_perf_counter%time(1) = cel_perf_counter%time(1) + (get_time() - time)
    cel_perf_counter%add= cel_perf_counter%add + vec_length
    cel_perf_counter%mult= cel_perf_counter%mult + vec_length
    cel_perf_counter%load= cel_perf_counter%load + 2_ik*vec_length*rk
    cel_perf_counter%store= cel_perf_counter%store + vec_length*rk
  end if
end subroutine cel_fop_vec_distr_scale

!x_(i)<-x_(i)-alpha*d_(i)
subroutine cel_fop_vec_distr_scale_minus(distr_vector_x,alpha, distr_vector_d,&
                                cel_omp,cel_omp_shared_work,&
                                cel_perf_counter, output_on,err_code)
  type(cel_sp_mat_distr_vec_type), intent(inout)                    :: distr_vector_x
  real(kind=rk),intent(in)                                          :: alpha
  type(cel_sp_mat_distr_vec_type), intent(in)                       :: distr_vector_d
  type(cel_omp_type), intent(in)                                    :: cel_omp
  type(cel_omp_shared_work_type), intent(inout)                     :: cel_omp_shared_work
  type(cel_perf_counter_type), intent(inout)                        :: cel_perf_counter
  logical(kind=lk)  , intent(in)                                    :: output_on
  integer(kind=ik), intent(out)                                     :: err_code
  integer(kind=ik)   :: vec_length
  integer(kind=ik) :: omp_start_idx,omp_end_idx,num_threads
  real(kind=rk)     :: time

  err_code = 0_ik
  if(cel_omp%worker_num .eq. 1_ik) then
   time = get_time()
  end if
  vec_length = distr_vector_x%num_elements
  call cel_omp_shared_work_distr(cel_omp_shared_work,&
                                 omp_start_idx,omp_end_idx,num_threads,&
                                 vec_length,vec_length,&
                                 cel_omp,output_on,err_code)

  if(cel_omp%worker_num .le. num_threads) then
      distr_vector_x%values(omp_start_idx:omp_end_idx,1_ik) = &
        distr_vector_x%values(omp_start_idx:omp_end_idx,1_ik) - &
        alpha*distr_vector_d%values(omp_start_idx:omp_end_idx,1_ik)
  end if
  !$omp barrier
  if(cel_omp%worker_num .eq. 1_ik) then
    cel_perf_counter%time(1) = cel_perf_counter%time(1) + (get_time() - time)
    cel_perf_counter%add= cel_perf_counter%add + vec_length
    cel_perf_counter%mult= cel_perf_counter%mult + vec_length
    cel_perf_counter%load= cel_perf_counter%load + 2_ik*vec_length*rk
    cel_perf_counter%store= cel_perf_counter%store + vec_length*rk
  end if
end subroutine cel_fop_vec_distr_scale_minus

!d_(i)=r_(i)+betta*d_(i) VecAYPX
subroutine cel_fop_vec_distr_scale_1(distr_vector_d,betta, distr_vector_r,&
                                cel_omp,cel_omp_shared_work,&
                                cel_perf_counter, output_on,err_code)
  type(cel_sp_mat_distr_vec_type), intent(inout)                    :: distr_vector_d
  real(kind=rk),intent(in)                                          :: betta
  type(cel_sp_mat_distr_vec_type), intent(in)                       :: distr_vector_r
  type(cel_omp_type), intent(in)                                    :: cel_omp
  type(cel_omp_shared_work_type), intent(inout)                     :: cel_omp_shared_work
  type(cel_perf_counter_type), intent(inout)                        :: cel_perf_counter
  logical(kind=lk)  , intent(in)                                    :: output_on
  integer(kind=ik), intent(out)                                     :: err_code
  integer(kind=ik)   :: vec_length
  real(kind=rk)     :: time
  integer(kind=ik) :: omp_start_idx,omp_end_idx,num_threads
  
  err_code = 0_ik
  if(cel_omp%worker_num .eq. 1_ik) then
   time = get_time()
  end if
  vec_length = distr_vector_d%num_elements
  call cel_omp_shared_work_distr(cel_omp_shared_work,&
                                 omp_start_idx,omp_end_idx,num_threads,&
                                 vec_length,vec_length,&
                                 cel_omp,output_on,err_code)

  if(cel_omp%worker_num .le. num_threads) then
      distr_vector_d%values(omp_start_idx:omp_end_idx,1_ik) = &
        betta*distr_vector_d%values(omp_start_idx:omp_end_idx,1_ik) + &
        distr_vector_r%values(omp_start_idx:omp_end_idx,1_ik)
  end if
  !$omp barrier
  if(cel_omp%worker_num .eq. 1_ik) then
    cel_perf_counter%time(1) = cel_perf_counter%time(1) + (get_time() - time)
    cel_perf_counter%add= cel_perf_counter%add + vec_length
    cel_perf_counter%mult= cel_perf_counter%mult + vec_length
    cel_perf_counter%load= cel_perf_counter%load + 2_ik*vec_length*rk
    cel_perf_counter%store= cel_perf_counter%store + vec_length*rk
  end if
end subroutine cel_fop_vec_distr_scale_1



end module cel_fop_vec_distr_module

