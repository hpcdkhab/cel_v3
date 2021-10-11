! Definition of the diagonal scaling of the matrix
! File:   cel_sp_mat_distr_scale_module.f90
! Project CRESTA Exascale library (see details on https://www.cresta-project.eu)
! Send you email about the bugs to
! The High Performance Computing Center Stuttgart (HLRS) of the University of Stuttgart
! Dmitry Khabi,  (email: khabi@hlrs.de)
! Created on Jun 14, 2013
!

module cel_sp_mat_distr_format_module
use cel_types_module
use cel_timer_interface_module
use cel_sp_mat_module
use cel_sp_mat_distr_module
use cel_perf_module
use cel_base_print_module
use cel_omp_module
use cel_omp_shared_work_module
use cel_sp_mat_distr_vec_module
use cel_sp_mat_jad_converter_module
use MPI
use OMP_LIB
implicit none
contains

!Format matrix in various formats
!coo and csr must be already done
subroutine cel_sp_mat_distr_format(a_sp_mats,vtxdist,cel_omp_shared_work,&
                               cel_omp,cel_perf_counter,output_on,err_code)
  type(cel_sp_mat_type), dimension(:), allocatable,intent(inout)    :: a_sp_mats
  integer(kind=ik), dimension(:), allocatable, intent(in)           :: vtxdist 
  type(cel_omp_shared_work_type), intent(inout)                     :: cel_omp_shared_work
  type(cel_omp_type), intent(in)                                    :: cel_omp
  type(cel_perf_counter_type), intent(inout)                        :: cel_perf_counter
  logical(kind=lk)  , intent(in)                                    :: output_on
  integer(kind=ik), intent(out)                                     :: err_code
  integer(kind=ik) :: ii
 
  err_code = 0_ik

  do ii=1_ik, size(a_sp_mats)
      call cel_sp_mat_distr_format_type(a_sp_mats(ii),vtxdist,cel_omp,cel_omp_shared_work,&
                               cel_perf_counter,output_on,err_code)
  end do
  
end subroutine cel_sp_mat_distr_format

!format  local sub matrices if needed
subroutine cel_sp_mat_distr_format_type(a_sp_mat,vtxdist,&
                            cel_omp,cel_omp_shared_work,cel_perf_counter, output_on,err_code)
  type(cel_sp_mat_type) ,intent(inout)                              :: a_sp_mat
  integer(kind=ik), dimension(:), allocatable, intent(in)           :: vtxdist 
  type(cel_omp_type), intent(in)                                    :: cel_omp
  type(cel_omp_shared_work_type), intent(inout)                     :: cel_omp_shared_work
  type(cel_perf_counter_type), intent(inout)                        :: cel_perf_counter
  logical(kind=lk)  , intent(in)                                    :: output_on
  integer(kind=ik), intent(out)                                     :: err_code

  err_code=0_ik
  
  select case (a_sp_mat%mv_algorithm) 
     case(MV_CSR)
        if(.not. a_sp_mat%is_in_csr) then
               call cel_error("Error cel_sp_mat_distr_format_type: csr matrix format is empty", &
                   err_code,cel_is_in_debug,.TRUE._lk)
        end if
     case(MV_JAD)
        if(.not. a_sp_mat%is_in_jad) then
          if(cel_omp%is_master) then
            call cel_sp_mat_jad_from_sp_mat_local(a_sp_mat,vtxdist,&
                                         cel_perf_counter,err_code)
            call cel_error("Error cel_sp_mat_distr_format_type: couldn't convert the matrix in jad format", &
                   err_code,cel_is_in_debug,.TRUE._lk)
            a_sp_mat%is_in_jad = .true._lk
          end if
          !$omp barrier
        end if
     case(MV_COO)
        if(.not. a_sp_mat%is_in_csr) then
          call cel_error("Error cel_sp_mat_distr_format_type: coo matrix format is empty", &
                   err_code,cel_is_in_debug,.TRUE._lk)
        end if
     case default
        call cel_error("Error cel_sp_mat_distr_format_type: no format defined", &
                   err_code,cel_is_in_debug,.TRUE._lk)
      end select
 
end subroutine cel_sp_mat_distr_format_type

!scale local sub matrices
subroutine cel_sp_mat_distr_scale(a_sp_mats,&
                            cel_omp,cel_omp_shared_work,cel_perf_counter, output_on,err_code)
  type(cel_sp_mat_type), dimension(:), allocatable,intent(inout)    :: a_sp_mats
  type(cel_omp_type), intent(inout)                                 :: cel_omp
  type(cel_omp_shared_work_type), intent(inout)                     :: cel_omp_shared_work
  type(cel_perf_counter_type), intent(inout)                        :: cel_perf_counter
  logical(kind=lk)  , intent(in)                                    :: output_on
  integer(kind=ik), intent(out)                                     :: err_code
  integer(kind=ik) :: start_idx,end_idx,omp_start_idx,omp_end_idx,num_threads
  integer(kind=ik) :: ii, length,sys_stat
  err_code = 0_ik

  length = size(a_sp_mats(1_ik)%diagonal)
  call cel_omp_shared_work_start(cel_omp_shared_work,cel_omp%num_threads)
  if(cel_omp%is_master) then
    if(.not. allocated(a_sp_mats(1_ik)%scaled_diagonal)) then
      allocate(a_sp_mats(1_ik)%scaled_diagonal(length), stat=sys_stat)  
      call cel_error("Error cel_sp_mat_distr_scale alloc 1", sys_stat,&
                        cel_is_in_debug,.TRUE._lk)
    end if
  end if
  call cel_omp_shared_work_end_and_wait(cel_omp_shared_work,cel_omp,LOOP_SLEEP_DEFLENGTH)
  !$omp barrier
  call cel_omp_shared_work_start(cel_omp_shared_work,cel_omp%num_threads)
  call cel_omp_shared_work_distr(cel_omp_shared_work,&
                                 omp_start_idx,omp_end_idx,num_threads,&
                                 length,length,&
                                 cel_omp,output_on,err_code)
  if(cel_omp%worker_num .le. num_threads) then
    a_sp_mats(1_ik)%scaled_diagonal(omp_start_idx:omp_end_idx)=&
      sqrt(abs(a_sp_mats(1_ik)%diagonal(omp_start_idx:omp_end_idx)))
  end if
  call cel_omp_shared_work_end_and_wait(cel_omp_shared_work,cel_omp,LOOP_SLEEP_DEFLENGTH)
  do ii=1_ik, size(a_sp_mats)
     select case (a_sp_mats(ii)%mv_algorithm)
       case(MV_CSR)
         call cel_sp_mat_distr_scale_type_csr(a_sp_mats(ii),a_sp_mats(1_ik)%scaled_diagonal,&
                          cel_omp,cel_omp_shared_work,&
                          cel_perf_counter,output_on,err_code)
       case(MV_JAD)
         call cel_sp_mat_distr_scale_type_csr(a_sp_mats(ii),a_sp_mats(1_ik)%scaled_diagonal,&
                          cel_omp,cel_omp_shared_work,&
                          cel_perf_counter,output_on,err_code)
       case(MV_COO)
        call cel_sp_mat_distr_scale_type_coo(a_sp_mats(ii),a_sp_mats(1_ik)%scaled_diagonal,&
                          cel_omp,cel_omp_shared_work,&
                          cel_perf_counter,output_on,err_code)
       case default
        call cel_error("Error cel_sp_mat_distr_scale: no mv_algorithm defined", &
                   err_code,cel_is_in_debug,.TRUE._lk)
    end select
  end do
  
  
end subroutine cel_sp_mat_distr_scale


!set jacobi diagonal matrix
subroutine cel_sp_mat_distr_set_jacobi(a_sp_mats,is_to_diag_scale,&
                            cel_omp,cel_omp_shared_work,cel_perf_counter, output_on,err_code)
  type(cel_sp_mat_type), dimension(:), allocatable,intent(inout)    :: a_sp_mats
  logical(kind=lk),intent(in)                                       :: is_to_diag_scale
  type(cel_omp_type), intent(inout)                                 :: cel_omp
  type(cel_omp_shared_work_type), intent(inout)                     :: cel_omp_shared_work
  type(cel_perf_counter_type), intent(inout)                        :: cel_perf_counter
  logical(kind=lk)  , intent(in)                                    :: output_on
  integer(kind=ik), intent(out)                                     :: err_code
  integer(kind=ik) :: start_idx,end_idx,omp_start_idx,omp_end_idx,num_threads
  integer(kind=ik) :: length,sys_stat
  err_code = 0_ik
  

  length = size(a_sp_mats(1_ik)%diagonal)
  call cel_omp_shared_work_start(cel_omp_shared_work,cel_omp%num_threads)
  if(cel_omp%is_master) then
    if(.not. allocated(a_sp_mats(1_ik)%jacobi)) then
      allocate(a_sp_mats(1_ik)%jacobi(length), stat=sys_stat)  
      call cel_error("Error cel_sp_mat_distr_set_jacobi alloc 1", sys_stat,&
                        cel_is_in_debug,.TRUE._lk)
      !write(*,'(A,I0)')"length:",length
    end if
  end if
  call cel_omp_shared_work_end_and_wait(cel_omp_shared_work,cel_omp,LOOP_SLEEP_DEFLENGTH)
   !write(*,'(A,I0)')"go thread: ",cel_omp%thread_num
  !!$omp barrier
  call cel_omp_shared_work_distr(cel_omp_shared_work,&
                                 omp_start_idx,omp_end_idx,num_threads,&
                                 length,length,&
                                 cel_omp,output_on,err_code)
  if(cel_omp%worker_num .le. num_threads) then
    if(is_to_diag_scale) then
        a_sp_mats(1_ik)%jacobi(omp_start_idx:omp_end_idx)=&
        a_sp_mats(1_ik)%scaled_diagonal(omp_start_idx:omp_end_idx) / &
        a_sp_mats(1_ik)%diagonal(omp_start_idx:omp_end_idx) 
    else
        a_sp_mats(1_ik)%jacobi(omp_start_idx:omp_end_idx)=&
        1.0_rk / &
        a_sp_mats(1_ik)%diagonal(omp_start_idx:omp_end_idx) 
    end if
  end if
  if(cel_omp%worker_num .eq. 1_ik) then
      cel_perf_counter%divide = cel_perf_counter%divide+length
      if(is_to_diag_scale) then
         cel_perf_counter%load = cel_perf_counter%load+&
        length*2_ik*rk
      else
        cel_perf_counter%load = cel_perf_counter%load+&
        length*1_ik*rk
      end if
      cel_perf_counter%store = cel_perf_counter%store+&
        length*rk
  end if
end subroutine cel_sp_mat_distr_set_jacobi


!scale boundaries
subroutine cel_sp_mat_distr_boundary_diag_scale(boundary_distr_vector, scaled_diagonal,&
                            cel_omp,cel_omp_shared_work,cel_perf_counter, output_on,err_code)
  type(cel_sp_mat_distr_vec_type), intent(inout)                    :: boundary_distr_vector
  real(kind=rk), dimension(:), intent(in)                           :: scaled_diagonal
  type(cel_omp_type), intent(inout)                                 :: cel_omp
  type(cel_omp_shared_work_type), intent(inout)                     :: cel_omp_shared_work
  type(cel_perf_counter_type), intent(inout)                        :: cel_perf_counter
  logical(kind=lk)  , intent(in)                                    :: output_on
  integer(kind=ik), intent(out)                                     :: err_code
  integer(kind=ik) :: length, ii
  real(kind=rk) :: scale_fct
  integer(kind=ik) :: start_idx,end_idx,omp_start_idx,omp_end_idx,num_threads
  err_code = 0_ik
  
  length = size(scaled_diagonal)
  call cel_omp_shared_work_distr(cel_omp_shared_work,&
                                 omp_start_idx,omp_end_idx,num_threads,&
                                 length,length,&
                                 cel_omp,output_on,err_code)
  if(cel_omp%worker_num .le. num_threads) then
    do ii=omp_start_idx, omp_end_idx
      scale_fct = scaled_diagonal(ii)
      boundary_distr_vector%values(ii,1_ik) = boundary_distr_vector%values(ii,1_ik) / scale_fct
    end do
  end if
  if(cel_omp%worker_num .eq. 1_ik) then
      cel_perf_counter%divide = cel_perf_counter%divide+length
      cel_perf_counter%load = cel_perf_counter%load+&
        length*1_ik*rk
      cel_perf_counter%store = cel_perf_counter%store+&
        length*rk
  end if
end subroutine cel_sp_mat_distr_boundary_diag_scale

!scale boundaries
subroutine cel_sp_mat_distr_vector_diag_scale_inv(distr_vector, scaled_diagonal,&
                            cel_omp,cel_omp_shared_work,cel_perf_counter, output_on,err_code)
  type(cel_sp_mat_distr_vec_type), intent(inout)                    :: distr_vector
  real(kind=rk), dimension(:), intent(in)                           :: scaled_diagonal
  type(cel_omp_type), intent(inout)                                 :: cel_omp
  type(cel_omp_shared_work_type), intent(inout)                     :: cel_omp_shared_work
  type(cel_perf_counter_type), intent(inout)                        :: cel_perf_counter
  logical(kind=lk)  , intent(in)                                    :: output_on
  integer(kind=ik), intent(out)                                     :: err_code
  integer(kind=ik) :: length, ii
  real(kind=rk) :: scale_fct
  integer(kind=ik) :: start_idx,end_idx,omp_start_idx,omp_end_idx,num_threads
  err_code = 0_ik
  
  length = size(scaled_diagonal)
  call cel_omp_shared_work_distr(cel_omp_shared_work,&
                                 omp_start_idx,omp_end_idx,num_threads,&
                                 length,length,&
                                 cel_omp,output_on,err_code)
  if(cel_omp%worker_num .le. num_threads) then
    do ii=omp_start_idx, omp_end_idx
      scale_fct = scaled_diagonal(ii)
      distr_vector%values(ii,1_ik) = distr_vector%values(ii,1_ik) * scale_fct
    end do
  end if
  if(cel_omp%worker_num .eq. 1_ik) then
      cel_perf_counter%mult = cel_perf_counter%mult+length
      cel_perf_counter%load = cel_perf_counter%load+&
        length*1_ik*rk
      cel_perf_counter%store = cel_perf_counter%store+&
        length*rk
  end if
end subroutine cel_sp_mat_distr_vector_diag_scale_inv

subroutine cel_sp_mat_distr_scale_type_csr(a_sp_mat,scaled_diagonal,&
                            cel_omp,cel_omp_shared_work,cel_perf_counter, output_on,err_code)
  type(cel_sp_mat_type), intent(inout)                              :: a_sp_mat
  real(kind=rk), dimension(:), allocatable, intent(in)              :: scaled_diagonal
  type(cel_omp_type), intent(inout)                                 :: cel_omp
  type(cel_omp_shared_work_type), intent(inout)                     :: cel_omp_shared_work
  type(cel_perf_counter_type), intent(inout)                        :: cel_perf_counter
  logical(kind=lk)  , intent(in)                                    :: output_on
  integer(kind=ik), intent(out)                                     :: err_code
  integer(kind=ik) :: ii, jj, nrow, col_idx
  integer(kind=ik) :: start_idx,end_idx,omp_start_idx,omp_end_idx,num_threads
  integer(kind=ik) :: col_offset
  real(kind=rk) :: scale_fct
  
  err_code=0_ik
  nrow = size(a_sp_mat%xadj) - 1_ik
  call cel_omp_shared_work_distr(cel_omp_shared_work,&
                                 omp_start_idx,omp_end_idx,num_threads,&
                                 nrow,a_sp_mat%elements_num,&
                                 cel_omp,output_on,err_code)
  if(cel_omp%worker_num .le. num_threads) then
    col_offset = a_sp_mat%col_offset
    do ii=omp_start_idx, omp_end_idx
      start_idx = a_sp_mat%xadj(ii)
      end_idx = a_sp_mat%xadj(ii+1_ik) - 1_ik
      scale_fct = scaled_diagonal(ii)
      a_sp_mat%values(start_idx:end_idx) = a_sp_mat%values(start_idx:end_idx)/scale_fct
    end do
  end if
  if(cel_omp%worker_num .eq. 1_ik) then
      cel_perf_counter%divide = cel_perf_counter%divide+&
        a_sp_mat%elements_num 
      cel_perf_counter%load = cel_perf_counter%load+&
        a_sp_mat%elements_num*rk+nrow*3_ik
      cel_perf_counter%store = cel_perf_counter%store+&
        a_sp_mat%elements_num*rk
      cel_perf_counter%sqr = cel_perf_counter%sqr+nrow
      cel_perf_counter%abso = cel_perf_counter%abso+nrow
  end if
end subroutine cel_sp_mat_distr_scale_type_csr


subroutine cel_sp_mat_distr_scale_type_coo(a_sp_mat,scaled_diagonal,&
                            cel_omp,cel_omp_shared_work,cel_perf_counter, output_on,err_code)
  type(cel_sp_mat_type),intent(inout)                               :: a_sp_mat
  real(kind=rk), dimension(:), allocatable, intent(in)              :: scaled_diagonal
  type(cel_omp_type), intent(inout)                                 :: cel_omp
  type(cel_omp_shared_work_type), intent(inout)                     :: cel_omp_shared_work
  type(cel_perf_counter_type), intent(inout)                        :: cel_perf_counter
  logical(kind=lk)  , intent(in)                                    :: output_on
  integer(kind=ik), intent(out)                                     :: err_code
  integer(kind=ik) :: ii
  integer(kind=ik) :: omp_start_idx,omp_end_idx,num_threads
  integer(kind=ik) :: row_offset, row_idx
  real(kind=rk) :: scale_fct

  err_code = 0_ik
  
  call cel_omp_shared_work_distr(cel_omp_shared_work,&
                                 omp_start_idx,omp_end_idx,num_threads,&
                                 a_sp_mat%elements_num,a_sp_mat%elements_num,&
                                 cel_omp,output_on,err_code)
  if(cel_omp%worker_num .le. num_threads) then
    row_offset =  a_sp_mat%row_offset
    do ii = omp_start_idx, omp_end_idx
      row_idx = a_sp_mat%rows(ii) - row_offset
      scale_fct = scaled_diagonal(row_idx)
      a_sp_mat%values(ii) = a_sp_mat%values(ii) / scale_fct
    end do
  end if
  if(cel_omp%worker_num .eq. 1_ik) then
      cel_perf_counter%divide = cel_perf_counter%divide+&
        a_sp_mat%elements_num 
      cel_perf_counter%load = cel_perf_counter%load+&
        a_sp_mat%elements_num*2_ik*rk+a_sp_mat%elements_num*ik
      cel_perf_counter%store = cel_perf_counter%store+&
        a_sp_mat%elements_num*rk
      cel_perf_counter%sqr = cel_perf_counter%sqr+a_sp_mat%elements_num
      cel_perf_counter%abso = cel_perf_counter%abso+a_sp_mat%elements_num
  end if
end subroutine cel_sp_mat_distr_scale_type_coo


end module cel_sp_mat_distr_format_module

