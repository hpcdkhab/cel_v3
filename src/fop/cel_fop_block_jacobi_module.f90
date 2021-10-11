! Definition of the distributed matrix vector multiplication
! File:   cel_fop_vec_distr_module.f90
! Project CRESTA Exascale library (see details on https://www.cresta-project.eu)
! Send you email about the bugs to
! The High Performance Computing Center Stuttgart (HLRS) of the University of Stuttgart
! Dmitry Khabi,  (email: khabi@hlrs.de)
! Created on Jun 14, 2013
!

module cel_fop_block_jacobi_module
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



!vector_p(i) = vector_r(i)*block_jacobi(i)
subroutine cel_fop_block_jacobi_vec_mult(vector_p,vector_r,diagonal_inv_blocks,&
                                num_last_diag_block_rows,&
                                cel_omp,cel_omp_shared_work,&
                                cel_perf_counter, output_on,err_code)
  type(cel_sp_mat_distr_vec_type), intent(inout)                    :: vector_p
  type(cel_sp_mat_distr_vec_type), intent(in)                       :: vector_r
  real(kind=rk), dimension(:,:,:), intent(in)                       :: diagonal_inv_blocks
   integer(kind=ik), intent(in)                                     :: num_last_diag_block_rows
  type(cel_omp_type), intent(in)                                    :: cel_omp
  type(cel_omp_shared_work_type), intent(inout)                     :: cel_omp_shared_work
  type(cel_perf_counter_type), intent(inout)                        :: cel_perf_counter
  logical(kind=lk)  , intent(in)                                    :: output_on
  integer(kind=ik), intent(out)                                     :: err_code
  integer(kind=ik) :: vec_length
  real(kind=trk) :: time
  integer(kind=ik) :: omp_start_idx,omp_end_idx,num_threads, start_idx,end_idx
  integer(kind=ik) :: dim1, dim2, dim3, parallel_blocks, ii

  err_code = 0_ik
  if(cel_omp%worker_num .eq. 1_ik) then
   time = get_time()
  end if
  dim1 = size(diagonal_inv_blocks,dim=1)
  dim3 = size(diagonal_inv_blocks,dim=3)
  vec_length = vector_p%num_elements
  
  if(num_last_diag_block_rows .gt. 0) then 
    parallel_blocks = dim3
  else
    parallel_blocks = dim3 - 1_ik
  end if
  call cel_omp_shared_work_distr(cel_omp_shared_work,&
                                 omp_start_idx,omp_end_idx,num_threads,&
                                 parallel_blocks,parallel_blocks,&
                                 cel_omp,output_on,err_code)

  if(cel_omp%worker_num .le. num_threads) then
    call cel_fop_vec_mult_distr_jacobi_block( vector_p%values(1_ik:vec_length,1_ik),&
       vector_r%values(1_ik:vec_length,1_ik), &
       diagonal_inv_blocks,&
       dim1,dim1,omp_start_idx, omp_end_idx)

  end if
  if(cel_omp%is_master .and. num_last_diag_block_rows .gt. 0_ik) then
     call cel_fop_vec_mult_distr_jacobi_block(vector_p%values(1_ik:vec_length,1_ik),&
      vector_r%values(1_ik:vec_length,1_ik),&
      diagonal_inv_blocks,dim1,&
       num_last_diag_block_rows,dim3,dim3)
  end if
  if(cel_omp%worker_num .eq. 1_ik) then
    cel_perf_counter%time(1_ik) = cel_perf_counter%time(1_ik) + (get_time() - time)
    cel_perf_counter%store= cel_perf_counter%store + &
    dim3*parallel_blocks*rk+num_last_diag_block_rows*rk
    cel_perf_counter%load= cel_perf_counter%store + dim3*dim3*parallel_blocks*rk+&
      num_last_diag_block_rows*num_last_diag_block_rows*rk
    cel_perf_counter%mult= cel_perf_counter%mult + (dim3*dim3)*parallel_blocks*rk+&
      (num_last_diag_block_rows*num_last_diag_block_rows)*rk
    cel_perf_counter%add= cel_perf_counter%add + (dim3*dim3-1_ik)*parallel_blocks*rk+&
      (num_last_diag_block_rows*num_last_diag_block_rows-1_ik)*rk
  end if
  
end subroutine cel_fop_block_jacobi_vec_mult

!subroutine for  one block multiplication 
subroutine cel_fop_vec_mult_distr_jacobi_block(vector_p,vector_r,diagonal_inv_blocks,&
  dim1,num_rows,omp_start_idx,omp_end_idx)

  real(kind=rk), dimension(:), intent(inout)                 :: vector_p
  real(kind=rk), dimension(:), intent(in)                    :: vector_r
  real(kind=rk), dimension(:,:,:), intent(in)                :: diagonal_inv_blocks
  integer(kind=ik), intent(in)                               :: dim1
  integer(kind=ik), intent(in)                               :: num_rows
  integer(kind=ik), intent(in)                               :: omp_start_idx
  integer(kind=ik), intent(in)                               :: omp_end_idx
  integer(kind=ik) :: ii, start_idx, end_idx
   select case (num_rows)
     case (1_ik)
       do ii = omp_start_idx,omp_end_idx
         start_idx = (ii-1_ik)*dim1 + 1_ik
         end_idx = ii*dim1
         call cel_fop_vec_mult_distr_jacobi_block_1x1(&
           vector_p(start_idx:end_idx),&
           vector_r(start_idx:end_idx),&
           diagonal_inv_blocks(start_idx:end_idx,start_idx:end_idx,ii))
       end do
     case (2_ik)
       do ii = omp_start_idx,omp_end_idx
         start_idx = (ii-1_ik)*dim1 + 1_ik
         end_idx = ii*dim1
         call cel_fop_vec_mult_distr_jacobi_block_2x2(&
           vector_p(start_idx:end_idx),&
           vector_r(start_idx:end_idx),&
           diagonal_inv_blocks(start_idx:end_idx,start_idx:end_idx,ii))
       end do
      case (3_ik)
        do ii = omp_start_idx,omp_end_idx
          start_idx = (ii-1_ik)*dim1 + 1_ik
          end_idx = ii*dim1
          call cel_fop_vec_mult_distr_jacobi_block_3x3(&
            vector_p(start_idx:end_idx),&
            vector_r(start_idx:end_idx),&
            diagonal_inv_blocks(start_idx:end_idx,start_idx:end_idx,ii))
       end do
      case (4_ik)
        do ii = omp_start_idx,omp_end_idx
          start_idx = (ii-1_ik)*dim1 + 1_ik
          end_idx = ii*dim1
          call cel_fop_vec_mult_distr_jacobi_block_4x4(&
            vector_p(start_idx:end_idx),&
            vector_r(start_idx:end_idx),&
            diagonal_inv_blocks(start_idx:end_idx,start_idx:end_idx,ii))
       end do
      case (5_ik)
        do ii = omp_start_idx,omp_end_idx
          start_idx = (ii-1_ik)*dim1 + 1_ik
          end_idx = ii*dim1
          call cel_fop_vec_mult_distr_jacobi_block_5x5(&
            vector_p(start_idx:end_idx),&
            vector_r(start_idx:end_idx),&
            diagonal_inv_blocks(start_idx:end_idx,start_idx:end_idx,ii))
       end do
      case (6_ik)
        do ii = omp_start_idx,omp_end_idx
          start_idx = (ii-1_ik)*dim1 + 1_ik
          end_idx = ii*dim1
          call cel_fop_vec_mult_distr_jacobi_block_6x6(&
            vector_p(start_idx:end_idx),&
            vector_r(start_idx:end_idx),&
            diagonal_inv_blocks(start_idx:end_idx,start_idx:end_idx,ii))
       end do
      case (7_ik)
        do ii = omp_start_idx,omp_end_idx
          start_idx = (ii-1_ik)*dim1 + 1_ik
          end_idx = ii*dim1
          call cel_fop_vec_mult_distr_jacobi_block_7x7(&
            vector_p(start_idx:end_idx),&
            vector_r(start_idx:end_idx),&
            diagonal_inv_blocks(start_idx:end_idx,start_idx:end_idx,ii))
       end do
      case (8_ik)
        do ii = omp_start_idx,omp_end_idx
          start_idx = (ii-1_ik)*dim1 + 1_ik
          end_idx = ii*dim1
          call cel_fop_vec_mult_distr_jacobi_block_8x8(&
            vector_p(start_idx:end_idx),&
            vector_r(start_idx:end_idx),&
            diagonal_inv_blocks(start_idx:end_idx,start_idx:end_idx,ii))
       end do
     case default
        do ii = omp_start_idx,omp_end_idx
          start_idx = (ii-1_ik)*dim1 + 1_ik
          end_idx = ii*dim1
          call cel_fop_vec_mult_distr_jacobi_block_nxn(&
            vector_p(start_idx:end_idx),&
            vector_r(start_idx:end_idx),&
            diagonal_inv_blocks(start_idx:end_idx,start_idx:end_idx,ii),num_rows)
       end do
    end select
  
end subroutine cel_fop_vec_mult_distr_jacobi_block
!subroutine for general (any size nxn) one block multiplication 
!vector_p = vector_r*block_jacobi
subroutine cel_fop_vec_mult_distr_jacobi_block_nxn(vector_p,vector_r,diagonal_inv_block,num_rows)  
  real(kind=rk), dimension(:), intent(inout)                 :: vector_p
  real(kind=rk), dimension(:), intent(in)                    :: vector_r
  real(kind=rk), dimension(:,:), intent(in)                  :: diagonal_inv_block
  integer(kind=ik), intent(in)                               :: num_rows
  integer(kind=ik) :: ii, jj

  vector_p(1_ik:num_rows) = 0_rk
  do ii=1_ik,num_rows
    do jj=1_ik,num_rows
      vector_p(ii) = vector_p(ii) + vector_r(ii)*diagonal_inv_block(ii,jj)
    end do
  end do

end subroutine cel_fop_vec_mult_distr_jacobi_block_nxn

!subroutine for  one block (size 1x1) multiplication 
!vector_p = vector_r*block_jacobi
subroutine cel_fop_vec_mult_distr_jacobi_block_1x1(vector_p,vector_r,diagonal_inv_block)   
  real(kind=rk), dimension(1_ik), intent(inout)                 :: vector_p
  real(kind=rk), dimension(1_ik), intent(in)                    :: vector_r
  real(kind=rk), dimension(1_ik,1_ik), intent(in)                  :: diagonal_inv_block

  vector_p(1_ik) = vector_r(1_ik)*diagonal_inv_block(1_ik,1_ik)

end subroutine cel_fop_vec_mult_distr_jacobi_block_1x1

!subroutine for  one block (size 2x2) multiplication 
!vector_p = vector_r*block_jacobi
subroutine cel_fop_vec_mult_distr_jacobi_block_2x2(vector_p,vector_r,diagonal_inv_block)   
  real(kind=rk), dimension(2_ik), intent(inout)                 :: vector_p
  real(kind=rk), dimension(2_ik), intent(in)                    :: vector_r
  real(kind=rk), dimension(2_ik,2_ik), intent(in)                  :: diagonal_inv_block

  vector_p(1_ik) = vector_r(1_ik)*diagonal_inv_block(1_ik,1_ik)+&
    vector_r(2_ik)*diagonal_inv_block(1_ik,2_ik)
  vector_p(2_ik) = vector_r(1_ik)*diagonal_inv_block(2_ik,1_ik)+&
    vector_r(2_ik)*diagonal_inv_block(2_ik,2_ik)
end subroutine cel_fop_vec_mult_distr_jacobi_block_2x2

!subroutine for  one block (size 3x3) multiplication 
!vector_p = vector_r*block_jacobi
subroutine cel_fop_vec_mult_distr_jacobi_block_3x3(vector_p,vector_r,diagonal_inv_block)   
  real(kind=rk), dimension(3_ik), intent(inout)                 :: vector_p
  real(kind=rk), dimension(3_ik), intent(in)                    :: vector_r
  real(kind=rk), dimension(3_ik,3_ik), intent(in)                  :: diagonal_inv_block

  vector_p(1_ik) = vector_r(1_ik)*diagonal_inv_block(1_ik,1_ik)+&
    vector_r(2_ik)*diagonal_inv_block(1_ik,2_ik)+vector_r(3_ik)*diagonal_inv_block(1_ik,3_ik)
  vector_p(2_ik) = vector_r(1_ik)*diagonal_inv_block(2_ik,1_ik)+&
    vector_r(2_ik)*diagonal_inv_block(2_ik,2_ik)+vector_r(3_ik)*diagonal_inv_block(2_ik,3_ik)
  vector_p(3_ik) = vector_r(1_ik)*diagonal_inv_block(3_ik,1_ik)+&
    vector_r(2_ik)*diagonal_inv_block(3_ik,2_ik)+vector_r(3_ik)*diagonal_inv_block(3_ik,3_ik)
end subroutine cel_fop_vec_mult_distr_jacobi_block_3x3

!subroutine for  one block (size 4x4) multiplication 
!vector_p = vector_r*block_jacobi
subroutine cel_fop_vec_mult_distr_jacobi_block_4x4(vector_p,vector_r,diagonal_inv_block)   
  real(kind=rk), dimension(4_ik), intent(inout)                 :: vector_p
  real(kind=rk), dimension(4_ik), intent(in)                    :: vector_r
  real(kind=rk), dimension(4_ik,4_ik), intent(in)                  :: diagonal_inv_block

  vector_p(1_ik) = vector_r(1_ik)*diagonal_inv_block(1_ik,1_ik)+&
    vector_r(2_ik)*diagonal_inv_block(1_ik,2_ik)+vector_r(3_ik)*diagonal_inv_block(1_ik,3_ik)+&
    vector_r(4_ik)*diagonal_inv_block(1_ik,4_ik)
  vector_p(2_ik) = vector_r(1_ik)*diagonal_inv_block(2_ik,1_ik)+&
    vector_r(2_ik)*diagonal_inv_block(2_ik,2_ik)+vector_r(3_ik)*diagonal_inv_block(2_ik,3_ik)+&
    vector_r(4_ik)*diagonal_inv_block(2_ik,4_ik)
  vector_p(3_ik) = vector_r(1_ik)*diagonal_inv_block(3_ik,1_ik)+&
    vector_r(2_ik)*diagonal_inv_block(3_ik,2_ik)+vector_r(3_ik)*diagonal_inv_block(3_ik,3_ik)+&
    vector_r(4_ik)*diagonal_inv_block(3_ik,4_ik)
  vector_p(4_ik) = vector_r(1_ik)*diagonal_inv_block(4_ik,1_ik)+&
    vector_r(2_ik)*diagonal_inv_block(4_ik,2_ik)+vector_r(3_ik)*diagonal_inv_block(4_ik,3_ik)+&
    vector_r(4_ik)*diagonal_inv_block(4_ik,4_ik)
    
end subroutine cel_fop_vec_mult_distr_jacobi_block_4x4

!subroutine for  one block (size 5x5) multiplication 
!vector_p = vector_r*block_jacobi
subroutine cel_fop_vec_mult_distr_jacobi_block_5x5(vector_p,vector_r,diagonal_inv_block)   
  real(kind=rk), dimension(5_ik), intent(inout)                 :: vector_p
  real(kind=rk), dimension(5_ik), intent(in)                    :: vector_r
  real(kind=rk), dimension(5_ik,5_ik), intent(in)                  :: diagonal_inv_block

  vector_p(1_ik) = vector_r(1_ik)*diagonal_inv_block(1_ik,1_ik)+&
    vector_r(2_ik)*diagonal_inv_block(1_ik,2_ik)+vector_r(3_ik)*diagonal_inv_block(1_ik,3_ik)+&
    vector_r(4_ik)*diagonal_inv_block(1_ik,4_ik)+vector_r(5_ik)*diagonal_inv_block(1_ik,5_ik)
  vector_p(2_ik) = vector_r(1_ik)*diagonal_inv_block(1_ik,1_ik)+&
    vector_r(2_ik)*diagonal_inv_block(2_ik,2_ik)+vector_r(3_ik)*diagonal_inv_block(2_ik,3_ik)+&
    vector_r(4_ik)*diagonal_inv_block(2_ik,4_ik)+vector_r(5_ik)*diagonal_inv_block(2_ik,5_ik)
  vector_p(3_ik) = vector_r(1_ik)*diagonal_inv_block(3_ik,1_ik)+&
    vector_r(2_ik)*diagonal_inv_block(3_ik,2_ik)+vector_r(3_ik)*diagonal_inv_block(3_ik,3_ik)+&
    vector_r(4_ik)*diagonal_inv_block(3_ik,4_ik)+vector_r(5_ik)*diagonal_inv_block(3_ik,5_ik)
  vector_p(4_ik) = vector_r(1_ik)*diagonal_inv_block(4_ik,1_ik)+&
    vector_r(2_ik)*diagonal_inv_block(4_ik,2_ik)+vector_r(3_ik)*diagonal_inv_block(4_ik,3_ik)+&
    vector_r(4_ik)*diagonal_inv_block(4_ik,4_ik)+vector_r(5_ik)*diagonal_inv_block(4_ik,5_ik)
  vector_p(5_ik) = vector_r(1_ik)*diagonal_inv_block(5_ik,1_ik)+&
    vector_r(2_ik)*diagonal_inv_block(5_ik,2_ik)+vector_r(3_ik)*diagonal_inv_block(5_ik,3_ik)+&
    vector_r(4_ik)*diagonal_inv_block(5_ik,4_ik)+vector_r(5_ik)*diagonal_inv_block(5_ik,5_ik)
    
end subroutine cel_fop_vec_mult_distr_jacobi_block_5x5

!subroutine for  one block (size 6x6) multiplication 
!vector_p = vector_r*block_jacobi
subroutine cel_fop_vec_mult_distr_jacobi_block_6x6(vector_p,vector_r,diagonal_inv_block)   
  real(kind=rk), dimension(6_ik), intent(inout)                 :: vector_p
  real(kind=rk), dimension(6_ik), intent(in)                    :: vector_r
  real(kind=rk), dimension(6_ik,6_ik), intent(in)                  :: diagonal_inv_block

  vector_p(1_ik) = vector_r(1_ik)*diagonal_inv_block(1_ik,1_ik)+&
    vector_r(2_ik)*diagonal_inv_block(1_ik,2_ik)+vector_r(3_ik)*diagonal_inv_block(1_ik,3_ik)+&
    vector_r(4_ik)*diagonal_inv_block(1_ik,4_ik)+vector_r(5_ik)*diagonal_inv_block(1_ik,5_ik)+&
    vector_r(6_ik)*diagonal_inv_block(1_ik,6_ik)
  vector_p(2_ik) = vector_r(1_ik)*diagonal_inv_block(1_ik,1_ik)+&
    vector_r(2_ik)*diagonal_inv_block(2_ik,2_ik)+vector_r(3_ik)*diagonal_inv_block(2_ik,3_ik)+&
    vector_r(4_ik)*diagonal_inv_block(2_ik,4_ik)+vector_r(5_ik)*diagonal_inv_block(2_ik,5_ik)+&
    vector_r(6_ik)*diagonal_inv_block(2_ik,6_ik)
  vector_p(3_ik) = vector_r(1_ik)*diagonal_inv_block(3_ik,1_ik)+&
    vector_r(2_ik)*diagonal_inv_block(3_ik,2_ik)+vector_r(3_ik)*diagonal_inv_block(3_ik,3_ik)+&
    vector_r(4_ik)*diagonal_inv_block(3_ik,4_ik)+vector_r(5_ik)*diagonal_inv_block(3_ik,5_ik)+&
    vector_r(6_ik)*diagonal_inv_block(3_ik,6_ik)
  vector_p(4_ik) = vector_r(1_ik)*diagonal_inv_block(4_ik,1_ik)+&
    vector_r(2_ik)*diagonal_inv_block(4_ik,2_ik)+vector_r(3_ik)*diagonal_inv_block(4_ik,3_ik)+&
    vector_r(4_ik)*diagonal_inv_block(4_ik,4_ik)+vector_r(5_ik)*diagonal_inv_block(4_ik,5_ik)+&
    vector_r(6_ik)*diagonal_inv_block(4_ik,6_ik)
  vector_p(5_ik) = vector_r(1_ik)*diagonal_inv_block(5_ik,1_ik)+&
    vector_r(2_ik)*diagonal_inv_block(5_ik,2_ik)+vector_r(3_ik)*diagonal_inv_block(5_ik,3_ik)+&
    vector_r(4_ik)*diagonal_inv_block(5_ik,4_ik)+vector_r(5_ik)*diagonal_inv_block(5_ik,5_ik)+&
    vector_r(6_ik)*diagonal_inv_block(5_ik,6_ik)
  vector_p(6_ik) = vector_r(1_ik)*diagonal_inv_block(6_ik,1_ik)+&
    vector_r(2_ik)*diagonal_inv_block(6_ik,2_ik)+vector_r(3_ik)*diagonal_inv_block(6_ik,3_ik)+&
    vector_r(4_ik)*diagonal_inv_block(6_ik,4_ik)+vector_r(5_ik)*diagonal_inv_block(6_ik,5_ik)+&
    vector_r(6_ik)*diagonal_inv_block(6_ik,6_ik)
    
end subroutine cel_fop_vec_mult_distr_jacobi_block_6x6

!subroutine for  one block (size 7x7) multiplication 
!vector_p = vector_r*block_jacobi
subroutine cel_fop_vec_mult_distr_jacobi_block_7x7(vector_p,vector_r,diagonal_inv_block)   
  real(kind=rk), dimension(7_ik), intent(inout)                 :: vector_p
  real(kind=rk), dimension(7_ik), intent(in)                    :: vector_r
  real(kind=rk), dimension(7_ik,7_ik), intent(in)                  :: diagonal_inv_block

  vector_p(1_ik) = vector_r(1_ik)*diagonal_inv_block(1_ik,1_ik)+&
    vector_r(2_ik)*diagonal_inv_block(1_ik,2_ik)+vector_r(3_ik)*diagonal_inv_block(1_ik,3_ik)+&
    vector_r(4_ik)*diagonal_inv_block(1_ik,4_ik)+vector_r(5_ik)*diagonal_inv_block(1_ik,5_ik)+&
    vector_r(6_ik)*diagonal_inv_block(1_ik,6_ik)+vector_r(7_ik)*diagonal_inv_block(1_ik,7_ik)
  vector_p(2_ik) = vector_r(1_ik)*diagonal_inv_block(1_ik,1_ik)+&
    vector_r(2_ik)*diagonal_inv_block(2_ik,2_ik)+vector_r(3_ik)*diagonal_inv_block(2_ik,3_ik)+&
    vector_r(4_ik)*diagonal_inv_block(2_ik,4_ik)+vector_r(5_ik)*diagonal_inv_block(2_ik,5_ik)+&
    vector_r(6_ik)*diagonal_inv_block(2_ik,6_ik)+vector_r(7_ik)*diagonal_inv_block(2_ik,7_ik)
  vector_p(3_ik) = vector_r(1_ik)*diagonal_inv_block(3_ik,1_ik)+&
    vector_r(2_ik)*diagonal_inv_block(3_ik,2_ik)+vector_r(3_ik)*diagonal_inv_block(3_ik,3_ik)+&
    vector_r(4_ik)*diagonal_inv_block(3_ik,4_ik)+vector_r(5_ik)*diagonal_inv_block(3_ik,5_ik)+&
    vector_r(6_ik)*diagonal_inv_block(3_ik,6_ik)+vector_r(7_ik)*diagonal_inv_block(3_ik,7_ik)
  vector_p(4_ik) = vector_r(1_ik)*diagonal_inv_block(4_ik,1_ik)+&
    vector_r(2_ik)*diagonal_inv_block(4_ik,2_ik)+vector_r(3_ik)*diagonal_inv_block(4_ik,3_ik)+&
    vector_r(4_ik)*diagonal_inv_block(4_ik,4_ik)+vector_r(5_ik)*diagonal_inv_block(4_ik,5_ik)+&
    vector_r(6_ik)*diagonal_inv_block(4_ik,6_ik)+vector_r(7_ik)*diagonal_inv_block(4_ik,7_ik)
  vector_p(5_ik) = vector_r(1_ik)*diagonal_inv_block(5_ik,1_ik)+&
    vector_r(2_ik)*diagonal_inv_block(5_ik,2_ik)+vector_r(3_ik)*diagonal_inv_block(5_ik,3_ik)+&
    vector_r(4_ik)*diagonal_inv_block(5_ik,4_ik)+vector_r(5_ik)*diagonal_inv_block(5_ik,5_ik)+&
    vector_r(6_ik)*diagonal_inv_block(5_ik,6_ik)+vector_r(7_ik)*diagonal_inv_block(5_ik,7_ik)
  vector_p(6_ik) = vector_r(1_ik)*diagonal_inv_block(6_ik,1_ik)+&
    vector_r(2_ik)*diagonal_inv_block(6_ik,2_ik)+vector_r(3_ik)*diagonal_inv_block(6_ik,3_ik)+&
    vector_r(4_ik)*diagonal_inv_block(6_ik,4_ik)+vector_r(5_ik)*diagonal_inv_block(6_ik,5_ik)+&
    vector_r(6_ik)*diagonal_inv_block(6_ik,6_ik)+vector_r(7_ik)*diagonal_inv_block(6_ik,7_ik)
  vector_p(7_ik) = vector_r(1_ik)*diagonal_inv_block(7_ik,1_ik)+&
    vector_r(2_ik)*diagonal_inv_block(7_ik,2_ik)+vector_r(3_ik)*diagonal_inv_block(7_ik,3_ik)+&
    vector_r(4_ik)*diagonal_inv_block(7_ik,4_ik)+vector_r(5_ik)*diagonal_inv_block(7_ik,5_ik)+&
    vector_r(6_ik)*diagonal_inv_block(7_ik,6_ik)+vector_r(7_ik)*diagonal_inv_block(7_ik,7_ik)
end subroutine cel_fop_vec_mult_distr_jacobi_block_7x7

!subroutine for  one block (size 7x7) multiplication 
!vector_p = vector_r*block_jacobi
subroutine cel_fop_vec_mult_distr_jacobi_block_8x8(vector_p,vector_r,diagonal_inv_block)   
  real(kind=rk), dimension(8_ik), intent(inout)                 :: vector_p
  real(kind=rk), dimension(8_ik), intent(in)                    :: vector_r
  real(kind=rk), dimension(8_ik,8_ik), intent(in)                  :: diagonal_inv_block

  vector_p(1_ik) = vector_r(1_ik)*diagonal_inv_block(1_ik,1_ik)+&
    vector_r(2_ik)*diagonal_inv_block(1_ik,2_ik)+vector_r(3_ik)*diagonal_inv_block(1_ik,3_ik)+&
    vector_r(4_ik)*diagonal_inv_block(1_ik,4_ik)+vector_r(5_ik)*diagonal_inv_block(1_ik,5_ik)+&
    vector_r(6_ik)*diagonal_inv_block(1_ik,6_ik)+vector_r(7_ik)*diagonal_inv_block(1_ik,7_ik)+&
    vector_r(8_ik)*diagonal_inv_block(1_ik,8_ik)
  vector_p(2_ik) = vector_r(1_ik)*diagonal_inv_block(1_ik,1_ik)+&
    vector_r(2_ik)*diagonal_inv_block(2_ik,2_ik)+vector_r(3_ik)*diagonal_inv_block(2_ik,3_ik)+&
    vector_r(4_ik)*diagonal_inv_block(2_ik,4_ik)+vector_r(5_ik)*diagonal_inv_block(2_ik,5_ik)+&
    vector_r(6_ik)*diagonal_inv_block(2_ik,6_ik)+vector_r(7_ik)*diagonal_inv_block(2_ik,7_ik)+&
    vector_r(8_ik)*diagonal_inv_block(2_ik,8_ik)
  vector_p(3_ik) = vector_r(1_ik)*diagonal_inv_block(3_ik,1_ik)+&
    vector_r(2_ik)*diagonal_inv_block(3_ik,2_ik)+vector_r(3_ik)*diagonal_inv_block(3_ik,3_ik)+&
    vector_r(4_ik)*diagonal_inv_block(3_ik,4_ik)+vector_r(5_ik)*diagonal_inv_block(3_ik,5_ik)+&
    vector_r(6_ik)*diagonal_inv_block(3_ik,6_ik)+vector_r(7_ik)*diagonal_inv_block(3_ik,7_ik)+&
    vector_r(8_ik)*diagonal_inv_block(3_ik,8_ik)
  vector_p(4_ik) = vector_r(1_ik)*diagonal_inv_block(4_ik,1_ik)+&
    vector_r(2_ik)*diagonal_inv_block(4_ik,2_ik)+vector_r(3_ik)*diagonal_inv_block(4_ik,3_ik)+&
    vector_r(4_ik)*diagonal_inv_block(4_ik,4_ik)+vector_r(5_ik)*diagonal_inv_block(4_ik,5_ik)+&
    vector_r(6_ik)*diagonal_inv_block(4_ik,6_ik)+vector_r(7_ik)*diagonal_inv_block(4_ik,7_ik)+&
    vector_r(8_ik)*diagonal_inv_block(4_ik,8_ik)
  vector_p(5_ik) = vector_r(1_ik)*diagonal_inv_block(5_ik,1_ik)+&
    vector_r(2_ik)*diagonal_inv_block(5_ik,2_ik)+vector_r(3_ik)*diagonal_inv_block(5_ik,3_ik)+&
    vector_r(4_ik)*diagonal_inv_block(5_ik,4_ik)+vector_r(5_ik)*diagonal_inv_block(5_ik,5_ik)+&
    vector_r(6_ik)*diagonal_inv_block(5_ik,6_ik)+vector_r(7_ik)*diagonal_inv_block(5_ik,7_ik)+&
    vector_r(8_ik)*diagonal_inv_block(5_ik,8_ik)
  vector_p(6_ik) = vector_r(1_ik)*diagonal_inv_block(6_ik,1_ik)+&
    vector_r(2_ik)*diagonal_inv_block(6_ik,2_ik)+vector_r(3_ik)*diagonal_inv_block(6_ik,3_ik)+&
    vector_r(4_ik)*diagonal_inv_block(6_ik,4_ik)+vector_r(5_ik)*diagonal_inv_block(6_ik,5_ik)+&
    vector_r(6_ik)*diagonal_inv_block(6_ik,6_ik)+vector_r(7_ik)*diagonal_inv_block(6_ik,7_ik)+&
    vector_r(8_ik)*diagonal_inv_block(6_ik,8_ik)
  vector_p(7_ik) = vector_r(1_ik)*diagonal_inv_block(7_ik,1_ik)+&
    vector_r(2_ik)*diagonal_inv_block(7_ik,2_ik)+vector_r(3_ik)*diagonal_inv_block(7_ik,3_ik)+&
    vector_r(4_ik)*diagonal_inv_block(7_ik,4_ik)+vector_r(5_ik)*diagonal_inv_block(7_ik,5_ik)+&
    vector_r(6_ik)*diagonal_inv_block(7_ik,6_ik)+vector_r(7_ik)*diagonal_inv_block(7_ik,7_ik)+&
    vector_r(8_ik)*diagonal_inv_block(7_ik,8_ik)
  vector_p(8_ik) = vector_r(1_ik)*diagonal_inv_block(8_ik,1_ik)+&
    vector_r(2_ik)*diagonal_inv_block(8_ik,2_ik)+vector_r(3_ik)*diagonal_inv_block(8_ik,3_ik)+&
    vector_r(4_ik)*diagonal_inv_block(8_ik,4_ik)+vector_r(5_ik)*diagonal_inv_block(8_ik,5_ik)+&
    vector_r(6_ik)*diagonal_inv_block(8_ik,6_ik)+vector_r(7_ik)*diagonal_inv_block(8_ik,7_ik)+&
    vector_r(8_ik)*diagonal_inv_block(8_ik,8_ik)
end subroutine cel_fop_vec_mult_distr_jacobi_block_8x8

end module cel_fop_block_jacobi_module

