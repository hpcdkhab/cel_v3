!     
! File:  cel_blspmt_cel.f90
! Project CRESTA (see details on https://cresta-project.eu) Exascale library
! Send you email about the bugs to
! The High Performance Computing Center Stuttgart (HLRS) of the University of Stuttgart
! Dmitry Khabi,  (email: khabi@hlrs.de)
! Created on Aug 27, 2013
!

module cel_prec_module

use cel_types_module
use cel_input_parameters_module
use cel_omp_module
use cel_omp_shared_work_module
use cel_timer_interface_module
use cel_vec_module
use cel_sp_mat_module
use cel_input_parameters_module
!external libraries
use MPI
implicit none

contains


subroutine cel_prec_build_diagonal(sp_mat_local_array, &
   diagonal_size,cel_omp_shared_work,&
   cel_omp, input_parameters,err_code) 
  type(cel_sp_mat_type), dimension(:), allocatable, intent(inout) :: sp_mat_local_array
  integer(kind=ik), intent(in)                                    :: diagonal_size
  type(cel_omp_shared_work_type), intent(inout)                   :: cel_omp_shared_work
  type(cel_omp_type), intent(in)                                  :: cel_omp
  type(cel_input_parameters_type), intent(in)                     :: input_parameters
  integer (kind=ik), intent(out)                                  :: err_code
  integer(kind=ik) :: num_rows, num_diag_blocks, num_last_diag_block_rows
  integer(kind=ik) :: num_diag_block_el
  integer(kind=ik) :: num_last_diag_block_el, num_matrix_rows
  integer(kind=ik) :: omp_start_idx,omp_end_idx,num_threads, ii, jj, temp_int
  integer(kind=ik) :: start_idx,end_idx,col,row,col_offset,row_offset, cur_block_num, left_bound
  integer(kind=ik) :: row_num_in_block, col_in_block, diagonal_size_local, right_bound
  real(kind=rk)    :: value_tmp, block_num_corr
  logical(kind=lk) :: output_on
  err_code = 0_ik
  output_on = input_parameters%verbosity
  diagonal_size_local = diagonal_size
  if(diagonal_size_local .le. 0_ik) then
     call cel_error("Error cel_prec_build_diagonal diagonal_size_local <0 0 ", 1_ik,&
               cel_is_in_debug,.FALSE._lk)
      err_code = 1_ik
      return
  end if
  !how much blocks ?
  num_rows = sp_mat_local_array(1_ik)%rows_num
  if(diagonal_size_local .gt. num_rows) then
     call cel_error("Error cel_prec_build_diagonal diagonal_size_local > num_rows ", 2_ik,&
               cel_is_in_debug,.FALSE._lk)
      err_code = 2_ik
      return
  end if
  block_num_corr = 0_ik
  if(diagonal_size_local .eq. 1_ik)  block_num_corr = 1_ik
  num_diag_blocks = num_rows / diagonal_size_local
  num_last_diag_block_rows = num_rows-num_diag_blocks*diagonal_size_local
  if(num_last_diag_block_rows .gt. 0_ik) num_diag_blocks = num_diag_blocks + 1_ik
  !allocate  diagonal blocks
  if(cel_omp%is_master) then
    if(.not. sp_mat_local_array(1_ik)%is_in_csr) then
      call cel_error("Error cel_prec_build_diagonal sparse matrix is not in csr format", 3_ik,&
               cel_is_in_debug,.FALSE._lk)
      err_code = 3_ik
      return
    end if
    allocate(sp_mat_local_array(1_ik)%diagonal_blocks(diagonal_size_local,diagonal_size_local,num_diag_blocks))
    sp_mat_local_array(1_ik)%num_last_diag_block_rows = num_last_diag_block_rows
  end if
  
  !$omp barrier
  call cel_omp_shared_work_start(cel_omp_shared_work,cel_omp%num_threads)
  call cel_omp_shared_work_distr(cel_omp_shared_work,&
                                 omp_start_idx,omp_end_idx,num_threads,&
                                 num_diag_blocks,num_diag_blocks,&
                                 cel_omp,output_on,err_code)
  if(cel_omp%worker_num .le. num_threads) then
    do ii=omp_start_idx, omp_end_idx
      sp_mat_local_array(1_ik)%diagonal_blocks(1_ik:diagonal_size_local,1_ik:diagonal_size_local,ii) = 0.0_rk
    end do
  end if
  call cel_omp_shared_work_end_and_wait(cel_omp_shared_work,cel_omp,LOOP_SLEEP_DEFLENGTH)
  col_offset = sp_mat_local_array(1_ik)%col_offset
  row_offset = sp_mat_local_array(1_ik)%row_offset
  call cel_omp_shared_work_start(cel_omp_shared_work,cel_omp%num_threads)
  call cel_omp_shared_work_distr(cel_omp_shared_work,&
                                 omp_start_idx,omp_end_idx,num_threads,&
                                 num_rows,&
                                 num_rows,&
                                 cel_omp,output_on,err_code)
      !if(cel_omp%worker_num .eq. 1_ik .and. cel_omp%master_num .eq. 0_ik) then
      !    write(*,'(3(A,I0))') "omp_start_idx:",omp_start_idx,"; omp_end_idx:",omp_end_idx, &
      !    "; block_num_corr:",block_num_corr
      !end if
   if(cel_omp%worker_num .le. num_threads) then
   
    do ii=omp_start_idx, omp_end_idx
      start_idx = sp_mat_local_array(1_ik)%xadj(ii)
      end_idx = sp_mat_local_array(1_ik)%xadj(ii+1_ik) - 1_ik
      temp_int = mod(ii-1_ik,diagonal_size_local)
      row_num_in_block = temp_int + 1_ik
      cur_block_num = (ii-temp_int)/diagonal_size_local - block_num_corr + 1_ik
      right_bound = diagonal_size_local*cur_block_num
      left_bound =  diagonal_size_local*cur_block_num-diagonal_size_local+1_ik
      do jj=start_idx, end_idx
        col = sp_mat_local_array(1_ik)%cols(jj) - col_offset
        !if(cel_omp%worker_num .eq. 1_ik) then
        !  write(*,'(3(A,I0))') "col:",col,"; ge:",temp_int,"; le:",right_bound
        !end if
       ! if(cel_omp%master_num .eq. 0_ik  .and. &
       !    cur_block_num .eq. 1_ik ) then
       !   write(*,'(3(A,I0))') "row:",ii,"; col:",col, "; cur_block_num:", cur_block_num
       ! end if
        if(  col .ge. left_bound ) then
           if( col .le. right_bound ) then
             col_in_block = col - ii + row_num_in_block
             sp_mat_local_array(1_ik)%diagonal_blocks(row_num_in_block,col_in_block,cur_block_num) = &
               sp_mat_local_array(1_ik)%values(jj)
         !    if(cel_omp%master_num .eq. 0_ik .and. cur_block_num .eq. 1_ik) then
        !       write(*,'((3(A,I0)),(A,E16.6))') "row_num_in_block:",row_num_in_block,"; idx:",&
        !        idx,"; cur_block_num:",cur_block_num, "; value:",&
         !       sp_mat_local_array(1_ik)%diagonal_blocks(row_num_in_block,idx,cur_block_num)
         !    end if
           else
             exit
          end if
        end if
      end do
    end do
  end if
  call cel_omp_shared_work_end_and_wait(cel_omp_shared_work,cel_omp,LOOP_SLEEP_DEFLENGTH)

end subroutine cel_prec_build_diagonal

subroutine cel_prec_build_diagonal_inv(sp_mat_local_array, &
   cel_omp_shared_work,&
   cel_omp, input_parameters,err_code) 
  type(cel_sp_mat_type), dimension(:), allocatable, intent(inout) :: sp_mat_local_array
  type(cel_omp_shared_work_type), intent(inout)                   :: cel_omp_shared_work
  type(cel_omp_type), intent(in)                                  :: cel_omp
  type(cel_input_parameters_type), intent(in)                     :: input_parameters
  integer (kind=ik), intent(out)                                  :: err_code
  logical (kind=lk) :: output_on
  integer (kind=ik) :: diagonal_size_local
  integer (kind=ik) :: num_diag_blocks, dim1, dim2, reduce
  integer(kind=ik) :: omp_start_idx,omp_end_idx,num_threads, ii, sys_stat
  real(kind=rk), dimension(:,:), allocatable :: temp_block
  
  
  err_code = 0_ik
  output_on = input_parameters%verbosity
  dim1 = size(sp_mat_local_array(1_ik)%diagonal_blocks,dim=1_ik)
  dim2 = size(sp_mat_local_array(1_ik)%diagonal_blocks,dim=2_ik)
  num_diag_blocks = size(sp_mat_local_array(1_ik)%diagonal_blocks,dim=3)
  
  if(sp_mat_local_array(1_ik)%num_last_diag_block_rows .gt. 0_ik) then
    reduce = 1_ik
  else
    reduce = 0_ik
  end if
  if(dim1 .ne. dim2) then
  call cel_error("Error cel_prec_build_diagonal sparse matrix diagonal_blocks is not quadratic", 1_ik,&
               output_on,.FALSE._lk)
      err_code = 1_ik
      return
  end if
  if(cel_omp%is_master) then
    if(.not. allocated(sp_mat_local_array(1_ik)%diagonal_blocks)) then
      call cel_error("Error cel_prec_build_diagonal sparse matrix diagonal_blocks is not allocated", 2_ik,&
               output_on,.FALSE._lk)
      err_code = 2_ik
      return
    end if
    allocate(sp_mat_local_array(1_ik)%diagonal_inv_blocks(dim1,dim2,num_diag_blocks), stat=sys_stat)
    call cel_error("Error cel_sp_mat_distr_from_sp_mat_local  alloc 2", sys_stat,&
         output_on,.FALSE._lk)
    if(sys_stat .gt. 0) then
      err_code = 3_ik
      return
    end if
  end if
  !$omp barrier
  call cel_omp_shared_work_distr(cel_omp_shared_work,&
                                 omp_start_idx,omp_end_idx,num_threads,&
                                 num_diag_blocks-reduce,&
                                 num_diag_blocks-reduce,&
                                 cel_omp,output_on,err_code)
  if(cel_omp%worker_num .le. num_threads) then
    call cel_error("Error cel_sp_mat_distr_from_sp_mat_local  alloc 3", sys_stat,&
         output_on,.TRUE._lk)
    call cel_prec_inverse_matrix_serial(sp_mat_local_array,&
         omp_start_idx,omp_end_idx,dim1,dim1)
  end if  

  if(cel_omp%is_master) then
    if(sp_mat_local_array(1_ik)%num_last_diag_block_rows .gt. 0) then
     call cel_prec_inverse_matrix_serial(sp_mat_local_array,&
         num_diag_blocks,num_diag_blocks,dim1,&
         sp_mat_local_array(1_ik)%num_last_diag_block_rows)
    end if
  end if

  !$omp barrier

end subroutine cel_prec_build_diagonal_inv

subroutine cel_prec_inverse_matrix_serial(sp_mat_local_array,omp_start_idx,omp_end_idx,dim1,nn)
  type(cel_sp_mat_type), dimension(:), allocatable, intent(inout) :: sp_mat_local_array
  integer(kind=ik), intent(in)                :: omp_start_idx
  integer(kind=ik), intent(in)                :: omp_end_idx
  integer(kind=ik), intent(in)                :: dim1
  integer(kind=ik), intent(in)                :: nn
  integer(kind=ik) :: ii,mp_start_idx,omp_end_idx,start_idx,end_idx
  real(kind=rk),dimension(:,:), allocatable:: temp_block

  allocate(temp_block(dim1,dim1))
  select case (nn)
    case (1_ik)
      do ii = omp_start_idx,omp_end_idx
        start_idx = (ii-1_ik) + 1_ik
        end_idx = ii
        temp_block(start_idx:end_idx,start_idx:end_idx) = &
          sp_mat_local_array(1_ik)%diagonal_blocks(start_idx:end_idx,&
          start_idx:end_idx,ii)
        call cel_prec_inverse_matrix_serial1_1x1(temp_block,&
          sp_mat_local_array(1_ik)%diagonal_inv_blocks(start_idx:end_idx,&
          start_idx:end_idx,ii),1_ik)
      end do
   case default
     do ii = omp_start_idx,omp_end_idx
        start_idx = (ii-1_ik)*dim1 + 1_ik
        end_idx = ii*dim1
        temp_block(start_idx:end_idx,start_idx:end_idx) = &
          sp_mat_local_array(1_ik)%diagonal_blocks(start_idx:end_idx,&
          start_idx:end_idx,ii)
        call cel_prec_inverse_matrix_serial1_nxn(temp_block,&
          sp_mat_local_array(1_ik)%diagonal_inv_blocks(start_idx:end_idx,&
          start_idx:end_idx,ii),nn)
     end do
 end select

end subroutine cel_prec_inverse_matrix_serial

subroutine cel_prec_inverse_matrix_serial1_1x1(aa,cc,nn)
real(kind=rk),dimension(1_ik,1_ik), intent(inout) :: aa
real(kind=rk),dimension(1_ik,1_ik), intent(inout) :: cc
integer(kind=ik), intent(in)                :: nn

 cc(1_ik,1_ik)=1.0_rk/aa(1_ik,1_ik)

end subroutine cel_prec_inverse_matrix_serial1_1x1

subroutine cel_prec_inverse_matrix_serial1_nxn(a,c,n)

!============================================================
! Inverse matrix
! Method: Based on Doolittle LU factorization for Ax=b
! Alex G. December 2009
!-----------------------------------------------------------
! input ...
! a(n,n) - array of coefficients for matrix A
! n      - dimension
! output ...
! c(n,n) - inverse matrix of A
! comments ...
! the original matrix a(n,n) will be destroyed 
! during the calculation
!===========================================================
implicit none 
integer(kind=ik) n
real(kind=rk) :: a(n,n), c(n,n)
real(kind=rk) :: L(n,n), U(n,n), b(n), d(n), x(n)
real(kind=rk) :: coeff
integer(kind=ik) :: i, j, k

! step 0: initialization for matrices L and U and b
! Fortran 90/95 aloows such operations on matrices
L=0.0
U=0.0
b=0.0

! step 1: forward elimination
do k=1, n-1
   do i=k+1,n
      coeff=a(i,k)/a(k,k)
      L(i,k) = coeff
      do j=k+1,n
         a(i,j) = a(i,j)-coeff*a(k,j)
      end do
   end do
end do

! Step 2: prepare L and U matrices 
! L matrix is a matrix of the elimination coefficient
! + the diagonal elements are 1.0
do i=1,n
  L(i,i) = 1.0
end do
! U matrix is the upper triangular part of A
do j=1,n
  do i=1,j
    U(i,j) = a(i,j)
  end do
end do

! Step 3: compute columns of the inverse matrix C
do k=1,n
  b(k)=1.0
  d(1) = b(1)
! Step 3a: Solve Ld=b using the forward substitution
  do i=2,n
    d(i)=b(i)
    do j=1,i-1
      d(i) = d(i) - L(i,j)*d(j)
    end do
  end do
! Step 3b: Solve Ux=d using the back substitution
  x(n)=d(n)/U(n,n)
  do i = n-1,1,-1
    x(i) = d(i)
    do j=n,i+1,-1
      x(i)=x(i)-U(i,j)*x(j)
    end do
    x(i) = x(i)/u(i,i)
  end do
! Step 3c: fill the solutions x(n) into column k of C
  do i=1,n
    c(i,k) = x(i)
  end do
  b(k)=0.0
end do
end subroutine cel_prec_inverse_matrix_serial1_nxn

end module cel_prec_module
