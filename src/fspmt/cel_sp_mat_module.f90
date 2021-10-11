!     
! File:   cel_sp_mat_module.f90
! Project CRESTA (see details on https://www.cresta-project.eu) Exascale library
! Send you email about the bugs to
! The High Performance Computing Center Stuttgart (HLRS) of the University of Stuttgart
! Dmitry Khabi,  (email: khabi@hlrs.de)
! Created on May 29, 2013
!
! Distributed Sparse matrix in Coordinate format.
! One process owns one block 
! Distribution is done by rows

module cel_sp_mat_module
use cel_types_module
use cel_algorithms_param_module
use cel_base_print_module
use cel_error_module
use cel_timer_interface_module
use OMP_LIB
use, intrinsic :: ISO_C_BINDING
implicit none

! Interface for  struct "sparse matrix" in c library libcelbase.a
type, bind(c) :: cel_c_sp_mat_type
  integer(kind=c_ik) :: block_indx ! block_indx block index of this part of matrix
  integer(kind=c_ik) ::  elements_num! elements_num number of elements
  integer(kind=c_ik) ::  rows_num! rows_num number of rows
  integer(kind=c_ik) ::  cols_num ! cols_num number of columns
  integer(kind=c_ik) ::  block_rows_num! rows_num number of rows in block 
  integer(kind=c_ik) ::  block_cols_num! cols_num number of columns in block
  type(c_ptr) :: rows! row index array of matrix
  type(c_ptr) :: cols! cols index array of matrix
  type(c_ptr) :: values! data array of matrix
end type cel_c_sp_mat_type


!Sparse matrix type 
type cel_sp_mat_type
  integer(kind=ik):: block_indx ! block_indx block index of this part of matrix
  integer(kind=ik):: elements_num ! elements_num number of elements
  integer(kind=ik):: rows_num; ! rows_num number of rows
  integer(kind=c_ik) ::  cols_num ! cols_num number of columns
  integer(kind=c_ik) ::  block_rows_num! rows_num number of rows in block 
  integer(kind=c_ik) ::  block_cols_num! cols_num number of columns in block
  integer(kind=ik) :: col_offset ! offset to the begin of the referenced part of the ditributed vector X (Ax=y)
  integer(kind=ik) :: row_offset ! offset to the begin of the referenced part of the ditributed vector Y (Ax=y)
  integer(kind=ik), dimension(:), allocatable :: rows ! row index array of matrix
  integer(kind=ik), dimension(:), allocatable :: cols ! col index array of matrix
  real(kind=rk), dimension(:), allocatable :: values ! data array of matrix
  real(kind=rk), dimension(:), allocatable :: jacobi ! data array of matrix
  real(kind=rk), dimension(:), allocatable :: diagonal! columnsvector of A given by the max norm
  real(kind=rk), dimension(:), allocatable :: scaled_diagonal! diagonal matrix for matrix diagonal scaling
  real(kind=rk), dimension(:,:,:), allocatable :: diagonal_blocks
  integer(kind=ik)                             :: num_last_diag_block_rows
  real(kind=rk), dimension(:,:,:), allocatable :: diagonal_inv_blocks
  integer(kind=ik), dimension(:), allocatable :: xadj !array with row pointers in csr format
  integer(kind=ik) :: jad_nmax !jad format
  integer(kind=ik) :: jad_mmax !jad format
  integer(kind=ik), allocatable, dimension(:) :: jad_begin !jad format
  integer(kind=ik), allocatable, dimension(:) :: jad_cols !jad format
  integer(kind=ik), allocatable, dimension(:) :: jad_perm !jad format
  real(kind=rk), allocatable, dimension(:) :: jad_values !jad format
  real(kind=rk), allocatable, dimension(:) :: jad_tmp_values !jad format
  real(kind=rk), allocatable, dimension(:) :: jad_tmp_reductions_values !jad format
  integer(kind=ik):: mv_algorithm ! MV_CSR or MV_JAD
  integer(kind=ik):: mv_communicator ! MV_COMM_ISEND or MV_COMM_IBCAST
  logical(kind=lk):: is_in_coo
  logical(kind=lk):: is_in_csr
  logical(kind=lk):: is_in_jad
  logical(kind=lk):: is_scaled
end type cel_sp_mat_type

interface cel_convert_c
  module procedure cel_sp_mat_c_to_f
end interface cel_convert_c

interface print
  module procedure cel_sp_mat_print
end interface print

interface new
  module procedure cel_sp_mat_new
end interface new

interface del
  module procedure cel_sp_mat_del
  module procedure cel_sp_mat_array_del
end interface del

interface del_jad
  module procedure cel_sp_mat_del_jad
end interface del_jad

interface object_size_byte
  module procedure cel_sp_mat_size_bytes
end interface object_size_byte

interface copy
  module procedure cel_sp_mat_copy
end interface copy


contains

subroutine cel_sp_mat_print(string,sp_mat,full_out)
  character(len=*) :: string
  type(cel_sp_mat_type), intent(inout)  :: sp_mat
  logical(kind=lk), intent(in)          :: full_out
  integer(kind=ik) :: ii
  integer(kind=ik) :: el_counter_diag, el_counter_up, el_counter_under
  integer(kind=ik) :: row_num, elements_num
  row_num = sp_mat%rows_num
  elements_num = sp_mat%rows_num
  el_counter_diag = 1_ik
  el_counter_up = row_num+1_ik
  el_counter_under = row_num + &
    (elements_num)/2_ik+1_ik

  write(cel_output_unit,'(2a)') trim(string)
  write(cel_output_unit,'(a,i0)') "cols_num:", sp_mat%cols_num
  write(cel_output_unit,'(a,i0)') "rows_num:", sp_mat%rows_num
  write(cel_output_unit,'(a,i0)') "elements_num:", sp_mat%elements_num
  write(cel_output_unit,'(a,i0)') "block_indx:", sp_mat%block_indx
  write(cel_output_unit,'(a,i0)') "block_rows_num:", sp_mat%block_rows_num
  write(cel_output_unit,'(a,i0)') "block_cols_num:", sp_mat%block_cols_num
  print*
  if(full_out) then
    print*,'idx,row,column,value'
    do ii=el_counter_diag,sp_mat%elements_num
       if(ii .EQ. el_counter_up) then
        write(cel_output_unit,'(A)') "=========start second part==========="
       end if
       if(ii .EQ. el_counter_under) then
        write(cel_output_unit,'(A)') "=========start third part==========="
       end if
       write(cel_output_unit,'(1x,i0,1x,i0,1x,i0,1x,E13.6)') ii, sp_mat%rows(ii),&
       sp_mat%cols(ii), sp_mat%values(ii)
    enddo
    print*
    print*
    
    call cel_sp_mat_print_jad( "=========jad matrix data===========",sp_mat)
  end if
 
end subroutine cel_sp_mat_print

subroutine cel_sp_mat_print_jad(string,sp_mat)
  character(len=*) :: string
  type(cel_sp_mat_type), intent(inout)  :: sp_mat
  integer(kind=ik) :: mm,nn
  integer(kind=ik),allocatable,dimension(:)          :: length
  
  write(cel_output_unit,'(2a)') trim(string)
  if(sp_mat%jad_nmax .gt. 0_ik .and. sp_mat%jad_mmax .gt. 0_ik .and. &
    allocated(sp_mat%jad_perm) .and. allocated(sp_mat%jad_begin) .and. &
    allocated(sp_mat%jad_cols) .and. allocated(sp_mat%jad_values)) then
    write(cel_output_unit,'(a,i0)') "nmax:", sp_mat%jad_nmax
    write(cel_output_unit,'(a,i0)') "mmax:", sp_mat%jad_mmax
    call print("perm",sp_mat%jad_perm)
    call print("begin",sp_mat%jad_begin)
    call print("cols",sp_mat%jad_cols)
    !call print("values",sparse_matrix_jad%values)
    allocate(length(sp_mat%jad_mmax))
    do mm=1,sp_mat%jad_mmax
       length(mm) = sp_mat%jad_begin(mm+1) - sp_mat%jad_begin(mm)
    enddo
    write (cel_output_unit,*)
    write (cel_output_unit,'(A)') 'begin,length'
    do mm=1,sp_mat%jad_mmax
       write(cel_output_unit,'(1x,i0,1x,i0)') sp_mat%jad_begin(mm),length(mm)
    enddo
    print*
    print*
    do mm=1,sp_mat%jad_mmax
      write(cel_output_unit,'(1x,i0,"|")',advance='no') mm
      do nn=sp_mat%jad_begin(mm),sp_mat%jad_begin(mm+1)-1_ik
        !if(nn < sparse_matrix_jad%begin(mm+1) .AND. &
        !  nn >= sparse_matrix_jad%begin(mm)) then
        write(cel_output_unit,'(1x,i0)',advance='no') sp_mat%jad_cols(nn)
        !endif
      enddo
      write(cel_output_unit,*)
    enddo
    deallocate(length)
  else
    write(cel_output_unit,'(A)') "!!!!!!!!no jad matrix data!!!!!!!!!!!"
  end if

end subroutine cel_sp_mat_print_jad

function cel_sp_mat_jad_size_bytes(sp_mat) result(bytes)
  type(cel_sp_mat_type), intent(inout)  :: sp_mat
  integer(kind=ik) :: bytes
  bytes = 0_ik
  if(sp_mat%jad_nmax .gt. 0_ik .and. sp_mat%jad_mmax .gt. 0_ik .and. &
    allocated(sp_mat%jad_perm) .and. allocated(sp_mat%jad_begin) .and. &
    allocated(sp_mat%jad_cols) .and. allocated(sp_mat%jad_values)) then
      bytes = bytes+size(sp_mat%jad_perm)*ik
      bytes = bytes+size(sp_mat%jad_begin)*ik
      bytes = bytes+size(sp_mat%jad_cols)*ik
      bytes = bytes+size(sp_mat%jad_values)*rk
      if(allocated(sp_mat%jad_tmp_values)) then
        bytes = bytes+size(sp_mat%jad_tmp_values,dim=1)*rk
      end if
   else
     write(cel_output_unit,'(A)') "!!!!!!!!no jad matrix data!!!!!!!!!!!"
  end if
end function cel_sp_mat_jad_size_bytes

function cel_sp_mat_size_bytes(sp_mat) result(bytes)
  type(cel_sp_mat_type), intent(inout)  :: sp_mat
  integer(kind=ik) :: bytes

  bytes = 0_ik
  if( allocated(sp_mat%rows) ) bytes = bytes + size(sp_mat%rows)*ik
  if( allocated(sp_mat%cols) ) bytes = bytes + size(sp_mat%cols)*ik
  if( allocated(sp_mat%values) ) bytes = bytes + size(sp_mat%values)*rk
  if( allocated(sp_mat%diagonal) ) bytes = bytes + size(sp_mat%diagonal)*rk
  if( allocated(sp_mat%scaled_diagonal) ) bytes = bytes + size(sp_mat%scaled_diagonal)*rk
  if( allocated(sp_mat%diagonal_inv_blocks) ) then
   bytes = bytes + size(sp_mat%diagonal_inv_blocks,dim=1_ik)*&
    size(sp_mat%diagonal_inv_blocks,dim=2_ik)*size(sp_mat%diagonal_inv_blocks,dim=3_ik)
  end if
  if( allocated(sp_mat%diagonal_blocks) ) then
    bytes = bytes + size(sp_mat%diagonal_blocks,dim=1_ik)*&
    size(sp_mat%diagonal_blocks,dim=2_ik)*size(sp_mat%diagonal_blocks,dim=3_ik)
  end if
  bytes = bytes + cel_sp_mat_jad_size_bytes(sp_mat)
  
end function cel_sp_mat_size_bytes

! convert c object to f object
subroutine cel_sp_mat_c_to_f(c_sp_mat, sp_mat, to_deallocate_c)
  type(cel_c_sp_mat_type), intent(inout) :: c_sp_mat
  type(cel_sp_mat_type), intent(inout) :: sp_mat
  logical(kind=lk), intent(in) :: to_deallocate_c
  real(kind=rk), pointer :: c_real_array(:)
  integer(kind=ik), pointer :: c_int_array(:)
  integer(kind=c_ik) :: ierr
  integer(kind=ik) :: sys_stat
  
  sp_mat%block_indx=c_sp_mat%block_indx
  sp_mat%elements_num=c_sp_mat%elements_num
  sp_mat%rows_num=c_sp_mat%rows_num
  sp_mat%cols_num=c_sp_mat%cols_num
  sp_mat%block_rows_num=c_sp_mat%block_rows_num
  sp_mat%block_cols_num=c_sp_mat%block_cols_num
  call c_f_pointer(c_sp_mat%values,c_real_array,[c_sp_mat%elements_num])
  allocate(sp_mat%values(c_sp_mat%elements_num), stat=sys_stat)
  call cel_error("cel_sp_mat_c_to_f c_real_array", sys_stat, .TRUE._lk,.TRUE._lk)
  if (associated(c_real_array)) sp_mat%values = c_real_array
  if(to_deallocate_c) then
    ierr = cel_free_double_c(c_real_array(1))
  endif

  call c_f_pointer(c_sp_mat%rows,c_int_array,[c_sp_mat%elements_num])
  allocate(sp_mat%rows(c_sp_mat%elements_num), stat=sys_stat)
  call cel_error("cel_sp_mat_c_to_f c_int_array 1", sys_stat, .TRUE._lk,.TRUE._lk)
  sp_mat%rows = c_int_array
  if(to_deallocate_c) then
    ierr = cel_free_int_c(c_int_array(1))
  endif
       
  call c_f_pointer(c_sp_mat%cols,c_int_array,[c_sp_mat%elements_num])
  allocate(sp_mat%cols(c_sp_mat%elements_num), stat=sys_stat)
  call cel_error("cel_sp_mat_c_to_f c_int_array 2", sys_stat, .TRUE._lk,.TRUE._lk)
  sp_mat%cols = c_int_array
  if(to_deallocate_c)  then
    ierr = cel_free_int_c(c_int_array(1))
  endif 

end subroutine cel_sp_mat_c_to_f

!initializing of new cel_sp_mat_type
subroutine cel_sp_mat_new(sp_mat, rows_num, cols_num,&
  block_rows_num, block_cols_num ,element_num)
  type(cel_sp_mat_type), intent(inout)  :: sp_mat
  integer(kind=ik), intent(in)  :: rows_num
  integer(kind=ik), intent(in)  :: cols_num
  integer(kind=ik), intent(in)  :: block_rows_num
  integer(kind=ik), intent(in)  :: block_cols_num
  integer(kind=ik), intent(in)  :: element_num
  integer (kind=ik):: err_code
  
  if( allocated(sp_mat%rows) ) deallocate(sp_mat%rows)
  if( allocated(sp_mat%cols) ) deallocate(sp_mat%cols)
  if( allocated(sp_mat%values) ) deallocate(sp_mat%values)
  if( allocated(sp_mat%diagonal) ) deallocate(sp_mat%diagonal)

  if( allocated(sp_mat%jacobi) ) deallocate(sp_mat%jacobi)
  if( allocated(sp_mat%xadj) ) deallocate(sp_mat%xadj)
  if( allocated(sp_mat%jad_begin) ) deallocate(sp_mat%jad_begin)
  if( allocated(sp_mat%jad_cols) ) deallocate(sp_mat%jad_cols)
  if( allocated(sp_mat%jad_perm) ) deallocate(sp_mat%jad_perm)
  if( allocated(sp_mat%jad_values) ) deallocate(sp_mat%jad_values)
  if( allocated(sp_mat%jad_tmp_values) ) deallocate(sp_mat%jad_tmp_values)
  if( allocated(sp_mat%jad_tmp_reductions_values) ) deallocate(sp_mat%jad_tmp_reductions_values)
  if( allocated(sp_mat%diagonal_inv_blocks) ) deallocate(sp_mat%diagonal_inv_blocks)
  if( allocated(sp_mat%diagonal_blocks) ) deallocate(sp_mat%diagonal_blocks)

  sp_mat%block_indx = 0_ik
  sp_mat%elements_num = element_num
  sp_mat%rows_num = rows_num
  sp_mat%cols_num = cols_num
  sp_mat%block_rows_num = block_rows_num
  sp_mat%block_cols_num = block_cols_num
  sp_mat%col_offset = 0_ik
  sp_mat%row_offset = 0_ik
  sp_mat%is_in_coo = .FALSE._lk
  sp_mat%is_in_csr = .FALSE._lk
  sp_mat%is_in_jad = .FALSE._lk
  sp_mat%is_scaled = .FALSE._lk
  sp_mat%num_last_diag_block_rows = 0_ik
  allocate(sp_mat%rows(element_num), stat=err_code)
  call cel_error("cel_sp_mat_new rows", err_code, cel_is_in_debug, .TRUE._lk)
  allocate(sp_mat%cols(element_num), stat=err_code)
  call cel_error("cel_sp_mat_new cols", err_code, cel_is_in_debug, .TRUE._lk)
  allocate(sp_mat%values(element_num), stat=err_code)
  call cel_error("cel_sp_mat_new values", err_code, cel_is_in_debug, .TRUE._lk)
end subroutine cel_sp_mat_new

!copy cel_sp_mat_type 
subroutine cel_sp_mat_copy(sp_mat_copy, sp_mat_orig)
  type(cel_sp_mat_type), intent(inout)  :: sp_mat_copy
  type(cel_sp_mat_type), intent(in)  :: sp_mat_orig
  integer(kind=ik) :: dim1, dim2,dim3
  call del(sp_mat_copy)
  call new(sp_mat_copy,&
    sp_mat_orig%rows_num,&
    sp_mat_orig%cols_num,& 
    sp_mat_orig%block_rows_num,&
    sp_mat_orig%block_cols_num,&
    sp_mat_orig%elements_num)
  sp_mat_copy%rows = sp_mat_orig%rows
  sp_mat_copy%cols = sp_mat_orig%cols
  sp_mat_copy%values = sp_mat_orig%values
  sp_mat_copy%diagonal = sp_mat_orig%diagonal
  sp_mat_copy%num_last_diag_block_rows = sp_mat_orig%num_last_diag_block_rows
  if (allocated(sp_mat_orig%jacobi)) then
    if (.not. allocated(sp_mat_copy%jacobi)) then
      allocate(sp_mat_copy%scaled_diagonal(size(sp_mat_orig%jacobi)))
    end if
    sp_mat_copy%jacobi = sp_mat_orig%jacobi
  end if
  if (allocated(sp_mat_orig%scaled_diagonal)) then
    if (.not. allocated(sp_mat_copy%scaled_diagonal)) then
      allocate(sp_mat_copy%scaled_diagonal(size(sp_mat_orig%scaled_diagonal)))
    end if
    sp_mat_copy%scaled_diagonal = sp_mat_orig%scaled_diagonal
  end if
  if (allocated(sp_mat_orig%xadj)) then
    allocate(sp_mat_copy%xadj(size(sp_mat_orig%xadj)))
     if (.not. allocated(sp_mat_copy%xadj)) then
      allocate(sp_mat_copy%xadj(size(sp_mat_orig%xadj)))
    end if
    sp_mat_copy%xadj = sp_mat_orig%xadj
  end if
  if( allocated(sp_mat_orig%diagonal_inv_blocks) ) then
    if (.not. allocated(sp_mat_copy%diagonal_inv_blocks)) then
      dim1 = size(sp_mat_orig%diagonal_inv_blocks, dim=1_ik)
      dim2 = size(sp_mat_orig%diagonal_inv_blocks, dim=2_ik)
      dim3 = size(sp_mat_orig%diagonal_inv_blocks, dim=3_ik)
      allocate(sp_mat_copy%diagonal_inv_blocks( &
        dim1,dim2,dim3))
    end if
    sp_mat_copy%diagonal_inv_blocks = sp_mat_orig%diagonal_inv_blocks
  end if
   if( allocated(sp_mat_orig%diagonal_blocks) ) then
    if (.not. allocated(sp_mat_copy%diagonal_blocks)) then
      dim1 = size(sp_mat_orig%diagonal_blocks, dim=1_ik)
      dim2 = size(sp_mat_orig%diagonal_blocks, dim=2_ik)
      dim3 = size(sp_mat_orig%diagonal_blocks, dim=3_ik)
      allocate(sp_mat_copy%diagonal_blocks(&
       dim1,dim2,dim3))
    end if
    sp_mat_copy%diagonal_blocks = sp_mat_orig%diagonal_blocks
  end if
end subroutine cel_sp_mat_copy





!deallocate cel_sp_mat_type' arrays
!all values are set to undefined
subroutine cel_sp_mat_del(sp_mat)
  type(cel_sp_mat_type), intent(inout)  :: sp_mat
  sp_mat%block_indx = cel_index_undef
  sp_mat%elements_num = cel_index_undef
  sp_mat%rows_num = cel_index_undef
  sp_mat%cols_num = cel_index_undef
  if( allocated(sp_mat%rows) ) deallocate(sp_mat%rows)
  if( allocated(sp_mat%cols) ) deallocate(sp_mat%cols)
  if( allocated(sp_mat%values) ) deallocate(sp_mat%values)
  if( allocated(sp_mat%values) ) deallocate(sp_mat%diagonal)
  if( allocated(sp_mat%scaled_diagonal) ) deallocate(sp_mat%scaled_diagonal)
  if( allocated(sp_mat%xadj) ) deallocate(sp_mat%xadj)
  if( allocated(sp_mat%jacobi) ) deallocate(sp_mat%jacobi)
  if( allocated(sp_mat%diagonal_inv_blocks) ) deallocate(sp_mat%diagonal_inv_blocks)
  if( allocated(sp_mat%diagonal_blocks) ) deallocate(sp_mat%diagonal_blocks)
  call cel_sp_mat_del_jad(sp_mat)
end subroutine cel_sp_mat_del

!deallocate cel_sp_mat_type' arrays
!all values are set to undefined
subroutine cel_sp_mat_del_jad(sp_mat)
  type(cel_sp_mat_type), intent(inout)  :: sp_mat
  sp_mat%jad_nmax = 0_ik
  sp_mat%jad_mmax = 0_ik

  if( allocated(sp_mat%jad_begin) ) deallocate(sp_mat%jad_begin)
  if( allocated(sp_mat%jad_cols) ) deallocate(sp_mat%jad_cols)
  if( allocated(sp_mat%jad_perm) ) deallocate(sp_mat%jad_perm)
  if( allocated(sp_mat%jad_values) ) deallocate(sp_mat%jad_values)
end subroutine cel_sp_mat_del_jad


subroutine cel_sp_mat_array_del(sp_mats)
  type(cel_sp_mat_type), dimension(:), allocatable,intent(inout)  :: sp_mats
  integer(kind=ik) :: ii
  
  if(allocated(sp_mats)) then
    do ii=1_ik, size(sp_mats)
      call del(sp_mats(ii))
    end do
    deallocate(sp_mats)
  end if
end subroutine cel_sp_mat_array_del

subroutine cel_sp_mat_to_csr(sp_mat)
  type(cel_sp_mat_type), intent(inout)  :: sp_mat
  integer(kind=ik) :: curr_row, last_row, idx, ii, curr_idx
  integer(kind=ik), dimension(:), allocatable :: num_elements_in_rows
  
  if(allocated(sp_mat%xadj)) then
     if(size(sp_mat%xadj) .lt. sp_mat%rows_num+1_ik) then
       deallocate(sp_mat%xadj)
     end if
  end if
  if(.not. allocated(sp_mat%xadj)) then
     allocate(sp_mat%xadj(sp_mat%rows_num+1_ik))
  end if
  allocate(num_elements_in_rows(sp_mat%rows_num))
  num_elements_in_rows = 0_ik
  do ii=1_ik,sp_mat%elements_num
    curr_row = sp_mat%rows(ii) - sp_mat%row_offset
    num_elements_in_rows(curr_row) = num_elements_in_rows(curr_row) + 1_ik
  end do
  
  sp_mat%xadj(1_ik)=1_ik
  do ii=2_ik,sp_mat%rows_num + 1_ik
    sp_mat%xadj(ii)=sp_mat%xadj(ii-1_ik)+num_elements_in_rows(ii - 1_ik)
  end do
  !sp_mat%cols=sp_mat%cols-sp_mat%col_offset

end subroutine cel_sp_mat_to_csr



!sparse_31 - FMPS Format for symm. matrix (quadratic)
!description of format
!as csr except that  the size of row_ptr is equal to the number of rows,
!because there aren't any element in the last row except on one diagonal element
!only the upper triangular part of the matrix is stored (see what is local variable element_num)
!the values in sp_mat are not sorted, the first part consists  of diagonal elements,
!the second part consists of upper triangular values
!the third part (see definition of el_counter_under) is the transposed part of upper triangulare values

subroutine cel_sp_mat_from_format_sparse_31_full(sp_mat, matrix_diagonal,&
    matrix_up_triangl, column_indices, row_ptr,&
    block_rows_num, block_cols_num, output_on, err_code)
  type(cel_sp_mat_type), intent(inout)  :: sp_mat
  REAL(kind = rk), ALLOCATABLE, DIMENSION(:), intent(in) :: matrix_diagonal
  REAL(kind = rk), ALLOCATABLE, DIMENSION(:), intent(in) :: matrix_up_triangl
  INTEGER(kind = ik), ALLOCATABLE, DIMENSION(:), intent(in) :: column_indices
  INTEGER(kind = ik), ALLOCATABLE, DIMENSION(:), intent(in) :: row_ptr
  integer (kind=ik), intent(in) :: block_rows_num
  integer (kind=ik), intent(in) :: block_cols_num
  logical (kind=lk), intent(in) :: output_on
  integer (kind=ik), intent(out) :: err_code
  integer (kind=ik) :: ii, jj, row_num, el_counter_up, el_counter_under, el_counter_diag
  integer (kind=ik) :: local_row_ptr_length_minus, first_index, last_index, column_index, element_num
  real(kind=rk) :: column_value, diagonal_value
  err_code = 0_ik
  local_row_ptr_length_minus = INT(size(row_ptr),kind=ik) - 1_ik
  row_num = size(matrix_diagonal,dim=1)
  element_num = size(column_indices,dim=1)*2_ik+row_num
  el_counter_diag = 1_ik
  el_counter_up = row_num+1_ik
  el_counter_under = row_num + size(column_indices,dim=1)+1_ik
  call new(sp_mat, row_num, row_num, block_rows_num, block_cols_num, element_num)
  DO ii = 1, local_row_ptr_length_minus
      first_index = row_ptr(ii)
      last_index = row_ptr(ii + 1) - 1
      diagonal_value = matrix_diagonal(ii)
      sp_mat%rows(el_counter_diag)=ii
      sp_mat%cols(el_counter_diag)=ii
      sp_mat%values(el_counter_diag)=diagonal_value
      el_counter_diag = el_counter_diag + 1_ik
      DO jj = first_index, last_index
        column_index = column_indices(jj)
        column_value = matrix_up_triangl(ii)
        sp_mat%rows(el_counter_up)=ii
        sp_mat%cols(el_counter_up)=column_index
        sp_mat%values(el_counter_up)=column_value
        el_counter_up = el_counter_up + 1_ik
        sp_mat%rows(el_counter_under)=column_index
        sp_mat%cols(el_counter_under)=ii
        sp_mat%values(el_counter_under)=column_value
        el_counter_under = el_counter_under + 1_ik
      END DO
  END DO
  sp_mat%elements_num = element_num
  sp_mat%rows_num = row_num
  sp_mat%cols_num = row_num
end subroutine cel_sp_mat_from_format_sparse_31_full

subroutine cel_sp_mat_format_sparse_31_sort(sp_mat, output_on, err_code)
  type(cel_sp_mat_type), intent(inout)  :: sp_mat
  logical (kind=lk), intent(in) :: output_on
  integer (kind=ik), intent(out) :: err_code
  integer(kind=ik) start_row,  last_row, cur_row, cur_index, row_counter
  integer(kind=ik) start_index, end_index
  integer(kind = ik) :: start_cycl_rows=0, end_cycl_rows=0
  integer(kind = ik) :: start_cycl_cols=0, end_cycl_cols=0
  
  err_code = 0_ik
  
  !sort nach rows
  start_cycl_rows = getrdtsc()
  
  call cel_sp_mat_merge_sort(sp_mat%rows,&
  sp_mat%cols,&
  sp_mat%values,&
  sp_mat%elements_num, 0_ik, output_on)
  end_cycl_rows  = getrdtsc()
  if(output_on) then
    write(*,'(A,E17.9,A)') "cel_sp_mat_format_sparse_31_sort(sort order rows): ", &
    real(end_cycl_rows-start_cycl_rows,kind=rk), " cycles;"
    write(*,'(A,E17.9,A)') "cel_sp_mat_format_sparse_31_sort(sort order rows): ", &
    real(end_cycl_rows-start_cycl_rows,kind=rk)*cel_tsc_cycles, &
     " seconds;"
  end if
  !sort nach cols
  cur_index = 1_ik
  last_row = 1_ik
  row_counter = 0_ik
  start_index = 1_ik
  end_index = 1_ik
  start_cycl_cols = getrdtsc()
  do while(row_counter .LT. sp_mat%rows_num-1_ik)
    cur_row = sp_mat%rows(cur_index)
    !write(*,'(2(A,I0))') "cur_row:",cur_row,"; last_row:",last_row
    if(last_row .NE. cur_row .OR. cur_index .EQ. sp_mat%elements_num) then
      if(last_row .GT. cur_row) then
        call cel_error("cel_sp_mat_format_sparse_31_sort rows not sorted",&
         1_ik, output_on, .TRUE._lk)
      end if 
      start_row  = last_row
      last_row = cur_row
      if(cur_index .EQ. sp_mat%elements_num) end_index = cur_index
      if(start_index .LT. end_index) then
        !write(*,'(2(A,I0))') "start_index:",start_index,"; end_index:",end_index
        call cel_sp_mat_merge_sort(sp_mat%cols(start_index:end_index),&
          sp_mat%rows(start_index:end_index),&
          sp_mat%values(start_index:end_index),&
          end_index-start_index + 1_ik, 0_ik, output_on)
      end if
      start_index = end_index + 1_ik
      row_counter = row_counter + 1_ik
      last_row = cur_row
    end if
    end_index = cur_index
    cur_index = cur_index + 1_ik
  end do
  end_cycl_cols = getrdtsc()
  if(output_on) then
    write(*,'(A,E17.9,A)') "cel_sp_mat_format_sparse_31_sort(sort order cols): ", &
    real(end_cycl_cols-start_cycl_cols,kind=rk), " cycles;"
    write(*,'(A,E17.9,A)') "cel_sp_mat_format_sparse_31_sort(sort order cols): ", &
    real(end_cycl_cols-start_cycl_cols,kind=rk)*cel_tsc_cycles, &
     " seconds;"
  end if
  
end subroutine


subroutine cel_sp_mat_bubble_sort(order_array,array_int,array_real,&
 local_length)
 
  integer(kind=ik), dimension(:),  intent(inout) :: order_array
  integer(kind=ik), dimension(:),  intent(inout) :: array_int
  real(kind=rk), dimension(:) :: array_real
  integer(kind=ik),intent(in)  :: local_length
  real(kind=rk) :: temp_real
  integer(kind=ik) :: temp, temp_int
  logical(kind=lk) :: swapped = .TRUE._lk
  integer(kind=ik) :: ii, jj

  do jj = local_length-1_ik, 1_ik, -1_ik
    swapped = .FALSE._lk
    do ii = 1, jj
      if(order_array(ii) > order_array(ii+1_ik)) then
        temp=order_array(ii)
        temp_int=array_int(ii)
        temp_real=array_real(ii)
        order_array(ii) = order_array(ii+1_ik)
        array_int(ii) = array_int(ii+1_ik) 
        array_real(ii) = array_real(ii+1_ik)
        order_array(ii+1_ik) = temp
        array_int(ii+1_ik) = temp_int
        array_real(ii+1_ik) = temp_real
        swapped = .TRUE._lk
      end if
    end do
    if(.NOT. swapped) exit
  end do
end subroutine cel_sp_mat_bubble_sort

recursive subroutine cel_sp_mat_merge_sort(order_array,array_int,array_real,&
 local_length, r_level, output_on)
 
  integer(kind=ik), dimension(:),  intent(inout) :: order_array
  integer(kind=ik), dimension(:),  intent(inout) :: array_int
  real(kind=rk), dimension(:),  intent(inout) :: array_real
  integer(kind=ik),intent(in)  :: local_length
  integer(kind=ik), dimension(:), allocatable:: int_1_tmp
  integer(kind=ik), dimension(:), allocatable:: int_2_tmp
  real(kind=rk), dimension(:), allocatable :: real_tmp
  integer(kind=ik),  intent(in) :: r_level
  logical (kind=lk), intent(in) :: output_on
  !real(kind=rk) :: swap_value_real
  !integer(kind=ik) :: swap_value_int,swap_value_int_2
  integer(kind=ik) :: left_part, right_part, err_code
  integer(kind=ik) :: next_r_level

  if (local_length<2) return
  
  if((r_level .GT. cel_algorithms_param_merge_max_r_level / 2) .OR. &
    local_length .LT. cel_algorithms_param_merge_max_array_length) then
    ! write(*,'(A,I0)') "--call left cel_sp_mat_bubble_sort left_part:", left_part
    call cel_sp_mat_bubble_sort(order_array(1:local_length),&
      array_int(1:local_length),array_real(1:local_length),&
      local_length)
    return
  endif
  !  if(local_length == 2) then
  !    if (order_array(1) > order_array(2)) then
  !       swap_value_int = order_array(1)
  !       order_array(1) = order_array(2)
  !       order_array(2) = swap_value_int
  !       swap_value_int_2 = array_int(1)
  !       array_int(1) = array_int(2)
  !       array_int(2) = swap_value_int_2
  !       swap_value_real = array_real(1)
  !       array_real(1) = array_real(2)
  !       array_real(2) = swap_value_real
  !    endif
  !   return
  !endif
  !write(*,'(A,I0)') "recursion level:", r_level
  left_part=(local_length+1)/2
  right_part=local_length-left_part
  next_r_level = r_level + 1
  
    !write(*,'(A,I0)') "--call left cel_sp_mat_merge_sort left_part:", left_part
  call cel_sp_mat_merge_sort(order_array(1:left_part),&
      array_int(1:left_part),array_real(1:left_part),&
      left_part, next_r_level, output_on)

   call cel_sp_mat_merge_sort(order_array(left_part+1_ik:local_length),&
      array_int(left_part+1_ik:local_length),array_real(left_part+1_ik:local_length),&
      right_part, next_r_level, output_on)


   if (order_array(left_part) > order_array(left_part+1)) then
      if(allocated(int_1_tmp))deallocate(int_1_tmp)
      allocate(int_1_tmp(left_part), stat=err_code)
      call cel_error("cel_sp_mat_merge_sort int_1_tmp", err_code, .TRUE._lk, .TRUE._lk)
      if(allocated(int_2_tmp))deallocate(int_2_tmp)
      allocate(int_2_tmp(left_part), stat=err_code)
      call cel_error("cel_sp_mat_merge_sort int_2_tmp", err_code, .TRUE._lk, .TRUE._lk)
      if(allocated(real_tmp))deallocate(real_tmp)
      allocate(real_tmp(left_part), stat=err_code)
      call cel_error("cel_sp_mat_merge_sort real_tmp", err_code, .TRUE._lk, .TRUE._lk)
      int_1_tmp(1:left_part)=order_array(1:left_part)
      int_2_tmp(1:left_part)=array_int(1:left_part)
      real_tmp(1:left_part)=array_real(1:left_part)
      call cel_sp_mat_merge(int_1_tmp,int_2_tmp,real_tmp,left_part,&
        order_array(left_part+1:local_length),&
        array_int(left_part+1:local_length),&
        array_real(left_part+1:local_length),&
        right_part,order_array,array_int,array_real,local_length)
   else
      allocate(int_1_tmp(0_ik))
      allocate(int_2_tmp(0_ik))
      allocate(real_tmp(0_ik))
   endif
   return
 
end subroutine cel_sp_mat_merge_sort




subroutine cel_sp_mat_merge(order_array_a,array_a_int,array_a_real,&
        length_a,&
        order_array_b,array_b_int,array_b_real,&
        length_b,&
        order_array,array_int,array_real,local_length)
  integer(kind=ik), dimension(:), intent(inout) :: order_array_a
  integer(kind=ik), dimension(:), intent(inout) :: array_a_int
  real(kind=rk), dimension(:), intent(inout) :: array_a_real
  integer(kind=ik), intent(in) ::length_a
  integer(kind=ik), dimension(:), intent(in) :: order_array_b
  integer(kind=ik), dimension(:), intent(in) :: array_b_int
  real(kind=rk), dimension(:), intent(in) :: array_b_real
  integer(kind=ik), intent(in) ::length_b
  integer(kind=ik), dimension(:), intent(inout) :: order_array
  integer(kind=ik), dimension(:), intent(inout) :: array_int
  real(kind=rk), dimension(:),  intent(inout) :: array_real
  integer(kind=ik), intent(in) ::local_length
  integer(kind=ik) :: ii,jj,kk
 
   ii = 1_ik; jj = 1_ik; kk = 1_ik;
   do while(ii <= length_a .and. jj <= length_b)
      if (order_array_a(ii) <= order_array_b(jj)) then
         order_array(kk) = order_array_a(ii)
         array_int(kk) = array_a_int(ii)
         array_real(kk) = array_a_real(ii)
         ii = ii+1_ik
      else
         order_array(kk) = order_array_b(jj)
         array_int(kk) = array_b_int(jj)
         array_real(kk) = array_b_real(jj)
         jj = jj+1_ik
      endif
      kk = kk + 1_ik
   enddo
   do while (ii <= length_a)
      order_array(kk) = order_array_a(ii)
      array_int(kk) = array_a_int(ii)
      array_real(kk) = array_a_real(ii)
      ii = ii + 1_ik
      kk = kk + 1_ik
   enddo
   return
 
end subroutine cel_sp_mat_merge
 
end module cel_sp_mat_module
