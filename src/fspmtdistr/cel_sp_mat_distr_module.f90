!     
! File:  cel_blspmt_cel.f90
! Project CRESTA (see details on https://cresta-project.eu) Exascale library
! Send you email about the bugs to
! The High Performance Computing Center Stuttgart (HLRS) of the University of Stuttgart
! Dmitry Khabi,  (email: khabi@hlrs.de)
! Created on Aug 27, 2013
!

module cel_sp_mat_distr_module
!library fbase
use cel_types_module
use cel_input_parameters_module
use cel_omp_module
use cel_timer_interface_module
!use cel_cg_module
!library fspmt
use cel_vec_module
use cel_sp_mat_module
!library fblspmt

!external libraries
use MPI
implicit none

contains

!subdivide the local part of the sparse matrix (sp_mat_local) in to the local sub matrices (or blocks) and save it in  sp_mat_local_array
!out sp_mat_local_array(i)  the local part of thesub matrix (output, coo format)
!procs_array_first_idx(i),procs_array_last_idx(i) - indices of the first and last columns with nonzero elements in the submatrix(i)
!in procs_array(i) id of the proc referenced by the columns in submatrix sp_mat_local_array(i)
!in vtxdist - distribution of the rows of the sparse matrix
!in sp_mat_local - the local part of the sparse matrix (input, coo format)
subroutine cel_sp_mat_distr_from_sp_mat_local(sp_mat_local_array, procs_array,&
  procs_array_first_idx,procs_array_last_idx,&
  com_buffer_perm,length_comm_buffer,& 
  vtxdist, sp_mat_local, omp, input_parameters,err_code) 
  type(cel_sp_mat_type), dimension(:), allocatable, intent(out) :: sp_mat_local_array
  integer(kind=ik), dimension(:), allocatable , intent(out)     :: procs_array ! array with proc_id for the communication elements in sp_mat_local_array
  integer(kind=ik), dimension(:), allocatable , intent(out)     :: procs_array_first_idx ! array with first idx of data needed from the process in procs_array
  integer(kind=ik), dimension(:), allocatable , intent(out)     :: procs_array_last_idx ! array with last idx of data needed from the process in procs_array
  integer(kind=ik), dimension(:,:), allocatable , intent(out)   :: com_buffer_perm
  integer(kind=ik), dimension(:), allocatable , intent(out)     :: length_comm_buffer
  integer(kind=ik), dimension(:), allocatable, intent(in)       :: vtxdist ! array of matrix distribution (i-th el. is the index of the first row on i-th process)
  type(cel_sp_mat_type), intent(in)                             :: sp_mat_local
  type(cel_omp_type), intent(in)                                :: omp
  type(cel_input_parameters_type), intent(in)                   :: input_parameters
  integer (kind=ik), intent(out)                                :: err_code
  integer (kind=ik), dimension(:), allocatable                  :: num_elements_per_proc 
  integer (kind=ik), dimension(:), allocatable                  :: num_elements_per_block 
  type(cel_sp_mat_type) :: smpt_coo_tmp
  integer(kind=ik) :: ii, num_blocks, idx,idx2, row, col, last_col, proc_id, last_proc_id, jj
  integer(kind=ik) :: coo_local_array_idx, cur_procs_array_idx, comm_buffer_num, length
  integer(kind=ik) :: block_offset, el_in_comm_buffer, sum_max
  real(kind=rk) :: el_value, diagonal_value
  integer(kind=ik), dimension(:), allocatable  :: cur_local_coo_idx, inv_procs_array
  integer(kind=ik), dimension(:,:), allocatable  :: vec_comm_voll_buffer 
  integer(kind=ik) :: sys_stat,row_local_idx
  integer(kind=ik) :: length_nieghbours, my_row_offset
  err_code=0_ik

  
  call cel_sp_mat_distr_get_num_neighbours(num_elements_per_proc,sp_mat_local, vtxdist, omp, err_code)
  call cel_error("Error in cel_sp_mat_distr_from_coo_global: cel_sp_mat_distr_get_num_neighbours", &
        &err_code,cel_is_in_debug,.TRUE._lk)
  !divide the local matrix sp_mat_local in the several sub matrices:
  !the first part consists of the elements, whose column indexes belong to the local part of the vector x (part of Ax  without communication)
  !the i-th part consists of the elements, whose column indexes belong to the part of the vector x on the proc(i)-th process  (part of Ax with communication)
  if(omp%is_master) then
    num_blocks = 0_ik
    num_blocks = 0_ik
    sum_max=0_ik
    
    !calculate size of the new matrix
    do ii=1, size(num_elements_per_proc)
      if(num_elements_per_proc(ii) .gt. 0_ik) num_blocks = num_blocks + 1_ik
    end do
    
    !if not empty sub matrix
    if(num_blocks .GT. 0) then
      !allocate memory and reset that for the new matrix
      length_nieghbours = 0_ik
      do ii=1,size(vtxdist)-1_ik
         length = vtxdist(ii+1_ik)-vtxdist(ii)
         length_nieghbours = max(length,length_nieghbours)
      end do
      length = vtxdist(omp%master_num+2_ik)-vtxdist(omp%master_num+1_ik)
      
      allocate(sp_mat_local_array(num_blocks), stat=sys_stat)       
        call cel_error("Error cel_sp_mat_distr_from_sp_mat_local alloc 1", sys_stat,&
                        cel_is_in_debug,.TRUE._lk)
      allocate(procs_array(num_blocks), stat=sys_stat)
        call cel_error("Error cel_sp_mat_distr_from_sp_mat_local alloc 3", sys_stat,&
                        cel_is_in_debug,.TRUE._lk)
      allocate(procs_array_first_idx(num_blocks), stat=sys_stat)
        call cel_error("Error cel_sp_mat_distr_from_sp_mat_local alloc 4", sys_stat,&
                        cel_is_in_debug,.TRUE._lk)
      allocate(procs_array_last_idx(num_blocks), stat=sys_stat)
        call cel_error("Error cel_sp_mat_distr_from_sp_mat_local alloc 5", sys_stat,&
                        cel_is_in_debug,.TRUE._lk)
      allocate(num_elements_per_block(num_blocks), stat=sys_stat)
        call cel_error("Error cel_sp_mat_distr_from_sp_mat_local alloc 6", sys_stat,&
                        cel_is_in_debug,.TRUE._lk)
      allocate(length_comm_buffer(num_blocks), stat=sys_stat)
        call cel_error("Error cel_sp_mat_distr_from_sp_mat_local alloc 7", sys_stat,&
                        cel_is_in_debug,.TRUE._lk)
      allocate(vec_comm_voll_buffer(length_nieghbours,num_blocks), stat=sys_stat)
        call cel_error("Error cel_sp_mat_distr_from_sp_mat_local alloc 8", sys_stat,&
                        cel_is_in_debug,.TRUE._lk)
      allocate(inv_procs_array(omp%num_masters), stat=sys_stat)
        call cel_error("Error cel_sp_mat_distr_from_sp_mat_local alloc 9", sys_stat,&
                        cel_is_in_debug,.TRUE._lk)                  
                        
      procs_array_first_idx=cel_index_max
      procs_array_last_idx=-1_ik
      num_elements_per_block=0_ik
      length_comm_buffer=0_ik
      vec_comm_voll_buffer=0_ik
      
      idx = 1_ik
      do ii=1, size(num_elements_per_proc)
        if(num_elements_per_proc(ii) .gt. 0_ik) then
          procs_array(idx) = ii-1_ik
          idx = idx + 1_ik
        end if
      end do
       !swap two local sub matrices, because the first must be the block with 
       !local (no mpi) communication 
       !this is a diagonal sub matrix with local rows and columns indexes)
       !one per mpi process
      if(size(procs_array) .gt. 1_ik) then
        if(procs_array(1_ik) .ne. omp%master_num) then
           do ii=1, size(procs_array)
            if(omp%master_num .eq. procs_array(ii)) then
                idx = ii
                exit
            end if
           end do
           do ii=idx, 2_ik, -1_ik
              procs_array(ii) = procs_array(ii-1_ik)
           end do
           procs_array(1_ik) = omp%master_num
         end if
      end if
      do ii=1, size(procs_array)
        inv_procs_array(procs_array(ii)+1_ik)=ii
      end do
      idx2 = 1_ik
      do ii=1, size(procs_array)
          !allocate the memory for a local sub matrix sp_mat_local_array(idx2)
          idx = procs_array(ii)+1_ik
          call new(sp_mat_local_array(idx2),&
            sp_mat_local%rows_num,&
            sp_mat_local%cols_num,&
            sp_mat_local%block_rows_num,&
            sp_mat_local%block_cols_num,&
            num_elements_per_proc(idx))
          num_elements_per_block(idx2)=num_elements_per_proc(idx)
          idx2 = idx2  + 1_ik
      end do
      allocate(sp_mat_local_array(1_ik)%diagonal(length), stat=sys_stat)
      call cel_error("Error cel_sp_mat_distr_from_sp_mat_local  alloc 2", sys_stat,&
         cel_is_in_debug,.TRUE._lk)
    else
      err_code = 1_ik
      call cel_error("Error in cel_sp_mat_distr_from_coo_global: matrix"//&
      " is empty, expected not empty! 10", &
          &err_code,cel_is_in_debug,.TRUE._lk)
    endif
    allocate(cur_local_coo_idx(num_blocks), stat=sys_stat)
        call cel_error("Error cel_sp_mat_distr_from_sp_mat_local" //&
        " cur_local_coo_idx alloc 11", sys_stat,&
         cel_is_in_debug,.TRUE._lk)
    cur_local_coo_idx = 1_ik
    last_proc_id = 0_ik
    last_col = -1_ik
    !distribute the all matrix elements between the local sub matrices
    my_row_offset = vtxdist(omp%master_num+1_ik) - 1_ik
    do ii=1, sp_mat_local%elements_num
      col = sp_mat_local%cols(ii)
      row = sp_mat_local%rows(ii)
      el_value = sp_mat_local%values(ii)
      if(last_col .ne. col) then
        call cel_sp_mat_distr_get_proc_id(proc_id, last_proc_id, omp%num_masters, col, vtxdist)
        last_col = col
      end if
      if (proc_id .NE. cel_index_undef) then
        cur_procs_array_idx = cel_index_undef
        cur_procs_array_idx = inv_procs_array(proc_id+1_ik)
        !do jj=1, num_blocks
        !  if(procs_array(jj) .eq. proc_id) then
        !    cur_procs_array_idx = jj
        !    exit
        !  end if
        !end do
          !if diagonal element, save it
          if(col .eq. row) then
            if(input_parameters%matrix_check_null_diagonal) then
              if(abs(el_value) .le. cel_real_small) then
                call cel_error("Error in cel_sp_mat_distr_from_coo_global:"//&
                " matrix diagonal is zero! 12", &
                   row,cel_is_in_debug,.TRUE._lk)
              end if
            end if
            row_local_idx = row - my_row_offset
            sp_mat_local_array(1_ik)%diagonal(row_local_idx)=el_value
          end if
      else
        err_code = ii
        call cel_error("Error in cel_sp_mat_distr_from_coo_global: the col index is wrong! 13", &
          err_code,cel_is_in_debug,.TRUE._lk)
        return
      end if
      !update the local sub matrix, which owns the matrix element
      if(cur_procs_array_idx .NE. cel_index_undef) then
        idx = cur_local_coo_idx(cur_procs_array_idx)
        sp_mat_local_array(cur_procs_array_idx)%rows(idx) = row
        sp_mat_local_array(cur_procs_array_idx)%cols(idx) = col
        sp_mat_local_array(cur_procs_array_idx)%values(idx) = el_value
        cur_local_coo_idx(cur_procs_array_idx) = cur_local_coo_idx(cur_procs_array_idx) + 1_ik
        if(procs_array_first_idx(cur_procs_array_idx) .gt. col) then
          procs_array_first_idx(cur_procs_array_idx) = col
        end if
        if(procs_array_last_idx(cur_procs_array_idx) .lt. col) then
          procs_array_last_idx(cur_procs_array_idx) = col
        end if
      
        block_offset = vtxdist(procs_array(cur_procs_array_idx)+1_ik)-1_ik!(procs_array(cur_procs_array_idx))*length
      !  if(omp%master_num .eq. 2) then
      !    write(*,"(5(A,I0))") "proc:",omp%master_num,&
      !       "; procs_array(cur_procs_array_idx)", procs_array(cur_procs_array_idx),&
      !      "; col:",col,&
      !      "; cur_procs_array_idx:",cur_procs_array_idx,&
      !      "; block offset:",block_offset
      !  end if
        vec_comm_voll_buffer(col-block_offset,cur_procs_array_idx) = 1_ik
      else
        err_code = 1
        call cel_error("Error in cel_sp_mat_distr_from_coo_global: the proc index is wrong! 14", &
          err_code,cel_is_in_debug,.TRUE._lk)
      end if
    end do
    !correction of the new matrix incl. new formating
    comm_buffer_num = 0_ik
    idx=1_ik
    do ii=1_ik,num_blocks
      sp_mat_local_array(ii)%row_offset = vtxdist(omp%master_num+1_ik) - 1_ik
      sp_mat_local_array(ii)%col_offset = vtxdist(procs_array(ii)+1_ik) - 1_ik
      call cel_sp_mat_to_csr(sp_mat_local_array(ii))
      sp_mat_local_array(ii)%is_in_csr = .true._lk
      sp_mat_local_array(ii)%is_in_coo = .true._lk
      procs_array_last_idx(ii)=procs_array_last_idx(ii)-sp_mat_local_array(ii)%col_offset 
      procs_array_first_idx(ii)= procs_array_first_idx(ii)-sp_mat_local_array(ii)%col_offset 
      !check whether it is better to use a buffer for the communication - under the development
      el_in_comm_buffer = SUM(vec_comm_voll_buffer(1_ik:length_nieghbours,ii))
!     if(omp%master_num .eq. 0) then
!       write(*,"(2(A,I0))") "el_in_comm_buffer:",el_in_comm_buffer,"; length;",length
!      end if
      if(ii .gt. 1_ik) sum_max=max(sum_max,el_in_comm_buffer)
      comm_buffer_num = comm_buffer_num +1_ik
      length_comm_buffer(ii)=el_in_comm_buffer
    end do
    if(sum_max .gt. 0) then
      allocate(com_buffer_perm(sum_max,num_blocks), stat=sys_stat)
        call cel_error("Error cel_sp_mat_distr_from_sp_mat_local alloc ", sys_stat,&
                        cel_is_in_debug,.TRUE._lk)
      !com_buffer_perm = 0_ik
      !allocate communication buffer for the matrix vector multiplication
      idx=1_ik
      do ii=2_ik,num_blocks
          idx2=1_ik
          idx = procs_array(ii)+1_ik
          length = vtxdist(idx+1_ik)-vtxdist(idx)
          do jj=1_ik,length
            if(vec_comm_voll_buffer(jj,ii) .gt. 0_ik) then
              com_buffer_perm(idx2,ii)=jj
              idx2=idx2+1_ik
            end if
          end do
      end do
    end if

   ! if(omp%master_num .eq. 0) then
   !   do ii=2_ik,num_blocks
   !     call print("vec_comm_voll_buffer",vec_comm_voll_buffer(1:length,ii))
   !   end do
   !   call print("length_comm_buffer",length_comm_buffer)
   !   do ii=2_ik,num_blocks
   !     call print("com_buffer_perm",com_buffer_perm(:,ii))
   !   end do
   ! end if
    deallocate(inv_procs_array)
  else
    allocate(num_elements_per_proc (0))
    allocate(num_elements_per_block(0))
    allocate(cur_local_coo_idx(0))
    allocate(vec_comm_voll_buffer(0,0))
  end if
  !$omp barrier

end subroutine cel_sp_mat_distr_from_sp_mat_local
!calculate vtxdist
!initialize sp_mat_local with local part of the matrix
subroutine cel_sp_mat_distr_local(sp_mat_local, vtxdist, sp_mat_global,&
  omp, err_code)
  type(cel_sp_mat_type), intent(inout)  :: sp_mat_local
  integer(kind=ik), dimension(:), allocatable, intent(inout)   :: vtxdist ! array of matrix distribution (i-th el. is the index of the first row on i-th process)
  type(cel_sp_mat_type), intent(in)  :: sp_mat_global
  type(cel_omp_type), intent(in) :: omp
  integer (kind=ik), intent(out)          :: err_code
  integer (kind=ik) :: num_blocks, num_blocks_one, num_blocks_last, last_row
  integer (kind=ik) :: ii
  integer (kind=ik), dimension(:), allocatable :: num_elements
  integer (kind=ik), dimension(:), allocatable :: distr_index
  integer (kind=mpik) :: ierr, count_mpik, dest_mpik, source_mpik, tag_mpik
  integer (kind=ik)   :: row, last_proc_id, proc_id, sum_el, start_idx, end_idx
  integer(kind=mpik),allocatable, dimension(:)  ::send_request, send_counts, displs
  integer(kind=mpik),allocatable, dimension(:) :: get_request
  integer(kind=ik) :: sys_stat
  integer(kind=mpik) :: get_request_1,get_request_2,get_request_3
  integer(kind=mpik) :: status_1(MPI_STATUS_SIZE),status_2(MPI_STATUS_SIZE),status_3(MPI_STATUS_SIZE)
  integer (kind=ik), dimension(0_ik) :: dummy_int
  integer (kind=ik), dimension(0_ik) :: dummy_real
  err_code=0
  if(omp%is_master) then
    allocate(vtxdist(omp%num_masters+1_ik), stat=sys_stat)
    call cel_error("Error cel_sp_mat_distr_local vtxdist alloc ", sys_stat,&
                        cel_is_in_debug,.TRUE._lk)

    num_blocks = sp_mat_global%rows_num / sp_mat_global%block_rows_num
    num_blocks_one = num_blocks / omp%num_masters
    num_blocks_last = num_blocks_one + ( num_blocks - num_blocks_one*omp%num_masters )
    vtxdist = num_blocks_one
    vtxdist(1)=1_ik
    vtxdist(omp%num_masters+1_ik)=num_blocks_last

    do ii=1,omp%num_masters
      vtxdist(ii+1_ik)=vtxdist(ii+1_ik)+vtxdist(ii)
    end do

    sp_mat_local%rows_num = vtxdist(omp%master_num+2)-vtxdist(omp%master_num+1)
    sp_mat_local%cols_num = sp_mat_global%cols_num
    sp_mat_local%block_rows_num = sp_mat_global%block_rows_num
    sp_mat_local%block_cols_num = sp_mat_global%block_cols_num
    !calculate number of elements per process
    allocate(num_elements(omp%num_masters))
    allocate(distr_index(omp%num_masters))

    if(omp%master_num .EQ. 0) then
      allocate(send_request(3*omp%num_masters))
      num_elements = 0_ik
      sum_el = 0_ik
      last_proc_id = 0_ik
      distr_index = 0_ik
      distr_index(1) = 1_ik
      last_row = -1_ik
      do ii=1_ik, sp_mat_global%elements_num
        row = sp_mat_global%rows(ii)
        if(last_row .ne. row) then
          call cel_sp_mat_distr_get_proc_id(proc_id, last_proc_id, omp%num_masters, row, vtxdist)
        end if
        if (proc_id .NE. cel_index_undef) then
          num_elements(proc_id+1_ik) = num_elements(proc_id+1_ik) + 1_ik
          if(last_proc_id .NE. proc_id) then
            if(distr_index(proc_id+1_ik) .ne. 0_ik .and. &
              proc_id .gt. 0_ik) then
              call cel_error("Error in cel_sp_mat_distr_distribute: matrix data is not sorted!", &
               &ii,cel_is_in_debug,.TRUE._lk)
            end if
            distr_index(proc_id+1_ik) = ii
            last_proc_id = proc_id
          end if
        else
          call cel_error("Error in cel_sp_mat_distr_distribute: the row index is wrong 1!", &
          &row,cel_is_in_debug,.TRUE._lk)
          err_code = row
          return
        end if
      end do
    end if
    !distribute the number of elemenets per process
    count_mpik = int(omp%num_masters,kind=mpik)
    call MPI_Bcast(num_elements(1_ik:omp%num_masters),count_mpik,omp%mpi_type_ik,&
                   0_mpik,omp%master_comm,ierr)
    !allocate memory for local part of the matix
    call new(sp_mat_local, sp_mat_local%rows_num, sp_mat_local%cols_num,&
      sp_mat_local%block_rows_num, sp_mat_local%block_cols_num ,num_elements(omp%master_num + 1_ik))
    !distribute the matrix
    allocate(send_counts(omp%num_masters))
    allocate(displs(omp%num_masters))
    do ii=1_ik, omp%num_masters
      send_counts(ii) = INT(num_elements(ii),kind=mpik)
    end do
    displs(1_ik)=0_mpik
    do ii=2_ik, omp%num_masters
      displs(ii)=displs(ii-1_ik)+send_counts(ii-1_ik)
    end do
    count_mpik = int(sp_mat_global%elements_num,kind=mpik)
    count_mpik = int(sp_mat_global%elements_num,kind=mpik)
    start_idx = 1_ik
    end_idx = num_elements(omp%master_num + 1_ik)
    if(omp%master_num .EQ. 0) then
      call MPI_Scatterv(sp_mat_global%rows(1_ik:sp_mat_global%elements_num), &
        send_counts(1_ik:omp%num_masters), displs(1_ik:omp%num_masters),omp%mpi_type_ik, &
        sp_mat_local%rows(start_idx:end_idx), send_counts(omp%master_num+1_ik),&
        omp%mpi_type_ik, 0_mpik,omp%master_comm,ierr)
      
      call MPI_Scatterv(sp_mat_global%cols(1_ik:sp_mat_global%elements_num), &
        send_counts(1_ik:omp%num_masters), displs(1_ik:omp%num_masters),omp%mpi_type_ik, &
        sp_mat_local%cols(start_idx:end_idx), send_counts(omp%master_num+1_ik),&
        omp%mpi_type_ik, 0_mpik,omp%master_comm,ierr)
        
      call MPI_Scatterv(sp_mat_global%values(1_ik:sp_mat_global%elements_num), &
        send_counts(1_ik:omp%num_masters), displs(1_ik:omp%num_masters),omp%mpi_type_rk, &
        sp_mat_local%values(start_idx:end_idx), send_counts(omp%master_num+1_ik),&
        omp%mpi_type_rk, 0_mpik,omp%master_comm,ierr)        
    else
     call MPI_Scatterv(dummy_int, &
        send_counts(1_ik:omp%num_masters), displs(1_ik:omp%num_masters),omp%mpi_type_ik, &
        sp_mat_local%rows(start_idx:end_idx), send_counts(omp%master_num+1_ik),&
        omp%mpi_type_ik, 0_mpik,omp%master_comm,ierr)
        
      call MPI_Scatterv(dummy_int, &
        send_counts(1_ik:omp%num_masters), displs(1_ik:omp%num_masters),omp%mpi_type_ik, &
        sp_mat_local%cols(start_idx:end_idx), send_counts(omp%master_num+1_ik),&
        omp%mpi_type_ik, 0_mpik,omp%master_comm,ierr)
        
      call MPI_Scatterv(dummy_real, &
        send_counts(1_ik:omp%num_masters), displs(1_ik:omp%num_masters),omp%mpi_type_rk, &
        sp_mat_local%values(start_idx:end_idx), send_counts(omp%master_num+1_ik),&
        omp%mpi_type_rk, 0_mpik,omp%master_comm,ierr)        
    end if
    if(.false.) then
      if(omp%master_num .EQ. 0) then
        do ii=1_ik, omp%num_masters
          start_idx = distr_index(ii)
          end_idx = start_idx + num_elements(ii) - 1_ik
          count_mpik = INT(num_elements(ii),kind=mpik)
          dest_mpik = INT(ii-1_ik,kind=mpik)
          tag_mpik = dest_mpik * 3_mpik
          call MPI_Isend(sp_mat_global%rows(start_idx:end_idx),count_mpik,omp%mpi_type_ik,&
                dest_mpik,tag_mpik,omp%master_comm,send_request((ii-1_ik)*3_ik+1_ik),ierr)
          call cel_error_proc(&
                 "Error cel_sp_mat_distr_local: MPI_Isend 1!", &
                  int(omp%master_num,kind=mpik), int(ierr,kind=mpik), cel_is_in_debug, .FALSE._lk)
          tag_mpik = dest_mpik * 5_mpik
          call MPI_Isend(sp_mat_global%cols(start_idx:end_idx),count_mpik,omp%mpi_type_ik,&
                dest_mpik,tag_mpik,omp%master_comm,send_request((ii-1_ik)*3_ik+2_ik),ierr)
          call cel_error_proc(&
                 "Error cel_sp_mat_distr_local: MPI_Isend 2!", &
                  int(omp%master_num,kind=mpik), int(ierr,kind=mpik), cel_is_in_debug, .FALSE._lk)
          tag_mpik = dest_mpik * 7_mpik
          call MPI_Isend(sp_mat_global%values(start_idx:end_idx),count_mpik,omp%mpi_type_rk,&
                dest_mpik,tag_mpik,omp%master_comm,send_request((ii-1_ik)*3_ik+3_ik),ierr)
          call cel_error_proc(&
                 "Error cel_sp_mat_distr_local: MPI_Isend 3!", &
                  int(omp%master_num,kind=mpik), int(ierr,kind=mpik), cel_is_in_debug, .FALSE._lk)
        end do
      end if
      allocate(get_request(4_ik))
      start_idx = 1_ik
      end_idx = num_elements(omp%master_num + 1_ik)
      count_mpik = INT(num_elements(omp%master_num+1_ik),kind=mpik)
      source_mpik = 0_mpik
      tag_mpik = int(omp%master_num,kind=mpik) * 3_mpik
      call MPI_Irecv(sp_mat_local%rows(start_idx:end_idx),count_mpik,omp%mpi_type_ik,&
            source_mpik,tag_mpik,omp%master_comm,get_request_1,ierr)
      call cel_error_proc(&
                 "Error cel_sp_mat_distr_local: MPI_Irecv 1!", &
                  int(omp%master_num,kind=mpik), int(ierr,kind=mpik), cel_is_in_debug, .FALSE._lk)
      tag_mpik = int(omp%master_num,kind=mpik) * 5_mpik
      call MPI_Irecv(sp_mat_local%cols(start_idx:end_idx),count_mpik,omp%mpi_type_ik,&
            source_mpik,tag_mpik,omp%master_comm,get_request_2,ierr)
      call cel_error_proc(&
                 "Error cel_sp_mat_distr_local: MPI_Irecv 2!", &
                  int(omp%master_num,kind=mpik), int(ierr,kind=mpik), cel_is_in_debug, .FALSE._lk)
      tag_mpik = int(omp%master_num,kind=mpik) * 7_mpik
      call MPI_Irecv(sp_mat_local%values(start_idx:end_idx),count_mpik,omp%mpi_type_rk,&
            source_mpik,tag_mpik,omp%master_comm,get_request_3,ierr)
      call cel_error_proc(&
                 "Error cel_sp_mat_distr_local: MPI_Irecv 3!", &
                  int(omp%master_num,kind=mpik), int(ierr,kind=mpik), cel_is_in_debug, .FALSE._lk)
      sp_mat_local%block_indx = omp%master_num
      call MPI_Wait(get_request_1,status_1,ierr)
      call MPI_Wait(get_request_2,status_2,ierr)
      call MPI_Wait(get_request_3,status_3,ierr)
      call cel_error_proc(&
                 "Error cel_sp_mat_distr_local: MPI_Waitall!", &
                  int(omp%master_num,kind=mpik), int(ierr,kind=mpik), cel_is_in_debug, .FALSE._lk)
    end if
  end if
  !$omp barrier

end subroutine cel_sp_mat_distr_local

!calculate how many elements of the local part of the matrix do belong to
!the local sub matrix, whose columns make it necessary to communicate in order to perform the matrix vector multiplication Ax
subroutine cel_sp_mat_distr_get_num_neighbours(num_elements, sp_mat_local, vtxdist, omp, err_code)
  integer (kind=ik), dimension(:), allocatable,intent(out) :: num_elements
  type(cel_sp_mat_type), intent(in)  :: sp_mat_local
  integer(kind=ik), dimension(:), allocatable, intent(in)   :: vtxdist ! array of matrix distribution (i-th el. is the index of the first row on i-th process)
  type(cel_omp_type), intent(in) :: omp
  integer (kind=ik), intent(out)          :: err_code
  integer (kind=ik)          :: sum_el, last_proc_id, col, proc_id
  integer (kind=ik) :: size_num_elements, ii, print_proc
  integer (kind=mpik) :: ierr
  err_code = 0_ik
  if(omp%is_master) then
    allocate(num_elements(omp%num_masters))
    num_elements = 0_ik
    sum_el = 0_ik
    last_proc_id = 0_ik
    do ii=1_ik, sp_mat_local%elements_num
      col = sp_mat_local%cols(ii)
      call cel_sp_mat_distr_get_proc_id(proc_id, last_proc_id, omp%num_masters, col, vtxdist)
      if (proc_id .NE. cel_index_undef) then
        num_elements(proc_id+1_ik) = num_elements(proc_id+1_ik) + 1_ik
        if(last_proc_id .NE. proc_id) then
          last_proc_id = proc_id
        end if
      else
        !$omp barrier
        call cel_error("Error cel_sp_mat_distr_get_num_neighbours the row index is wrong 2!", &
        &ii,cel_is_in_debug,.TRUE._lk)
        err_code = 1_ik
        EXIT
      end if
    end do
   end if
  !$omp barrier
  !if(omp%is_master) then
  ! print_proc = 0_ik
  ! do ii=1_ik, omp%num_masters
  !  if(print_proc .eq. omp%master_num) then
  !    write(cel_output_unit,'(A,I0,A,I0)') "(",omp%master_num,"): num_elements:",num_elements
  !  end if
  !  print_proc = print_proc + 1_ik
  !  call MPI_BARRIER(omp%master_comm, ierr)
  ! end do
  !end if
  !!$omp barrier
end subroutine cel_sp_mat_distr_get_num_neighbours

!get proc_id (mpi_rank) to which the matrix row belongs
subroutine cel_sp_mat_distr_get_proc_id(proc_id, last_proc_id, num_procs, row, vtxdist)
  integer(kind=ik), intent(out)               :: proc_id
  integer(kind=ik), intent(in)                :: last_proc_id
  integer(kind=ik), intent(in)                :: row
  integer(kind=ik), dimension(:),allocatable, intent(in)  :: vtxdist
  integer(kind=ik)                            :: num_procs
  integer(kind=ik)                            :: ii, proc_index
  integer(kind=ik)                            :: low_bound, upper_bound
  
  proc_index = last_proc_id + 1_ik
  proc_id = cel_index_undef
  upper_bound = 0_ik
  low_bound = 0_ik
  
  if(proc_index .LE. size(vtxdist)-1) then
    low_bound = vtxdist(proc_index)
    upper_bound = vtxdist(proc_index+1_ik)
  else
    call cel_error("Error cel_sp_mat_distr_get_proc_id 0 proc_index is wrong!", &
        proc_index,cel_is_in_debug,.TRUE._lk)
  end if
  if( ( low_bound .LE. row ) .AND. ( upper_bound .GT. row ) ) then
    proc_id = last_proc_id
  else 
    if(upper_bound .LE. row) then
      do ii=last_proc_id+1_ik, num_procs - 1_ik
        proc_index = ii + 1_ik
        if(proc_index .gt. size(vtxdist)-1_ik) then
          call cel_error("Error cel_sp_mat_distr_get_proc_id 1 proc_index is wrong!", &
        proc_index,cel_is_in_debug,.TRUE._lk)
        end if
        low_bound = vtxdist(proc_index)
        upper_bound = vtxdist(proc_index+1_ik)
        if( ( low_bound .LE. row ) .AND. ( upper_bound .GT. row) ) then
          proc_id = ii
          EXIT
        end if
      end do
    else
      do ii=last_proc_id, 0_ik, -1_ik
        proc_index = ii + 1_ik
        if(proc_index .gt. size(vtxdist)-1_ik) then
          call cel_error("Error cel_sp_mat_distr_get_proc_id 2 proc_index is wrong!", &
        proc_index,cel_is_in_debug,.TRUE._lk)
        end if
        low_bound = vtxdist(proc_index)
        upper_bound = vtxdist(proc_index+1_ik)
        if( ( low_bound .LE. row ) .AND. ( upper_bound .GT. row) ) then
          proc_id = ii
          EXIT
        end if
      end do
    end if
  end if
end subroutine cel_sp_mat_distr_get_proc_id


end module cel_sp_mat_distr_module
