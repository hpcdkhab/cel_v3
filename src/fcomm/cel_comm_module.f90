!     
! File:   cel_blsp_mat_connection_module.f90
! Project CRESTA (see details on https://cresta-project.eu) Exascale library
! Send you email about the bugs to
! The High Performance Computing Center Stuttgart (HLRS) of the University of Stuttgart
! Dmitry Khabi,  (email: khabi@hlrs.de)
! Created on Jun 14, 2013
!

module cel_comm_module
use cel_types_module
use cel_timer_interface_module
use cel_cpu_param_module
use cel_perf_module
use cel_base_print_module
use cel_omp_module
use MPI
use OMP_LIB
implicit none

interface del
  procedure cel_comm_del
  procedure cel_bcast_groups_del
end interface del

interface print
  procedure cel_comm_print
end interface print

type cel_bcast_group_type
  integer(kind=mpik), dimension(:), allocatable :: procs_array
end type cel_bcast_group_type

type cel_comm_type
  integer(kind=mpik)                                    :: comm
  integer(kind=mpik)                                    :: group
  integer(kind=mpik), dimension(:), allocatable         :: send_to
  integer(kind=mpik), dimension(:), allocatable         :: get_from
  integer(kind=ik),   dimension(:), allocatable         :: send_to_first_idx
  integer(kind=ik),   dimension(:), allocatable         :: send_to_last_idx
  integer(kind=ik),   dimension(:), allocatable         :: get_from_first_idx
  integer(kind=ik),   dimension(:), allocatable         :: get_from_last_idx
  integer(kind=ik),   dimension(:,:), allocatable       :: send_com_buffer_perm
  integer(kind=ik),   dimension(:), allocatable         :: send_length_comm_buffer
  integer(kind=ik),   dimension(:,:), allocatable       :: get_com_buffer_perm
  integer(kind=ik),   dimension(:), allocatable         :: get_length_comm_buffer
  integer(kind=mpik), dimension(:), allocatable         :: send_to_transl
  integer(kind=mpik), dimension(:), allocatable         :: get_from_transl
  integer(kind=mpik)                                    :: send_comm
  integer(kind=mpik)                                    :: send_group
  integer(kind=mpik), dimension(:), allocatable         :: get_comms
  integer(kind=mpik), dimension(:), allocatable         :: get_groups
  type(cel_bcast_group_type), allocatable, dimension(:) :: bcast_groups
  integer(kind=mpik), dimension(:), allocatable         :: send_requests
  integer(kind=mpik), dimension(:), allocatable         :: get_requests
  integer(kind=mpik), dimension(:), allocatable         :: vtxdist_displ !to use by scatter
  integer(kind=mpik), dimension(:), allocatable         :: vtxdist_send_counts !to use by scatter and gather
end type cel_comm_type


contains

subroutine cel_comm_print(string,cel_comm,omp)
  character(len=*), intent(in)    :: string
  type(cel_comm_type), intent(in) :: cel_comm
  type(cel_omp_type), intent(in)  :: omp
  integer (kind=ik) :: ii,idx
  character(len=64)  :: string_2

  if(omp%is_master) then
    write(cel_output_unit,'(2a)') trim(string)
    call print("send_to",cel_comm%send_to)
    call print("get_from",cel_comm%get_from)
    call print("send_length_comm_buffer",cel_comm%send_length_comm_buffer)
    call print("get_length_comm_buffer",cel_comm%get_length_comm_buffer)
    do ii=1_ik,size(cel_comm%send_to)
      write(string_2,"(A,I0,A)") "send_com_buffer_perm(:",cel_comm%send_to(ii),")"
      call print(trim(string_2),cel_comm%send_com_buffer_perm(:,ii))
    end do
    do ii=1_ik,size(cel_comm%get_from)
      write(string_2,"(A,I0,A)") "get_com_buffer_perm(:",cel_comm%get_from(ii),")"
      call print(trim(string_2),cel_comm%get_com_buffer_perm(:,ii))
    end do
  end if

end subroutine cel_comm_print

subroutine cel_comm_del(cel_comm)
  type(cel_comm_type), intent(inout) :: cel_comm
  integer(kind=ik) :: ii
  integer(kind=mpik) :: ierr
  if(allocated(cel_comm%send_to)) deallocate(cel_comm%send_to)
  if(allocated(cel_comm%get_from)) deallocate(cel_comm%get_from)
  if(allocated(cel_comm%send_to_transl)) deallocate(cel_comm%send_to_transl)
  if(allocated(cel_comm%get_from_transl)) deallocate(cel_comm%get_from_transl)
  if(allocated(cel_comm%get_comms)) then
    do ii=1,size(cel_comm%get_comms)
      call MPI_Comm_free(cel_comm%get_comms(ii),ierr)
    end do
    deallocate(cel_comm%get_comms)
  end if
  if(allocated(cel_comm%get_groups)) then
    do ii=1,size(cel_comm%get_groups)
      call MPI_Group_free(cel_comm%get_groups(ii),ierr)
    end do
    deallocate(cel_comm%get_groups)
    call MPI_Group_free(cel_comm%group,ierr)
  end if
  if(allocated(cel_comm%send_requests)) deallocate(cel_comm%send_requests)
  if(allocated(cel_comm%get_requests)) deallocate(cel_comm%get_requests)
  if(allocated(cel_comm%vtxdist_displ)) deallocate(cel_comm%vtxdist_displ)
  if(allocated(cel_comm%vtxdist_send_counts)) deallocate(cel_comm%vtxdist_send_counts)
  call del(cel_comm%bcast_groups)

end subroutine cel_comm_del

subroutine cel_bcast_groups_del(cel_bcast_groups)
  type(cel_bcast_group_type), allocatable, dimension(:) :: cel_bcast_groups
  integer(kind=ik) :: ii

  if(allocated(cel_bcast_groups)) then
    do ii=1,size(cel_bcast_groups)
      if(allocated(cel_bcast_groups(ii)%procs_array)) then
        deallocate(cel_bcast_groups(ii)%procs_array)
      end if
    end do
    deallocate(cel_bcast_groups)
  end if

end subroutine cel_bcast_groups_del

subroutine cel_comm_init_vtxdist_counts(cel_comm, vtxdist, cel_omp, err_code)
  type(cel_comm_type), intent(inout)                :: cel_comm
  integer(kind=ik),dimension(:), intent(in)         :: vtxdist
  type(cel_omp_type), intent(in)                    :: cel_omp
  integer(kind=ik), intent(inout)                   :: err_code
  integer(kind=ik):: sys_stat, ii, length
  
  err_code = 0_ik
  if(cel_omp%is_master) then
    if(allocated(cel_comm%vtxdist_send_counts)) then
      if(size(cel_comm%vtxdist_send_counts) .lt. cel_omp%num_masters) then
        deallocate(cel_comm%vtxdist_send_counts)
      end if
    end if
    if(.not. allocated(cel_comm%vtxdist_send_counts)) then
      allocate(cel_comm%vtxdist_send_counts(cel_omp%num_masters))
    end if
    do ii=1_ik,size(vtxdist)-1_ik
      length = vtxdist(ii+1)-vtxdist(ii)
      cel_comm%vtxdist_send_counts(ii)=int(length,kind=mpik)
    end do

    if(allocated(cel_comm%vtxdist_displ)) then
      if(size(cel_comm%vtxdist_displ) .lt. cel_omp%num_masters) then
        deallocate(cel_comm%vtxdist_displ)
      end if
    end if
    if(.not. allocated(cel_comm%vtxdist_displ)) then
      allocate(cel_comm%vtxdist_displ(cel_omp%num_masters))
    end if
    do ii=1_ik,size(vtxdist)-1_ik
      cel_comm%vtxdist_displ(ii)=int(vtxdist(ii)-1_ik,kind=mpik)
    end do
  endif

end subroutine cel_comm_init_vtxdist_counts

!initialize cel_comm%send_to and cel_comm%get_from
!initialize cel_comm%send_to_first_idx and cel_comm%send_to_last_idx
!initialize cel_comm%get_from_first_idx and cel_comm%get_from_last_idx
!initialize cel_comm%get_idx_with_comm_buffer cel_comm%send_idx_with_comm_buffer
!initialize cel_comm%get_com_buffer_perm cel_comm%send_com_buffer_perm
!initialize cel_comm%get_length_comm_buffer cel_comm%gend_length_comm_buffer
!procs_array containes the ranks of the procs having a halo data (inclusive itself at index 1)
!procs_array is initialized in the method cel_sp_mat_distr_from_sp_mat_local
subroutine cel_comm_get_connections(cel_comm, procs_array,&
  procs_array_first_idx, procs_array_last_idx,&
  com_buffer_perm,length_comm_buffer,&
  omp, err_code)
  type(cel_comm_type), intent(inout) :: cel_comm
  integer(kind=ik), dimension(:), allocatable, intent(in)  :: procs_array
  integer(kind=ik), dimension(:), allocatable , intent(in) :: procs_array_first_idx ! array with first idx of data needed from the process in procs_array
  integer(kind=ik), dimension(:), allocatable , intent(in) :: procs_array_last_idx ! array with last idx of data needed from the process in procs_array
  integer(kind=ik), dimension(:,:), allocatable , intent(in) :: com_buffer_perm
  integer(kind=ik), dimension(:), allocatable , intent(in) :: length_comm_buffer
  type(cel_omp_type), intent(in) :: omp
  integer (kind=ik), intent(out) :: err_code
  integer (kind=ik) :: ii, idx, num, max_size, num_procs_to_get, size_comm_buffer
  integer(kind=mpik), dimension(:), allocatable :: send_to_proc_array
  integer(kind=ik), dimension(:), allocatable :: send_to_proc_array_first_idx
  integer(kind=ik), dimension(:), allocatable :: send_to_proc_array_last_idx
  integer(kind=ik), dimension(:), allocatable :: send_length_comm_buffer
  integer(kind=mpik)  :: ierr, mpi_win_disp_unit, mpi_win_disp_unit_idx
  integer(kind=mpik)  :: mpi_win_handle, mpi_win_handle_first_idx,mpi_win_handle_last_idx
  integer(kind=mpik)  :: mpi_win_handle_comm_buffer
  integer(MPI_ADDRESS_KIND) :: target_displ, mpi_mem_size, mpi_mem_size_idx
  integer(kind=mpik) :: one
  err_code = 0_ik
  
  if(omp%is_master) then
    cel_comm%comm = omp%master_comm
    num_procs_to_get = size(procs_array)
    size_comm_buffer = size(com_buffer_perm,1_ik)
    mpi_mem_size = INT(omp%num_masters,kind=MPI_ADDRESS_KIND)*mpik
    mpi_mem_size_idx = INT(omp%num_masters,kind=MPI_ADDRESS_KIND)*ik
    mpi_win_disp_unit = INT(mpik,kind=MPI_ADDRESS_KIND)
    mpi_win_disp_unit_idx = INT(ik,kind=MPI_ADDRESS_KIND)
    allocate(send_to_proc_array(omp%num_masters))
    allocate(send_to_proc_array_first_idx(omp%num_masters))
    allocate(send_to_proc_array_last_idx(omp%num_masters))
    allocate(send_length_comm_buffer(omp%num_masters))
    allocate(cel_comm%get_from(num_procs_to_get))
    allocate(cel_comm%get_from_first_idx(num_procs_to_get))
    allocate(cel_comm%get_from_last_idx(num_procs_to_get))
    allocate(cel_comm%get_com_buffer_perm(size(com_buffer_perm,1),num_procs_to_get))
    allocate(cel_comm%get_length_comm_buffer(size(length_comm_buffer)))
    send_to_proc_array = 0_mpik
    send_length_comm_buffer = 0_mpik
    cel_comm%get_from = INT(procs_array,kind=mpik)
    cel_comm%get_from_first_idx = INT(procs_array_first_idx,kind=ik)
    cel_comm%get_from_last_idx = INT(procs_array_last_idx,kind=ik)
    cel_comm%get_length_comm_buffer(1:num_procs_to_get)=length_comm_buffer(1:num_procs_to_get)
    do ii=2_ik,num_procs_to_get
      num = length_comm_buffer(ii)
      cel_comm%get_com_buffer_perm(1:num,ii)=&
        com_buffer_perm(1:num,ii)
    end do

    one = 1_mpik
    call MPI_Win_create(send_to_proc_array, mpi_mem_size, mpi_win_disp_unit, MPI_INFO_NULL,&
      cel_comm%comm, mpi_win_handle,ierr)
    call MPI_Win_create(send_to_proc_array_first_idx, mpi_mem_size_idx, mpi_win_disp_unit_idx,&
      MPI_INFO_NULL,&
      cel_comm%comm, mpi_win_handle_first_idx,ierr)
    call MPI_Win_create(send_to_proc_array_last_idx, mpi_mem_size_idx, mpi_win_disp_unit_idx,&
      MPI_INFO_NULL,&
      cel_comm%comm, mpi_win_handle_last_idx,ierr)
    call MPI_Win_create(send_length_comm_buffer, mpi_mem_size_idx, mpi_win_disp_unit_idx,&
      MPI_INFO_NULL,&
      cel_comm%comm, mpi_win_handle_comm_buffer,ierr)
    
    call MPI_Win_Fence(0_mpik,mpi_win_handle,ierr)
    call MPI_Win_Fence(0_mpik,mpi_win_handle_first_idx,ierr)
    call MPI_Win_Fence(0_mpik,mpi_win_handle_last_idx,ierr)
    call MPI_Win_Fence(0_mpik,mpi_win_handle_comm_buffer,ierr)
    do ii=2, size(procs_array)
      target_displ = int(omp%master_num,kind=MPI_ADDRESS_KIND)
      CALL MPI_Put( one, 1_mpik, MPI_INT, &
        cel_comm%get_from(ii),target_displ,1_mpik,MPI_INT,mpi_win_handle,ierr )
      CALL MPI_Put( cel_comm%get_from_first_idx(ii), 1_mpik, omp%mpi_type_ik, &
        cel_comm%get_from(ii),target_displ,1_mpik,omp%mpi_type_ik,mpi_win_handle_first_idx,ierr)
      CALL MPI_Put( cel_comm%get_from_last_idx(ii), 1_mpik, omp%mpi_type_ik, &
        cel_comm%get_from(ii),target_displ,1_mpik,omp%mpi_type_ik,mpi_win_handle_last_idx,ierr)
      if( cel_comm%get_length_comm_buffer(ii) .gt. 0) then
        CALL MPI_Put( cel_comm%get_length_comm_buffer(ii), 1_mpik, omp%mpi_type_ik, &
          cel_comm%get_from(ii),target_displ,1_mpik,omp%mpi_type_ik,mpi_win_handle_comm_buffer,ierr )
      end if
    end do
    call MPI_Win_Fence(0_mpik,mpi_win_handle,ierr)
    call MPI_Win_Fence(0_mpik,mpi_win_handle_first_idx,ierr)
    call MPI_Win_Fence(0_mpik,mpi_win_handle_last_idx,ierr)
    call MPI_Win_Fence(0_mpik,mpi_win_handle_comm_buffer,ierr)
    call cel_error_proc("Error in cel_comm_get_connections: MPI_Win_Fence", &
      int(omp%master_num,kind=mpik), ierr, cel_is_in_debug, .FALSE._lk)
    call MPI_Win_Free(mpi_win_handle,ierr)
    call MPI_Win_Free(mpi_win_handle_first_idx,ierr)
    call MPI_Win_Free(mpi_win_handle_last_idx,ierr)
    call MPI_Win_Free(mpi_win_handle_comm_buffer,ierr)
    call cel_error_proc("Error incel_comm_get_connections: MPI_Win_Free", &
      int(omp%master_num,kind=mpik), ierr, cel_is_in_debug, .FALSE._lk)
    num = SUM(send_to_proc_array)
    allocate(cel_comm%send_to(num+1_ik))
    allocate(cel_comm%send_to_first_idx(num+1_ik))
    allocate(cel_comm%send_to_last_idx(num+1_ik))
    
    max_size = MAXVAL(send_length_comm_buffer)
    if(num .gt. 0) then
      allocate(cel_comm%send_com_buffer_perm(max_size,num+1_ik))
      allocate(cel_comm%send_length_comm_buffer(num+1_ik))
     ! allocate(cel_comm%send_idx_with_comm_buffer(num+1_ik))
    else
      allocate(cel_comm%send_com_buffer_perm(0_ik,num+1_ik))
      allocate(cel_comm%send_length_comm_buffer(0_ik))
     ! allocate(cel_comm%send_idx_with_comm_buffer(0_ik))
    end if
    cel_comm%send_to = -1_ik
    cel_comm%send_to_first_idx = -1_ik
    cel_comm%send_to_last_idx = -1_ik
    cel_comm%send_com_buffer_perm = -1_ik
   ! if( omp%master_num .eq. 0) then
     ! call print("send_length_comm_buffer:",send_length_comm_buffer)
   ! end if
    idx=1_ik
    cel_comm%send_to(idx) = int(omp%master_num,kind=mpik)
    cel_comm%send_length_comm_buffer=0_ik
    do ii=1, size(send_to_proc_array)
      if(send_to_proc_array(ii) .gt. 0 .and. ii .ne. omp%master_num+1_ik) then
        idx = idx + 1_ik
        cel_comm%send_to(idx)=int(ii-1_mpik,kind=mpik)
        cel_comm%send_to_first_idx(idx)=send_to_proc_array_first_idx(ii)
        cel_comm%send_to_last_idx(idx)=send_to_proc_array_last_idx(ii)
        cel_comm%send_length_comm_buffer(idx)=send_length_comm_buffer(ii)
      end if
    end do
    call MPI_BARRIER(cel_comm%comm, ierr)
    call cel_comm_get_comm_buffer_connections(cel_comm, omp, err_code)
  end if
  !$omp barrier
end subroutine cel_comm_get_connections

subroutine cel_comm_get_comm_buffer_connections(cel_comm, omp, err_code)
  type(cel_comm_type), intent(inout) :: cel_comm
  type(cel_omp_type), intent(in) :: omp
  integer (kind=ik), intent(out) :: err_code
  integer (kind=ik) :: ii,num_get,num_send
  integer (kind=mpik) :: dest_mpik, root,length, ierr, tag_mpik
  err_code = 0_ik
  if(omp%is_master) then
    !Send comm_buffer
    num_get = size(cel_comm%get_from)
    num_send = size(cel_comm%send_to)
    root = int(omp%master_num,kind=mpik)
    if(.not. allocated(cel_comm%send_requests)) allocate(cel_comm%send_requests(num_get-1_ik))
    cel_comm%send_requests = MPI_REQUEST_NULL
    do ii=2, num_get
      dest_mpik = int(cel_comm%get_from(ii),kind=mpik)
      length=int(cel_comm%get_length_comm_buffer(ii),kind=mpik)
      if(length .gt. 0_ik) then
        tag_mpik = root + 3_mpik*dest_mpik
        call MPI_Isend(cel_comm%get_com_buffer_perm(1:length,ii),&
               length,omp%mpi_type_ik,&
               dest_mpik,tag_mpik,cel_comm%comm,cel_comm%send_requests(ii-1_mpik),ierr)
      end if
    end do
    !Get comm_buffer
    if(.not. allocated(cel_comm%get_requests)) allocate(cel_comm%get_requests(num_send-1_ik))
    cel_comm%get_requests = MPI_REQUEST_NULL
    do ii=2, num_send
      length=int(cel_comm%send_length_comm_buffer(ii),kind=mpik)
      !cel_comm%send_com_buffer_perm(:,:) = 0_ik
      if(length .gt. 0_ik) then
        dest_mpik = int(cel_comm%send_to(ii),kind=mpik)
        tag_mpik = dest_mpik + 3_mpik*root
        call MPI_Irecv(cel_comm%send_com_buffer_perm(1:length,ii),&
              length,omp%mpi_type_ik,&
              dest_mpik,tag_mpik,cel_comm%comm,cel_comm%get_requests(ii-1_mpik),ierr)
     end if
    end do
    if (num_send .gt. 1) then
      call MPI_Waitall(num_send-1_ik, &
         cel_comm%get_requests(1_ik:num_send-1_ik), MPI_STATUSES_IGNORE, ierr)
    end if
  end if
  deallocate(cel_comm%get_requests)
  deallocate(cel_comm%send_requests)
end subroutine

!initialize bcast communicators/groups in the cel_comm for vector distribution
subroutine cel_comm_get_communicators(cel_comm, omp, err_code)
  type(cel_comm_type), intent(inout) :: cel_comm
  type(cel_omp_type), intent(in) :: omp
  integer (kind=ik), intent(out) :: err_code
  integer(kind=mpik)  :: ierr, num
  integer(kind=ik) ::ii, jj,  idx,idx_get_from
  logical(kind=lk) :: flag, flag_two
  integer(kind=mpik), dimension(:), allocatable :: group_proc_array_size
  integer(kind=mpik)  :: mpi_win_handl,  mpi_win_disp_unit, tmp_comm, group_size
  integer(MPI_ADDRESS_KIND) :: target_displ, mpi_mem_size
  integer(kind=mpik),allocatable, dimension(:)  ::send_request
  integer(kind=mpik),allocatable, dimension(:)  ::get_requests
  integer(kind=ik), dimension(:), allocatable :: get_from_sort_order
  err_code = 0_ik
  if(omp%is_master) then
    !send the size of the bcast group for x-vector distribution to all group members
    !
    !group_proc_array_size - array with the lengths of the bcast groups, 
    !                        in which the origin process is involved (if > 0) or not (if == 0)
    !                        All processes are the receiver in some bcast groups,
    !                        whose roots are the processes with the rank cel_comm%get_from(i)
    !                        In the cel_comm%get_from(1) is always saved the rank of the origin process self.
    allocate(group_proc_array_size(omp%num_masters))
    group_proc_array_size = 0_mpik
    mpi_mem_size = INT(size(group_proc_array_size)*mpik,kind=MPI_ADDRESS_KIND)
    mpi_win_disp_unit = INT(mpik,kind=MPI_ADDRESS_KIND)
    call MPI_Win_create(group_proc_array_size, mpi_mem_size, mpi_win_disp_unit, &
      MPI_INFO_NULL,&
      cel_comm%comm, mpi_win_handl,ierr)
    target_displ = int(omp%master_num,kind=MPI_ADDRESS_KIND)
    num = int(size(cel_comm%send_to),kind=mpik)
    call MPI_Win_Fence(0_mpik,mpi_win_handl,ierr)
    do ii = 1, size(cel_comm%send_to)
      !call MPI_Win_lock(MPI_LOCK_SHARED,INT(cel_comm%send_to(ii),kind=mpik),0_mpik,mpi_win,ierr)
      CALL MPI_Put( num, 1_mpik, MPI_INT, &
        cel_comm%send_to(ii),target_displ,1_mpik,MPI_INT,mpi_win_handl,ierr)
      !call MPI_Win_unlock(INT(cel_comm%send_to(ii),kind=mpik),mpi_win,ierr)
    end do
    call MPI_Win_Fence(0_mpik,mpi_win_handl,ierr)
    call MPI_Win_Free(mpi_win_handl,ierr)
    !bcast_groups - array wtih the bcast groups in which the current process is involved
    allocate(cel_comm%get_comms(size(cel_comm%get_from)))
    allocate(cel_comm%get_groups(size(cel_comm%get_from)))
    allocate(cel_comm%bcast_groups(size(cel_comm%get_from)))
    do ii=1, size(cel_comm%get_from)
      allocate(cel_comm%bcast_groups(ii)%&
        procs_array(group_proc_array_size(cel_comm%get_from(ii)+1_mpik)))
    end do
    allocate(send_request(size(cel_comm%send_to)))
    do ii=1, size(cel_comm%send_to)
      call MPI_Isend(cel_comm%send_to, size(cel_comm%send_to),&
        MPI_INT,cel_comm%send_to(ii),&
        int(cel_comm%send_to(ii)*1+omp%master_num*3,kind=mpik),&
        cel_comm%comm,send_request(ii),ierr)
    end do
    allocate(get_requests(size(cel_comm%get_from)))
    do ii=1, size(cel_comm%get_from)
      call MPI_Irecv(cel_comm%bcast_groups(ii)%procs_array, &
        size(cel_comm%bcast_groups(ii)%procs_array),MPI_INT,&
        cel_comm%get_from(ii),&
        int(omp%master_num*1+cel_comm%get_from(ii)*3,kind=mpik),&
        cel_comm%comm,get_requests(ii),ierr)
    end do
    call MPI_Waitall(size(cel_comm%get_from), &
      get_requests,MPI_STATUSES_IGNORE,ierr)
    !inizialize the bcast groups and communicators in mpi
    call cel_comm_bubble_sort_seq(cel_comm%get_from,get_from_sort_order)
    call MPI_Comm_group(cel_comm%comm,cel_comm%group,ierr)
    idx = 1_ik
    flag_two = .false._lk
    do ii=1_ik, omp%num_masters
      flag = .false._lk
      do jj=idx, size(get_from_sort_order)
        idx_get_from = get_from_sort_order(jj)
        if(cel_comm%get_from(idx_get_from) + 1_ik .eq. ii) then
          flag = .true._lk
          idx = jj
          if(flag_two) then
           call cel_error_proc("Error cel_comm_get_communicators: two different bcast_send_groups!", &
            int(omp%master_num,kind=mpik), ierr, cel_is_in_debug, .FALSE._lk)
            err_code = 1_ik
          end if
          flag_two = .true._lk
        end if
        if(flag) then
          if(cel_comm%bcast_groups(idx_get_from)%procs_array(1_ik) .eq. omp%master_num) then
            group_size = size(cel_comm%bcast_groups(idx_get_from)%procs_array)
            !call print("ranks",cel_comm%bcast_groups(idx_get_from)%procs_array(1:group_size))
            call MPI_Group_incl(cel_comm%group,&
              group_size,&
              cel_comm%bcast_groups(idx_get_from)%procs_array(1:group_size),&
              cel_comm%send_group,ierr)
              call MPI_Comm_create(cel_comm%comm, cel_comm%send_group,&
                cel_comm%send_comm,ierr)
              call cel_error_proc(&
               "Error cel_comm_get_communicators: could not create send_comm!", &
                int(omp%master_num,kind=mpik), int(ierr,kind=mpik), cel_is_in_debug, .FALSE._lk)
              err_code = 2_ik
              cel_comm%get_comms(1) = cel_comm%send_comm
              cel_comm%get_groups(1)  = cel_comm%send_group
          else
            group_size = size(cel_comm%bcast_groups(idx_get_from)%procs_array)
            !call print("ranks",cel_comm%bcast_groups(idx_get_from)%procs_array(1:group_size))
            call MPI_Group_incl(cel_comm%group,&
              group_size,&
              cel_comm%bcast_groups(idx_get_from)%procs_array(1:group_size),&
              cel_comm%get_groups(idx_get_from),ierr)
            call MPI_Comm_create(cel_comm%comm,cel_comm%get_groups(idx_get_from),&
              cel_comm%get_comms(idx_get_from),ierr)
            if(cel_comm%get_comms(idx_get_from) .eq. MPI_COMM_NULL) then
              call cel_error_proc(&
               "Error cel_comm_get_communicators: could not create bcast_get_comm!", &
                int(omp%master_num,kind=mpik), int(ierr,kind=mpik), cel_is_in_debug, .FALSE._lk)
              err_code = 3_ik
            end if
          end if
          exit
        end if
      end do
      if(.not. flag) then
        call MPI_Comm_create(cel_comm%comm,MPI_GROUP_EMPTY,&
              tmp_comm,ierr)
      end if
    end do
    call MPI_Barrier(cel_comm%comm,ierr)
    allocate(cel_comm%send_to_transl(size(cel_comm%send_to)))
    do ii=1_ik,size(cel_comm%send_to)
      call MPI_Group_Translate_ranks(cel_comm%group,1_mpik,cel_comm%send_to(ii:ii),&
        cel_comm%send_group,cel_comm%send_to_transl(ii:ii),ierr)
      call cel_error_proc(&
               "Error cel_comm_get_communicators: MPI_Group_Translate_ranks send!", &
                int(omp%master_num,kind=mpik), int(ierr,kind=mpik), cel_is_in_debug, .FALSE._lk)
    end do
    allocate(cel_comm%get_from_transl(size(cel_comm%get_from)))
    do ii=1_ik,size(cel_comm%get_from)
      call MPI_Group_Translate_ranks(cel_comm%group,1_mpik,cel_comm%get_from(ii:ii),&
      cel_comm%get_groups(ii), cel_comm%get_from_transl(ii:ii),ierr)
      call cel_error_proc(&
               "Error cel_comm_get_communicators: MPI_Group_Translate_ranks get!", &
                int(omp%master_num,kind=mpik), int(ierr,kind=mpik), cel_is_in_debug, .FALSE._lk)
    end do
    call cel_comm_check_communicators_master(cel_comm, omp, err_code)
     call cel_error_proc(&
               "Error cel_comm_get_communicators: test not passed!", &
                int(omp%master_num,kind=mpik), int(ierr,kind=mpik), cel_is_in_debug, .TRUE._lk)
  end if


end subroutine cel_comm_get_communicators
!verify the bcast communicators
subroutine cel_comm_check_communicators_master(cel_comm, omp, err_code)
  type(cel_comm_type), intent(in) :: cel_comm
  type(cel_omp_type), intent(in) :: omp
  integer (kind=ik), intent(out) :: err_code
  integer(kind=mpik), dimension(:), allocatable         :: bcast_comm_ranks
  integer(kind=mpik) :: master_comm_group, ierr, group, group_size
  integer (kind=ik) :: ii, jj
  err_code = 0_ik
  call MPI_Comm_group(cel_comm%comm,master_comm_group,ierr)
  !check bcast send communicator
  call MPI_Comm_group(cel_comm%send_comm,group,ierr)
  call MPI_Group_size(group,group_size,ierr)
  if(group_size .ne. size(cel_comm%send_to)) then
    call cel_error_proc(&
      "Error cel_comm_check_communicators_master: length of send_comm. is wrong!", &
            int(omp%master_num,kind=mpik), 1_mpik, cel_is_in_debug, .FALSE._lk)
    err_code = 1_ik
  else
    allocate(bcast_comm_ranks(size(cel_comm%send_to)))
    call MPI_Group_translate_ranks(master_comm_group,group_size,cel_comm%send_to,&
      group,bcast_comm_ranks,ierr)
    do ii = 1_ik, group_size
      if(bcast_comm_ranks(ii) .ne. ii - 1_mpik) then
        call cel_error_proc(&
          "Error cel_comm_check_communicators_master: ranks of send_comm. are wrong!", &
          int(omp%master_num,kind=mpik), 1_mpik, cel_is_in_debug, .FALSE._lk)
        err_code = 2_ik
      end if
    end do
  end if
  call MPI_Group_free(group,ierr)
  !check get communicators
  do jj = 1, size(cel_comm%get_comms)
    call MPI_Comm_group(cel_comm%get_comms(jj),group,ierr)
    call MPI_Group_size(group,group_size,ierr)
    if(group_size .ne. size(cel_comm%bcast_groups(jj)%procs_array)) then
      call cel_error_proc(&
        "Error cel_comm_check_communicators_master: length of bcast_get_comm. is wrong!", &
              int(omp%master_num,kind=mpik), 1_mpik, cel_is_in_debug, .FALSE._lk)
      err_code = 3_ik
    else
      if(allocated(bcast_comm_ranks)) then
        if(size(bcast_comm_ranks) .ne. group_size) then
          deallocate(bcast_comm_ranks)
          allocate(bcast_comm_ranks(group_size))
        end if
      else
        allocate(bcast_comm_ranks(group_size))
      end if
      call MPI_Group_translate_ranks(master_comm_group,&
        group_size,&
        cel_comm%bcast_groups(jj)%procs_array,&
        group,&
        bcast_comm_ranks,ierr)
      do ii = 1_ik, group_size
        if(bcast_comm_ranks(ii) .ne. ii - 1_mpik) then
          call cel_error_proc(&
            "Error cel_comm_check_communicators_master: ranks of bcast_get_comm. are wrong!", &
            int(omp%master_num,kind=mpik), 1_mpik, cel_is_in_debug, .FALSE._lk)
          err_code = 4_ik
        end if
      end do
    end if
    call MPI_Group_free(group,ierr)
  end do
  call MPI_Group_free(master_comm_group,ierr)
  
end subroutine cel_comm_check_communicators_master

subroutine cel_comm_bubble_sort_seq(array, order)
 
  integer(kind=mpik), dimension(:),allocatable,  intent(in) :: array
  integer(kind=ik), dimension(:),allocatable,  intent(out) :: order
  integer(kind=ik) :: local_length
  integer(kind=ik) :: temp
  logical(kind=lk) :: swapped = .TRUE._lk
  integer(kind=ik) :: ii, jj

  local_length = size(array)
  allocate(order(local_length))
  do ii=1,local_length
    order(ii)=ii
  end do
  do jj = local_length-1_ik, 1_ik, -1_ik
    swapped = .FALSE._lk
    do ii = 1, jj
      if(array(order(ii)) > array(order(ii+1_ik))) then
        temp=order(ii)
        order(ii) = order(ii+1_ik)
        order(ii+1_ik) = temp
        swapped = .TRUE._lk
      end if
    end do
    if(.NOT. swapped) exit
  end do
end subroutine cel_comm_bubble_sort_seq


subroutine cel_comm_show_send_group_graphvizio(dir, filename_prefix, filename_ext,&
  cel_comm, omp, err_code)
  character(len=*), intent(in) :: dir
  character(len=*), intent(in) :: filename_prefix
  character(len=*), intent(in) :: filename_ext
  type(cel_comm_type), intent(inout) :: cel_comm
  type(cel_omp_type), intent(in) :: omp
  integer (kind=ik), intent(out) :: err_code
  integer(kind=ik) :: ii
  integer(kind=mpik) :: ierr
  character(len=1024) :: path
  integer :: unit,iostat

  err_code = 0_ik
  if(omp%is_master) then
    unit=3021
    write(path,'(A,A,A,I0,A,A)') trim(dir),trim("/"), trim(filename_prefix), &
      omp%master_num, trim("."), trim(filename_ext)
    if(omp%master_num .eq. 0 ) then
      write(cel_error_unit,'(A,A)') "graph output directory:",trim(path)
      call system('mkdir -p ' // trim(dir) )
      call system('sleep 0.1 ')
    end if
    call MPI_BARRIER(omp%master_comm, ierr)
    open(unit=unit,iostat=iostat,file=path,status='REPLACE',action='READWRITE')
    if(iostat .gt. 0) then
      call cel_error_proc(&
       "Error cel_comm_show_send_group_graphvizio: could not open a file", &
       int(omp%master_num,kind=mpik), int(iostat,kind=mpik), cel_is_in_debug, .FALSE._lk)
    endif
    write(unit,'(A,I0)',advance='yes')"subgraph cluster",omp%master_num
    write(unit,'(A)',advance='yes')"{"
    if(size(cel_comm%send_to) .gt. 1) then
        write(unit,'(A,I0,A,A,A,I0,A,A)',advance='yes') "root",omp%master_num,&
            "[label=",CHAR(34),"proc ",omp%master_num,CHAR(34),"];"
        do ii=2,size(cel_comm%send_to)
          write(unit,'(A,I0,A,I0,A,A,A,I0,A,A)',advance='yes') "node",omp%master_num,"_",cel_comm%send_to(ii),&
            "[label=",CHAR(34),"proc ",cel_comm%send_to(ii),CHAR(34),"];"
        end do
        write(unit,'(A,I0,A,I0,A)',advance='no') "root",omp%master_num ," -> {"
        do ii=2,size(cel_comm%send_to)
          write(unit,'(A,I0,A,I0)',advance='no') "node",omp%master_num,"_",cel_comm%send_to(ii)
          if(ii .lt. size(cel_comm%send_to)) then
            write(unit,'(A)',advance='no') "; "
          else
            write(unit,'(A)',advance='yes') "} [arrowsize=1, penwidth=2,len=1.4];"
          end if
        end do
    end if
    write(unit,'(A)',advance='yes') "}"
    close(unit=unit)
  end if
end subroutine cel_comm_show_send_group_graphvizio

end module cel_comm_module
