!     
! File:   cel_sp_mat_distr_vec_module.f90
! Project CRESTA (see details on https://cresta-project.eu) Exascale library
! Send you email about the bugs to
! The High Performance Computing Center Stuttgart (HLRS) of the University of Stuttgart
! Dmitry Khabi,  (email: khabi@hlrs.de)
! Created on NOVEMBER 24, 2013
!

module cel_sp_mat_distr_vec_module
use cel_types_module
use cel_timer_interface_module
use cel_error_module
use cel_base_print_module
use cel_vec_module
use cel_omp_module
use cel_perf_module
use cel_comm_module
use cel_omp_shared_work_module
use OMP_LIB
use MPI
use, intrinsic :: ISO_C_BINDING
implicit none

interface new
  module procedure cel_sp_mat_distr_vec_new_master
  module procedure cel_sp_mat_distr_vec_new_from_procs
end interface new

interface del
  module procedure cel_sp_mat_distr_vec_del
  module procedure cel_sp_mat_distr_vec_array_del
end interface del

interface print
  module procedure cel_sp_mat_distr_vec_print
end interface print

type cel_sp_mat_distr_vec_type
  real(kind=rk), allocatable, dimension(:,:) ::  values! values array of the vec
  real(kind=rk), allocatable, dimension(:,:) ::  send_comm_buffer! values array for the communication
  real(kind=rk), allocatable, dimension(:,:) ::  get_comm_buffer! values array for the communication
  integer(kind=ik) :: num_elements
  integer(kind=ik) :: num_distr_vectors
  logical(kind=lk) :: is_comm_buffer_sync
end type cel_sp_mat_distr_vec_type

contains

subroutine cel_sp_mat_distr_vec_new_master(distr_vectors,num_vectors,&
  num_elements,length_vectors)
  type(cel_sp_mat_distr_vec_type), intent(inout)  :: distr_vectors
  integer(kind=ik), intent(in)    :: num_vectors
  integer(kind=ik), intent(in)    :: num_elements
  integer(kind=ik), intent(in)    :: length_vectors
  integer(kind=ik) ::sys_stat
  
  call del(distr_vectors)
  if(num_vectors .gt. 0_ik .and. length_vectors .gt. 0_ik) then
    allocate(distr_vectors%values(length_vectors,num_vectors), stat=sys_stat)
    distr_vectors%is_comm_buffer_sync = .true._lk
    call cel_error("cel_sp_mat_distr_vec_new distr_vectors%values ", sys_stat,&
     cel_is_in_debug,.TRUE._lk)
    distr_vectors%num_distr_vectors=num_vectors
    distr_vectors%num_elements = num_elements
  endif

end subroutine cel_sp_mat_distr_vec_new_master

subroutine cel_sp_mat_distr_vec_del(distr_vectors)
  type(cel_sp_mat_distr_vec_type), intent(inout)  :: distr_vectors
  integer(kind=ik) :: ii
  
  if(allocated(distr_vectors%values)) deallocate(distr_vectors%values)
  if(allocated(distr_vectors%send_comm_buffer)) deallocate(distr_vectors%send_comm_buffer)
  if(allocated(distr_vectors%get_comm_buffer)) deallocate(distr_vectors%get_comm_buffer)
  distr_vectors%num_distr_vectors=0_ik
  
end subroutine cel_sp_mat_distr_vec_del


subroutine cel_sp_mat_distr_vec_array_del(distr_vectors)
  type(cel_sp_mat_distr_vec_type),allocatable,dimension(:), intent(inout)  :: distr_vectors
  integer(kind=ik) :: ii, length
  
  if(allocated(distr_vectors)) then
    length = size(distr_vectors)
    do ii=1_ik, length
      if(allocated(distr_vectors(ii)%values)) deallocate(distr_vectors(ii)%values)
      if(allocated(distr_vectors(ii)%send_comm_buffer)) deallocate(distr_vectors(ii)%send_comm_buffer)
      if(allocated(distr_vectors(ii)%get_comm_buffer)) deallocate(distr_vectors(ii)%get_comm_buffer)
    end do
    deallocate(distr_vectors)
  end if
  
end subroutine cel_sp_mat_distr_vec_array_del

subroutine cel_sp_mat_distr_vec_print_sorted(string,vector,&
                                cel_omp,cel_comm,&
                                output_on,err_code)
  character(len=*), intent(in)                                   :: string
  type(cel_sp_mat_distr_vec_type), intent(inout)                    :: vector
  type(cel_omp_type), intent(in)                                    :: cel_omp
  type(cel_comm_type), intent(in)                                   :: cel_comm
  logical(kind=lk)  , intent(in)                                    :: output_on
  integer(kind=ik), intent(out)                                     :: err_code
  integer(kind=ik) :: ii, jj, kk, vec_length, num_vecs
  real(kind=trk) :: time
  integer(kind=mpik) :: ierr
  character(len=STRING_MAX_LENGTH) :: proc_string
  err_code = 0_ik
  vec_length = vector%num_elements
  num_vecs = size(vector%values,dim=2)
  if(cel_omp%is_master) then
    do ii=1,cel_omp%num_masters
      do jj=1,cel_omp%num_masters
        if(ii .eq. jj) then
          do kk=1_ik,num_vecs
            write(proc_string,"(I3,A,A,A,I3)") &
              cel_omp%master_num,trim(": "),trim(string),"; related proc:",&
              cel_comm%get_from(kk)
            call print(proc_string,vector%values(:,kk))
            call sleep(1)
          end do 
        else
          call MPI_Barrier(cel_omp%master_comm, ierr)
        end if
      end do
     end do
  end if
  !$omp barrier

end subroutine cel_sp_mat_distr_vec_print_sorted




subroutine cel_sp_mat_distr_vec_print(string,distr_vectors,omp,output_on)
  character(len=*) :: string
  type(cel_sp_mat_distr_vec_type), intent(inout) :: distr_vectors
  type(cel_omp_type), intent(in)                 :: omp
  logical(kind=lk)  ::output_on
  character(len=512) :: label
  integer(kind=ik) :: ii,num_buf
  
  if(omp%is_master .and. output_on) then
    write(label,'(A,I0,A,I0,A)') trim(string),omp%master_num,&
      ": distr_vectors [",distr_vectors%num_distr_vectors,"]"
    call print(label,distr_vectors%values)
    if(allocated (distr_vectors%send_comm_buffer)) then
      call print("send_comm_buffer",distr_vectors%send_comm_buffer)
    end if
    if(allocated (distr_vectors%send_comm_buffer)) then
      call print("get_comm_buffer",distr_vectors%get_comm_buffer)
    end if
  end if
    
end subroutine cel_sp_mat_distr_vec_print

subroutine cel_sp_mat_distr_local_vec(distr_vectors,local_vector,vtxdist,cel_omp)
  type(cel_sp_mat_distr_vec_type), intent(inout) :: distr_vectors
  real(kind=rk), dimension(:), intent(in)        :: local_vector
  integer(kind=ik),dimension(:), intent(in)      :: vtxdist
  type(cel_omp_type), intent(in)                 :: cel_omp
  integer(kind=ik)      :: length
  
  if(cel_omp%is_master) then
    length = vtxdist(cel_omp%master_num+2_ik)-vtxdist(cel_omp%master_num+1_ik)
    distr_vectors%values(1_ik:length,1_ik) = local_vector(1_ik:length)
  end if
end subroutine cel_sp_mat_distr_local_vec

subroutine cel_sp_mat_distr_global_vec(distr_vectors,global_vector,cel_comm,vtxdist,cel_omp)
  type(cel_sp_mat_distr_vec_type), intent(inout) :: distr_vectors
  real(kind=rk), dimension(:), intent(in)        :: global_vector
  type(cel_comm_type), intent(inout)             :: cel_comm
  integer(kind=ik),dimension(:), intent(in)      :: vtxdist
  type(cel_omp_type), intent(in)                 :: cel_omp
  integer(kind=ik)      :: length, gl_vector_length
  integer(kind=ik) :: ii
  integer(kind=mpik) :: ierr,length_mpik
  real(kind=rk) :: tmp_buff

  if(cel_omp%is_master) then
    gl_vector_length = size(global_vector)
    length = vtxdist(cel_omp%master_num+2_ik)-vtxdist(cel_omp%master_num+1_ik)
    length_mpik = int(length,kind=mpik)
    if(cel_omp%master_num .eq. 0_mpik) then
        call MPI_Scatterv(global_vector(1_ik:gl_vector_length),&
                          cel_comm%vtxdist_send_counts(1_ik:cel_omp%num_masters),&
                          cel_comm%vtxdist_displ(1_ik:cel_omp%num_masters),&
                          cel_omp%mpi_type_rk,&
                          distr_vectors%values(1_ik:length,1_ik),&
                          length_mpik,&
                          cel_omp%mpi_type_rk,&
                          0_mpik,&
                          cel_omp%master_comm,ierr)
    else
        call MPI_Scatterv(tmp_buff,&
                          cel_comm%vtxdist_send_counts(1_ik:cel_omp%num_masters),&
                          cel_comm%vtxdist_displ(1_ik:cel_omp%num_masters),&
                          cel_omp%mpi_type_rk,&
                          distr_vectors%values(1_ik:length,1_ik),&
                          length_mpik,&
                          cel_omp%mpi_type_rk,&
                          0_mpik,&
                          cel_omp%master_comm,ierr)
    end if
 end if

end subroutine cel_sp_mat_distr_global_vec

subroutine cel_sp_mat_distr_collect_global_vec(global_vector,distr_vectors,cel_comm,vtxdist,cel_omp,err_code)
  real(kind=rk), allocatable, dimension(:), intent(inout)  :: global_vector
  type(cel_sp_mat_distr_vec_type), intent(in)              :: distr_vectors
  type(cel_comm_type), intent(inout)                       :: cel_comm
  integer(kind=ik),dimension(:), intent(in)                :: vtxdist
  type(cel_omp_type), intent(in)                           :: cel_omp
  integer(kind=ik), intent(inout)                          :: err_code
  integer(kind=ik)      :: length, start_indx, end_indx, gl_vector_length
  integer(kind=ik) :: ii, sys_stat
  integer(kind=mpik) :: ierr,length_mpik
  real(kind=rk) :: tmp_buff
  
  err_code = 0_ik
  if(cel_omp%is_master) then
    gl_vector_length = vtxdist(cel_omp%num_masters+1_ik) - vtxdist(1_ik)
    length = vtxdist(cel_omp%master_num+2_ik)-vtxdist(cel_omp%master_num+1_ik)
    length_mpik = int(length,kind=mpik)
    if(cel_omp%master_num .eq. 0_mpik) then
      if(allocated(global_vector)) then
        if(size(global_vector) .lt. gl_vector_length) then
          deallocate(global_vector)
        end if
      end if
      if(.not. allocated(global_vector)) then
        allocate(global_vector(gl_vector_length), stat=sys_stat)
        call cel_error("Error cel_sp_mat_distr_collect_global_vec alloc ", sys_stat,&
                        .TRUE._lk,.TRUE._lk)
      end if
      call MPI_Gatherv(distr_vectors%values(1_ik:length,1_ik),&
                       length_mpik,&
                       cel_omp%mpi_type_rk,&
                       global_vector(1_ik:gl_vector_length),&
                       cel_comm%vtxdist_send_counts(1_ik:cel_omp%num_masters),&
                       cel_comm%vtxdist_displ(1_ik:cel_omp%num_masters),&
                       cel_omp%mpi_type_rk,&
                       0_mpik,&
                       cel_omp%master_comm,ierr)
    else
      call MPI_Gatherv(distr_vectors%values(1_ik:length,1_ik),&
                       length_mpik,&
                       cel_omp%mpi_type_rk,&
                       tmp_buff,&
                       cel_comm%vtxdist_send_counts(1_ik:cel_omp%num_masters),&
                       cel_comm%vtxdist_displ(1_ik:cel_omp%num_masters),&
                       cel_omp%mpi_type_rk,&
                       0_mpik,&
                       cel_omp%master_comm,ierr)
    end if
  end if

end subroutine cel_sp_mat_distr_collect_global_vec

subroutine cel_sp_mat_distr_vec_isend_start_all(distr_vectors,cel_comm,vtxdist,omp,cel_perf_counter,err_code)
  type(cel_sp_mat_distr_vec_type), intent(inout) :: distr_vectors
  type(cel_comm_type), intent(inout)             :: cel_comm
  integer(kind=ik),dimension(:), intent(in)      :: vtxdist
  integer(kind=ik), intent(out)                  :: err_code
  type(cel_omp_type), intent(in)                 :: omp
  type(cel_perf_counter_type), intent(inout)           :: cel_perf_counter
  integer(kind=mpik)      :: length, root, ierr, ii, dest_mpik, tag_mpik, num_receiver
  integer(kind=mpik)      :: num_sender, source_mpik

  err_code = 0_ik
  if(omp%is_master) then
    root = INT(omp%master_num,kind=mpik)
    length = INT(vtxdist(omp%master_num+2_ik)-vtxdist(omp%master_num+1_ik),kind=mpik)

    !SEND if distributed
    num_receiver = size(cel_comm%send_to)
    if(num_receiver .gt. 1_mpik) then
       if(.not. allocated(cel_comm%send_requests)) allocate(cel_comm%send_requests(num_receiver-1_mpik))
       do ii=2_mpik,num_receiver
        dest_mpik = int(cel_comm%send_to(ii),kind=mpik)
        tag_mpik = root + 3_mpik*dest_mpik
        call MPI_Isend(distr_vectors%values(1_mpik:length,1_mpik),length,omp%mpi_type_rk,&
              dest_mpik,tag_mpik,cel_comm%comm,cel_comm%send_requests(ii-1_mpik),ierr)
        call cel_error_proc(&
               "Error cel_sp_mat_distr_vec_start_all: MPI_Isend!", &
                int(omp%master_num,kind=mpik), int(ierr,kind=mpik), cel_is_in_debug, .FALSE._lk)
!        call MPI_Request_free(cel_comm%send_requests(ii-1_mpik),ierr)
        cel_perf_counter%send_bytes = cel_perf_counter%send_bytes + length*rk
        cel_perf_counter%send_messages_num = cel_perf_counter%send_messages_num + 1_perf_ik
       end do
    end if
   !GET
    num_sender = size(cel_comm%get_from)
    if(num_sender .gt. 1_mpik) then
      if(.not. allocated(cel_comm%get_requests)) allocate(cel_comm%get_requests(num_sender-1_mpik))
      do ii=2_ik, num_sender
        root = int(cel_comm%get_from(ii),kind=mpik)
        length = INT(vtxdist(root+2_ik)-vtxdist(root+1_ik),kind=mpik)
        tag_mpik = root + 3_mpik*int(omp%master_num,kind=mpik)
        call MPI_Irecv(distr_vectors%values(1_mpik:length,ii),length,omp%mpi_type_rk,&
                       root,tag_mpik,cel_comm%comm,cel_comm%get_requests(ii-1_mpik),ierr)
        call cel_error_proc(&
                           "Error cel_sp_mat_distr_vec_start_all: MPI_Irecv!", &
                           int(omp%master_num,kind=mpik), int(ierr,kind=mpik),&
                           cel_is_in_debug, .FALSE._lk)
        cel_perf_counter%get_bytes = cel_perf_counter%get_bytes + length*rk
        cel_perf_counter%get_messages_num = cel_perf_counter%get_messages_num + 1_perf_ik
      end do
    end if
   end if

end subroutine cel_sp_mat_distr_vec_isend_start_all

subroutine cel_sp_mat_distr_vec_isend_idx_start_all(distr_vectors,cel_comm,vtxdist,omp,cel_perf_counter,err_code)
  type(cel_sp_mat_distr_vec_type), intent(inout) :: distr_vectors
  type(cel_comm_type), intent(inout)             :: cel_comm
  integer(kind=ik),dimension(:), intent(in)      :: vtxdist
  integer(kind=ik), intent(out)                  :: err_code
  type(cel_omp_type), intent(in)                 :: omp
  type(cel_perf_counter_type), intent(inout)           :: cel_perf_counter
  integer(kind=mpik)      :: length, root, ierr, ii, dest_mpik, tag_mpik, num_receiver
  integer(kind=ik)      :: first_idx, last_idx
  integer(kind=mpik)      :: num_sender, source_mpik

  err_code = 0_ik
  if(omp%is_master) then
    root = INT(omp%master_num,kind=mpik)
    length = INT(vtxdist(omp%master_num+2_ik)-vtxdist(omp%master_num+1_ik),kind=mpik)

    !SEND if distributed
    num_receiver = size(cel_comm%send_to)
    if(num_receiver .gt. 1_mpik) then
       if(.not. allocated(cel_comm%send_requests)) allocate(cel_comm%send_requests(num_receiver-1_mpik))
       do ii=2_mpik,num_receiver
        dest_mpik = int(cel_comm%send_to(ii),kind=mpik)
        tag_mpik = root + 3_mpik*dest_mpik
        first_idx = cel_comm%send_to_first_idx(ii)
        last_idx = cel_comm%send_to_last_idx(ii)
        length = int(last_idx-first_idx,kind=mpik)+1_mpik
        call MPI_Isend(distr_vectors%values(first_idx:last_idx,1_ik),length,omp%mpi_type_rk,&
              dest_mpik,tag_mpik,cel_comm%comm,cel_comm%send_requests(ii-1_mpik),ierr)
        call cel_error_proc(&
               "Error cel_sp_mat_distr_vec_start_all: MPI_Isend!", &
                int(omp%master_num,kind=mpik), int(ierr,kind=mpik), cel_is_in_debug, .FALSE._lk)
!        call MPI_Request_free(cel_comm%send_requests(ii-1_mpik),ierr)
        cel_perf_counter%send_bytes = cel_perf_counter%send_bytes + int(length*rk,kind=ik)
        cel_perf_counter%send_messages_num = cel_perf_counter%send_messages_num + 1_perf_ik
       end do
    end if
   !GET
    num_sender = size(cel_comm%get_from)
    if(num_sender .gt. 1_mpik) then
      if(.not. allocated(cel_comm%get_requests)) allocate(cel_comm%get_requests(num_sender-1_mpik))
      do ii=2_ik, num_sender
        root = int(cel_comm%get_from(ii),kind=mpik)
        first_idx = cel_comm%get_from_first_idx(ii)
        last_idx = cel_comm%get_from_last_idx(ii)
        length = int(last_idx-first_idx,kind=mpik)+1_mpik
        tag_mpik = root + 3_mpik*int(omp%master_num,kind=mpik)
        call MPI_Irecv(distr_vectors%values(first_idx:last_idx,ii),length,omp%mpi_type_rk,&
                       root,tag_mpik,cel_comm%comm,cel_comm%get_requests(ii-1_mpik),ierr)
        call cel_error_proc(&
                           "Error cel_sp_mat_distr_vec_start_all: MPI_Irecv!", &
                           int(omp%master_num,kind=mpik), int(ierr,kind=mpik),&
                           cel_is_in_debug, .FALSE._lk)
        cel_perf_counter%get_bytes = cel_perf_counter%get_bytes + int(length*rk,kind=ik)
        cel_perf_counter%get_messages_num = cel_perf_counter%get_messages_num + 1_perf_ik
      end do
    end if
   end if

end subroutine cel_sp_mat_distr_vec_isend_idx_start_all

subroutine cel_sp_mat_distr_vec_isend_groups_start_all(distr_vectors,cel_comm,vtxdist,omp,cel_perf_counter,err_code)
  type(cel_sp_mat_distr_vec_type), intent(inout) :: distr_vectors
  type(cel_comm_type), intent(inout)             :: cel_comm
  integer(kind=ik),dimension(:), intent(in)      :: vtxdist
  integer(kind=ik), intent(out)                  :: err_code
  type(cel_omp_type), intent(in)                 :: omp
  type(cel_perf_counter_type), intent(inout)           :: cel_perf_counter
  integer(kind=mpik)      :: length, root, ierr, ii, dest_mpik, tag_mpik, num_receiver
  integer(kind=mpik)      :: num_sender, source_mpik
  integer(kind=mpik)      :: root_transl, dest_transl

  err_code = 0_ik
  if(omp%is_master) then
    root = INT(omp%master_num,kind=mpik)
    root_transl = 0_mpik
    length = INT(vtxdist(omp%master_num+2_ik)-vtxdist(omp%master_num+1_ik),kind=mpik)
   
    !SEND if distributed
    num_receiver = size(cel_comm%send_to)
    if(num_receiver .gt. 1_mpik) then
       if(.not. allocated(cel_comm%send_requests)) &
         allocate(cel_comm%send_requests(num_receiver-1_mpik))
       do ii=2_mpik,num_receiver
        dest_mpik = cel_comm%send_to(ii)
        tag_mpik = root + 3_mpik*dest_mpik
        dest_transl = cel_comm%send_to_transl(ii)
        call MPI_Isend(distr_vectors%values(1_mpik:length,1_mpik),length,omp%mpi_type_rk,&
              dest_transl,tag_mpik,cel_comm%send_comm,cel_comm%send_requests(ii-1_mpik),ierr)
        call cel_error_proc(&
               "Error cel_sp_mat_distr_vec_start_all: MPI_Isend!", &
                int(omp%master_num,kind=mpik), int(ierr,kind=mpik), cel_is_in_debug, .FALSE._lk)
!        call MPI_Request_free(cel_comm%send_requests(ii-1_mpik),ierr)
        cel_perf_counter%send_bytes = cel_perf_counter%send_bytes + length*rk
        cel_perf_counter%send_messages_num = cel_perf_counter%send_messages_num + 1_perf_ik
       end do
    end if
 !GET
    num_sender = size(cel_comm%get_from)
    if(num_sender .gt. 1_mpik) then
      if(.not. allocated(cel_comm%get_requests)) allocate(cel_comm%get_requests(num_sender-1_mpik))
      do ii=2_ik, num_sender
        root = cel_comm%get_from(ii)
        length = INT(vtxdist(root+2_ik)-vtxdist(root+1_ik),kind=mpik)
        tag_mpik = root + 3_mpik*int(omp%master_num,kind=mpik)
        call MPI_Irecv(distr_vectors%values(1_mpik:length,ii),length,omp%mpi_type_rk,&
                       root_transl,tag_mpik,cel_comm%get_comms(ii),cel_comm%get_requests(ii-1_mpik),ierr)
        call cel_error_proc(&
                           "Error cel_sp_mat_distr_vec_start_all: MPI_Irecv!", &
                           int(omp%master_num,kind=mpik), int(ierr,kind=mpik),&
                           cel_is_in_debug, .TRUE._lk)
        cel_perf_counter%get_bytes = cel_perf_counter%get_bytes + length*rk
        cel_perf_counter%get_messages_num = cel_perf_counter%get_messages_num + 1_perf_ik
      end do
    end if
   end if

end subroutine cel_sp_mat_distr_vec_isend_groups_start_all

subroutine cel_sp_mat_distr_vec_isend_comm_buffer_start_all(distr_vectors,cel_comm,vtxdist,omp,cel_perf_counter,err_code)
  type(cel_sp_mat_distr_vec_type), intent(inout) :: distr_vectors
  type(cel_comm_type), intent(inout)             :: cel_comm
  integer(kind=ik),dimension(:), intent(in)      :: vtxdist
  integer(kind=ik), intent(out)                  :: err_code
  type(cel_omp_type), intent(in)                 :: omp
  type(cel_perf_counter_type), intent(inout)           :: cel_perf_counter
  integer(kind=mpik)      :: length, root, ierr, ii, dest_mpik, tag_mpik, num_receiver
  integer(kind=mpik)      :: num_sender, source_mpik
  integer(kind=ik) :: exist_buffer_dim_1,exist_buffer_dim_2,need_buffer_dim_1,need_buffer_dim_2
  integer(kind=ik) :: current_buffer_length, jj, idx, count_bytes, sys_stat
  err_code = 0_ik
  if(omp%is_master) then
    root = INT(omp%master_num,kind=mpik)
    length = INT(vtxdist(omp%master_num+2_ik)-vtxdist(omp%master_num+1_ik),kind=mpik)
    need_buffer_dim_1=size(cel_comm%send_com_buffer_perm,1)
    need_buffer_dim_2=size(cel_comm%send_com_buffer_perm,2)
    if(allocated(distr_vectors%send_comm_buffer)) then
      exist_buffer_dim_1=size(distr_vectors%send_comm_buffer,1)
      exist_buffer_dim_2=size(distr_vectors%send_comm_buffer,2)
    else
      exist_buffer_dim_1=0_ik
      exist_buffer_dim_2=0_ik
    end if
    if(need_buffer_dim_1 .gt. exist_buffer_dim_1 .or. need_buffer_dim_2 .gt. exist_buffer_dim_2) then
      if(allocated(distr_vectors%send_comm_buffer)) deallocate(distr_vectors%send_comm_buffer)
      allocate(distr_vectors%send_comm_buffer(need_buffer_dim_1,need_buffer_dim_2), stat=sys_stat)
        call cel_error("Error cel_sp_mat_distr_vec_isend_comm_buffer_start_all alloc ", sys_stat,&
                        .TRUE._lk,.TRUE._lk)
    end if
    !SEND if distributed
    num_receiver = size(cel_comm%send_to)
    !if(omp%master_num .eq. 0) call print("distr_vectors", distr_vectors%values)
    !if(omp%master_num .eq. 0) write(cel_output_unit,'(A,I0)') "num_receiver:", num_receiver
    if(num_receiver .gt. 1_mpik) then
       if(.not. allocated(cel_comm%send_requests)) allocate(cel_comm%send_requests(num_receiver-1_mpik))
       !copy ins buffer
       count_bytes = 0_ik
       do ii=2_mpik,num_receiver
        dest_mpik = cel_comm%send_to(ii)
        tag_mpik = root + 3_mpik*dest_mpik
        current_buffer_length = cel_comm%send_length_comm_buffer(ii)
        count_bytes = count_bytes + current_buffer_length
        !if(omp%master_num .eq. 0) write(cel_output_unit,'(2(A,I0))') "current_buffer_length:", &
        ! current_buffer_length,"; ii:",ii
        do jj=1,current_buffer_length
          idx = cel_comm%send_com_buffer_perm(jj,ii)
          !if(omp%master_num .eq. 3) write(cel_output_unit,'(2(A,I0))') "jj:", jj,"; idx:",idx
          distr_vectors%send_comm_buffer(jj,ii) = distr_vectors%values(idx,1_ik)
        end do
     !   if(omp%master_num .eq. 0) call print("send_comm_buffer", distr_vectors%send_comm_buffer(:,ii))
        call MPI_Isend(distr_vectors%send_comm_buffer(1_mpik:current_buffer_length,ii),&
              current_buffer_length,omp%mpi_type_rk,&
              dest_mpik,tag_mpik,cel_comm%comm,cel_comm%send_requests(ii-1_mpik),ierr)
        call cel_error_proc(&
               "Error cel_sp_mat_distr_vec_isend_comm_buffer_start_all: MPI_Isend!", &
                int(omp%master_num,kind=mpik), int(ierr,kind=mpik), cel_is_in_debug, .FALSE._lk)
!        call MPI_Request_free(cel_comm%send_requests(ii-1_mpik),ierr)
       end do
       cel_perf_counter%send_bytes = cel_perf_counter%send_bytes + count_bytes*rk
       cel_perf_counter%send_messages_num = cel_perf_counter%send_messages_num + 1_perf_ik
       cel_perf_counter%load = cel_perf_counter%load+count_bytes*(ik+rk)
       cel_perf_counter%store = cel_perf_counter%store+count_bytes*rk
    end if
    !GET
    need_buffer_dim_1=size(cel_comm%get_com_buffer_perm,1)
    need_buffer_dim_2=size(cel_comm%get_com_buffer_perm,2)
    if(allocated(distr_vectors%get_comm_buffer)) then
      exist_buffer_dim_1=size(distr_vectors%get_comm_buffer,1)
      exist_buffer_dim_2=size(distr_vectors%get_comm_buffer,2)
    else
      exist_buffer_dim_1=0_ik
      exist_buffer_dim_2=0_ik
    end if
    if(need_buffer_dim_1 .gt. exist_buffer_dim_1 .or. need_buffer_dim_2 .gt. exist_buffer_dim_2) then
      if(allocated(distr_vectors%get_comm_buffer)) deallocate(distr_vectors%get_comm_buffer)
      allocate(distr_vectors%get_comm_buffer(need_buffer_dim_1,need_buffer_dim_2), stat=sys_stat)
        call cel_error("Error cel_sp_mat_distr_vec_isend_comm_buffer_start_all alloc ", sys_stat,&
                        .TRUE._lk,.TRUE._lk)
    end if
    num_sender = size(cel_comm%get_from)
    if(num_sender .gt. 1_mpik) then
      if(.not. allocated(cel_comm%get_requests)) allocate(cel_comm%get_requests(num_sender-1_mpik))
      distr_vectors%is_comm_buffer_sync = .false._lk
      count_bytes = 0_ik
      do ii=2_ik, num_sender
        root = cel_comm%get_from(ii)
        length = INT(vtxdist(root+2_ik)-vtxdist(root+1_ik),kind=mpik)
        tag_mpik = root + 3_mpik*int(omp%master_num,kind=mpik)
        current_buffer_length = cel_comm%get_length_comm_buffer(ii)
        count_bytes = count_bytes + current_buffer_length
!         if(omp%master_num .eq. 1) write(cel_output_unit,'(2(A,I0))') "get current_buffer_length:", &
!         current_buffer_length,"; ii:",ii
        call MPI_Irecv(distr_vectors%get_comm_buffer(1_mpik:current_buffer_length,ii),&
                 current_buffer_length,omp%mpi_type_rk,&
                 root,tag_mpik,cel_comm%comm,cel_comm%get_requests(ii-1_mpik),ierr)
        call cel_error_proc(&
                           "Error cel_sp_mat_distr_vec_isend_comm_buffer_start_all: MPI_Irecv!", &
                           int(omp%master_num,kind=mpik), int(ierr,kind=mpik),&
                           cel_is_in_debug, .TRUE._lk)
      end do
      cel_perf_counter%get_bytes = cel_perf_counter%get_bytes + count_bytes*rk
      cel_perf_counter%get_messages_num = cel_perf_counter%get_messages_num + 1_perf_ik
    end if
   end if

end subroutine cel_sp_mat_distr_vec_isend_comm_buffer_start_all


subroutine cel_sp_mat_distr_vec_synchronize_comm_buffer_worker(distr_vectors,cel_comm,&
   cel_omp_shared_work,cel_omp,cel_perf_counter,err_code)
  type(cel_sp_mat_distr_vec_type), intent(inout) :: distr_vectors
  type(cel_comm_type), intent(inout)             :: cel_comm
  type(cel_omp_shared_work_type), intent(inout)  :: cel_omp_shared_work
  integer(kind=ik), intent(out)                  :: err_code
  type(cel_omp_type), intent(in)                 :: cel_omp
  type(cel_perf_counter_type), intent(inout)     :: cel_perf_counter
  integer(kind=ik) :: omp_start_idx,omp_end_idx, num_sender,num_threads, ii, jj
  integer(kind=ik) :: current_buffer_length, idx, count_bytes
  
  err_code = 0_ik
  if(cel_omp%is_worker) then
    count_bytes=0_ik
    num_sender = size(cel_comm%get_from)
    call cel_omp_shared_work_distr(cel_omp_shared_work,&
        omp_start_idx,omp_end_idx,num_threads,&
        num_sender,num_sender,&
        cel_omp,cel_is_in_debug,err_code)
    if(cel_omp%worker_num .le. num_threads) then
       do ii=omp_start_idx,omp_end_idx
         current_buffer_length = cel_comm%get_length_comm_buffer(ii)
         count_bytes = count_bytes + current_buffer_length
         do jj=1_ik,current_buffer_length
           idx = cel_comm%get_com_buffer_perm(jj,ii)
           distr_vectors%values(idx,ii)=distr_vectors%get_comm_buffer(jj,ii)
         end do
       end do
       !$omp atomic update
       cel_perf_counter%load = cel_perf_counter%load+&
         count_bytes*ik+count_bytes*rk
       !$omp atomic update
       cel_perf_counter%store = cel_perf_counter%store+&
         count_bytes*rk
    else
      call cel_error(&
          "Error cel_sp_mat_distr_vec_synchronize_comm_buffer_worker:"// &
             "cel_omp_shared_work_distr - num_threads!", &
          num_threads, cel_is_in_debug, .TRUE._lk)
          err_code = num_threads
    end if
  end if
  
end subroutine cel_sp_mat_distr_vec_synchronize_comm_buffer_worker

subroutine cel_sp_mat_distr_vec_synchronize_comm_buffer_master(distr_vectors,cel_comm,&
             cel_omp,cel_perf_counter,err_code)
  type(cel_sp_mat_distr_vec_type), intent(inout) :: distr_vectors
  type(cel_comm_type), intent(inout)             :: cel_comm
  integer(kind=ik), intent(out)                  :: err_code
  type(cel_omp_type), intent(in)                 :: cel_omp
  type(cel_perf_counter_type), intent(inout)     :: cel_perf_counter
  integer(kind=ik) :: num_sender, ii, jj
  integer(kind=ik) :: current_buffer_length, idx, count_bytes

  err_code = 0_ik
  if(cel_omp%is_master) then
    count_bytes = 0_ik
    num_sender = size(cel_comm%get_from)
    do ii=2_ik,num_sender
      current_buffer_length = cel_comm%get_length_comm_buffer(ii)
      count_bytes=count_bytes+current_buffer_length
      do jj=1_ik,current_buffer_length
        idx = cel_comm%get_com_buffer_perm(jj,ii)
        distr_vectors%values(idx,ii)=distr_vectors%get_comm_buffer(jj,ii)
      end do
    end do
    distr_vectors%is_comm_buffer_sync = .true._lk
    cel_perf_counter%load = cel_perf_counter%load+&
       count_bytes*(ik+rk)
    cel_perf_counter%store = cel_perf_counter%store+&
       count_bytes*rk
  end if
end subroutine cel_sp_mat_distr_vec_synchronize_comm_buffer_master

subroutine cel_sp_mat_distr_vec_end_all(cel_comm,omp,err_code)
  type(cel_comm_type), intent(in)                :: cel_comm
  type(cel_omp_type), intent(in)                 :: omp
  integer(kind=ik), intent(out)                  :: err_code
  integer(kind=mpik)      :: root, ierr,  num_req, num_receiver
  integer(kind=ik) :: ii
  err_code = 0_ik
  if(omp%is_master) then
    !num_receiver = size(cel_comm%send_to)-1_mpik
    !if(num_receiver .gt. 0_mpik) then
      !call MPI_Waitall(num_receiver, cel_comm%send_requests(1_mpik:num_receiver), MPI_STATUSES_IGNORE, ierr)
    !end if
    num_req = size(cel_comm%get_from)-1_mpik
    if(num_req .gt. 0_mpik) then
      call MPI_Waitall(num_req, cel_comm%get_requests(1:num_req), MPI_STATUSES_IGNORE, ierr)
    end if
  end if
end subroutine cel_sp_mat_distr_vec_end_all

subroutine cel_sp_mat_distr_vec_ibcast_start_all(distr_vectors,cel_comm,vtxdist,omp,cel_perf_counter,err_code)
  type(cel_sp_mat_distr_vec_type), intent(inout) :: distr_vectors
  type(cel_comm_type), intent(inout)                :: cel_comm
  integer(kind=ik),dimension(:), intent(in)      :: vtxdist
  integer(kind=ik), intent(out)                  :: err_code
  type(cel_omp_type), intent(in)                 :: omp
  type(cel_perf_counter_type), intent(inout)           :: cel_perf_counter
  integer(kind=mpik)      :: length, root, ierr, ii, dest_mpik, tag_mpik, num_receiver
  integer(kind=mpik)      :: num_sender, source_mpik, root_transl

  err_code = 0_ik
  if(omp%is_master) then
    root = INT(omp%master_num,kind=mpik)
    length = INT(vtxdist(omp%master_num+2_ik)-vtxdist(omp%master_num+1_ik),kind=mpik)
    !GET
    num_sender = size(cel_comm%get_from)
    !SEND if distributed
    num_receiver = size(cel_comm%send_to)
    if(num_receiver .gt. 1_mpik) then
      root_transl = 0_mpik
      !proc is a root and bcasts the local part of the vector to other procs, which are in the group bcast_send_comm
      if(.not. allocated(cel_comm%get_requests)) allocate(cel_comm%send_requests(1_mpik))
      call MPI_Ibcast(distr_vectors%values(1_ik:length,1_ik), length, omp%mpi_type_rk,&
      root_transl, cel_comm%send_comm, cel_comm%send_requests(1_mpik), ierr)
      call cel_error_proc( "Error cel_sp_mat_distr_vec_ibcast_start_all: MPI_Ibcast send!", &
                           int(omp%master_num,kind=mpik), int(ierr,kind=mpik),&
                           cel_is_in_debug, .FALSE._lk)
      !call MPI_Request_free(cel_comm%send_requests(1_mpik),ierr)
      cel_perf_counter%send_bytes = cel_perf_counter%send_bytes + length*rk
      cel_perf_counter%send_messages_num = cel_perf_counter%send_messages_num + 1_perf_ik
    end if
    if(num_sender .gt. 1_mpik) then
      !proc is a receiver of the bcast operations
      if(.not. allocated(cel_comm%get_requests)) allocate(cel_comm%get_requests(num_sender-1_mpik))
      !write(*,'(I0,A,I0)') omp%master_num," - num_sender: ",num_sender
      do ii=2_ik, num_sender
        root = cel_comm%get_from(ii)
        length = INT(vtxdist(root+2_ik)-vtxdist(root+1_ik),kind=mpik)
        tag_mpik = root + 3_mpik*int(omp%master_num,kind=mpik)
        root_transl = 0_mpik
        call MPI_Ibcast(distr_vectors%values(1_ik:length,ii), length, omp%mpi_type_rk,&
        root_transl, cel_comm%get_comms(ii), cel_comm%get_requests(ii-1_ik), ierr)
        call cel_error_proc( "Error cel_sp_mat_distr_vec_ibcast_start_all: MPI_Ibcast get!", &
                           int(omp%master_num,kind=mpik), int(ierr,kind=mpik),&
                           cel_is_in_debug, .TRUE._lk)
        cel_perf_counter%get_bytes = cel_perf_counter%get_bytes + length*rk
        cel_perf_counter%get_messages_num = cel_perf_counter%get_messages_num + 1_perf_ik
      end do
    end if
   end if

end subroutine cel_sp_mat_distr_vec_ibcast_start_all

subroutine cel_sp_mat_distr_vec_ibcast_end_all(cel_comm,omp,err_code)
  type(cel_comm_type), intent(in)                :: cel_comm
  type(cel_omp_type), intent(in)                 :: omp
  integer(kind=ik), intent(out)                  :: err_code
  integer(kind=mpik)      :: root, ierr,  num_req, req_count
  integer(kind=ik) :: num_receiver, ii, count
  integer(kind=mpik) :: mpistatus(MPI_STATUS_SIZE)
  logical ::flag
  err_code = 0_ik
  if(omp%is_master) then
    req_count = 0_ik
    num_req = size(cel_comm%get_from)-1_mpik
    if(num_req .gt. 0_mpik) then
      call MPI_Waitall(num_req, cel_comm%get_requests(1:num_req), MPI_STATUSES_IGNORE, ierr)
!      do while(req_count .lt. num_req)
!        !call MPI_Testall(num_req, cel_comm%get_requests(1_mpik:num_req), flag, MPI_STATUSES_IGNORE, ierr)
        !call MPI_Waitall(num_req, cel_comm%get_requests(1_mpik:num_req),  MPI_STATUSES_IGNORE, ierr)
!        req_count = 0_ik
!        do ii=1, num_req
!          call MPI_Test(cel_comm%get_requests(ii),flag,mpistatus, ierr)
!          if(flag) req_count = req_count + 1
!        end do
!      end do
    end if
  end if
!  num_receiver = size(cel_comm%send_to)
!  if(num_receiver .gt. 1_mpik) then
!    call MPI_Waitall(1_mpik, cel_comm%send_requests(1_mpik), MPI_STATUSES_IGNORE, ierr)
!  end if

end subroutine cel_sp_mat_distr_vec_ibcast_end_all


subroutine cel_sp_mat_distr_vec_new_from_procs(distr_vectors,procs_array,vtxdist,cel_omp)
  type(cel_sp_mat_distr_vec_type), intent(inout) :: distr_vectors
  integer(kind=mpik),dimension(:), intent(in)    :: procs_array
  integer(kind=ik),dimension(:), intent(in)      :: vtxdist
  type(cel_omp_type), intent(in)      :: cel_omp
  integer(kind=ik) :: ii, length_procs, length_vtxdist, max_length, tmp_int
  integer(kind=ik) :: num_elements
  if(cel_omp%is_master) then
    length_procs=size(procs_array)
    length_vtxdist=size(vtxdist)
    num_elements = vtxdist(procs_array(1)+2_ik)-vtxdist(procs_array(1)+1_ik)
    max_length = num_elements
    do ii=2, length_procs
      tmp_int=vtxdist(procs_array(ii)+2_ik)-vtxdist(procs_array(ii)+1_ik)
      max_length=max(max_length,tmp_int)
    end do
    !if(cel_omp%is_master .and. cel_omp%master_num .eq. 0) then
     ! write(cel_output_unit,'(2(A,I0))') "length_procs:", length_procs,&
     !  "; max_length:", max_length
    !end if
    call cel_sp_mat_distr_vec_new_master(distr_vectors,length_procs,num_elements,&
     max_length)
  end if
end subroutine cel_sp_mat_distr_vec_new_from_procs
 


end module cel_sp_mat_distr_vec_module
