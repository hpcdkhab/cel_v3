!     
! File:  cel_sp_mat_jad_converter_module.f90
! Project CRESTA (see details on https://cresta-project.eu) Exascale library
! Send you email about the bugs to
! The High Performance Computing Center Stuttgart (HLRS) of the University of Stuttgart
! Modified by Dmitry Khabi,  (email: khabi@hlrs.de)
! Created on Aug 27, 2013
!
module cel_sp_mat_jad_converter_module
!library fbase
use cel_types_module
use cel_input_parameters_module
use cel_omp_module
use cel_comm_module
use cel_perf_module
use cel_timer_interface_module
use cel_omp_shared_work_module
!library fspmt
use cel_vec_module
use cel_sp_mat_module
!external libraries
use MPI
implicit none

contains 

subroutine cel_sp_mat_jad_from_sp_mat_distr(sp_mats,cel_comm,cel_omp_shared_work,vtxdist,cel_omp,&
                                cel_perf_counter,output_on,err_code)
  type(cel_sp_mat_type), dimension(:), allocatable,intent(inout)    :: sp_mats
  type(cel_comm_type), intent(inout)                                :: cel_comm
  type(cel_omp_shared_work_type), intent(inout)                     :: cel_omp_shared_work
  integer(kind=ik),dimension(:), intent(in)                         :: vtxdist
  type(cel_omp_type), intent(in)                                    :: cel_omp
  type(cel_perf_counter_type), intent(inout)                        :: cel_perf_counter
  logical(kind=lk)  , intent(in)                                    :: output_on
  integer(kind=ik), intent(out)                                     :: err_code
  integer(kind=ik) :: ii

   err_code = 0_ik  
  if(cel_omp%is_master) then
    do ii=1_ik, size(sp_mats)
      call cel_sp_mat_jad_from_sp_mat_local(sp_mats(ii),vtxdist,&
                                         cel_perf_counter,err_code)
      call cel_error("Error cel_cg_parameter_new cel_sp_mat_jad_from_sp_mat_distr ", err_code,&
                        cel_is_in_debug,.TRUE._lk)
    end do
  end if

end subroutine cel_sp_mat_jad_from_sp_mat_distr

!! DCSORT computes a sorting permutation for a vector.
!  Author:
!
!    Michael Heroux, Sandra Carney
!    Mathematical Software Research Group
!    Cray Research, Inc.
subroutine cel_sp_mat_jad_dcsort( ival, nn, icnt, index, ilo, ihi )

  integer ( kind = ik ),intent(in)  :: ihi
  integer ( kind = ik ),intent(in)  :: ilo
  integer ( kind = ik ),intent(in)  :: nn
  integer ( kind = ik ),intent(in)  :: ival(nn)
  integer ( kind = ik ),intent(inout) :: icnt(ilo:ihi)
  integer ( kind = ik ),intent(inout) :: index(nn)
  
  integer ( kind = ik ) :: i,j
  integer ( kind = ik ) :: ivali,ivalj

  icnt(ilo:ihi) = 0

  do i = 1, nn
    ivali = ival(i)
    icnt(ivali) = icnt(ivali) + 1
  end do

  do i = ihi-1, ilo, -1
    icnt(i) = icnt(i) + icnt(i+1)
  end do

  do j = nn, 1, -1
    ivalj = ival(j)
    index(icnt(ivalj)) = j
    icnt(ivalj) = icnt(ivalj) - 1_ik
  end do

  return

end subroutine cel_sp_mat_jad_dcsort

! CSRJAD converts Compressed Sparse Row to Jagged Diagonal storage.
! Author: Youcef Saad
! Modified: 19.06.2103
subroutine cel_sp_mat_jad_from_sp_mat_local(sp_mat,vtxdist,&
                                cel_perf_counter,err_code)
  type(cel_sp_mat_type), intent(inout)                              :: sp_mat
  integer(kind=ik),dimension(:), intent(in)                         :: vtxdist
  type(cel_perf_counter_type), intent(inout)                        :: cel_perf_counter
  integer(kind=ik), intent(out)                                     :: err_code
  integer(kind=ik) :: ii, jj, kk, k0, k1
  integer(kind=ik) :: nrow, el_num
  integer ( kind = ik ), allocatable, dimension(:) :: temp_array
  integer ( kind = ik ) :: idiag, ilo, length
  err_code = 0_ik

  if(sp_mat%is_in_csr) then
  
    !
    !  Define initial IPERM and get lengths of each row.
    !  JAO is used a work vector to store tehse lengths.
    !
    call del_jad(sp_mat)
    nrow = sp_mat%rows_num
    el_num = sp_mat%elements_num
    allocate(sp_mat%jad_values(el_num))
    allocate(sp_mat%jad_cols(el_num))
    allocate(sp_mat%jad_perm(nrow))
    allocate(temp_array(nrow))
    idiag = 0
    ilo = nrow
    do jj = 1, nrow
      sp_mat%jad_perm(jj) = jj
      length = sp_mat%xadj(jj+1) - sp_mat%xadj(jj)
      ilo = min ( ilo, length )
      idiag = max ( idiag, length )
      temp_array(jj) = length
    end do
   
    allocate(sp_mat%jad_begin(idiag+1))
    sp_mat%jad_mmax = idiag
    sp_mat%jad_nmax = nrow
   

  !  Call the sorter to get permutation
    call cel_sp_mat_jad_dcsort ( temp_array, nrow ,&
      sp_mat%jad_begin, sp_mat%jad_perm, ilo, idiag )
    

    !sp_mat%jad_begin = 0_ik

    !do kk = 1, nrow
    !  length = temp_array(sp_mat%jad_perm(kk))
    !  sp_mat%jad_begin(1:length) = sp_mat%jad_begin(1:length) + 1_ik
    !end do
   
    k1 = 1_ik
    k0 = k1
    do jj = 1, idiag
      length = sp_mat%jad_begin(jj)
      do kk = 1, length
        ii = sp_mat%xadj(sp_mat%jad_perm(kk)) + jj -1_ik
        sp_mat%jad_values(k1) = sp_mat%values(ii)
        sp_mat%jad_cols(k1) = sp_mat%cols(ii) - sp_mat%col_offset
        k1 = k1 + 1_ik
      end do
      sp_mat%jad_begin(jj) = k0
      k0 = k1
    end do

    sp_mat%jad_begin(idiag+1) = k1
    
    if(sp_mat%jad_nmax .GT. VRL*2_IK) then
      allocate(sp_mat%jad_tmp_values(sp_mat%jad_nmax))
      allocate(sp_mat%jad_tmp_reductions_values(sp_mat%jad_nmax))
      
    else
      allocate(sp_mat%jad_tmp_values(VRL*2_IK))
      allocate(sp_mat%jad_tmp_reductions_values(VRL*2_IK))
    end if

    deallocate(temp_array)
  else
    err_code=1_ik
    allocate(temp_array(0_ik))
  end if
end subroutine cel_sp_mat_jad_from_sp_mat_local

end module cel_sp_mat_jad_converter_module
