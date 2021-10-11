!     
! File:   cel_sp_mat_distr_gl_module.f90
! Project CRESTA (see details on https://www.cresta-project.eu) Exascale library
! Send you email about the bugs to
! The High Performance Computing Center Stuttgart (HLRS) of the University of Stuttgart
! Dmitry Khabi,  (email: khabi@hlrs.de)
! Created on Jun 14, 2013
!

module cel_sp_mat_distr_gl_module
use cel_types_module
use cel_input_parameters_module
use cel_omp_module
use cel_sp_mat_module
implicit none

contains

!set matrix vector and communication algorithms in the sparse matrix
subroutine cel_sp_mat_distr_set_mv_parameters(a_sp_mats,&
  cel_omp,input_parameters,err_code)
  
  type(cel_sp_mat_type), dimension(:), allocatable, intent(inout)  :: a_sp_mats
  type(cel_omp_type), intent(in)                                 :: cel_omp
  type(cel_input_parameters_type), intent(in)                    :: input_parameters
  integer(kind=ik), intent(out)                                  :: err_code
  integer(kind=ik) :: ii
  err_code = 0

  if(cel_omp%is_master) then
    select case (input_parameters%mv_algorithm) 
       case(MV_CSR)
        a_sp_mats(1_ik)%mv_algorithm=MV_CSR
      case(MV_JAD)
        a_sp_mats(1_ik)%mv_algorithm=MV_JAD
      case(MV_COO)
        a_sp_mats(1_ik)%mv_algorithm=MV_COO
      case default
        a_sp_mats(1_ik)%mv_algorithm=MV_CSR
    end select

    select case (input_parameters%mv_communicator) 
      case(MV_COMM_ISEND)
        a_sp_mats(1_ik)%mv_communicator=MV_COMM_ISEND
      case(MV_COMM_ISEND_GROUPS)
        a_sp_mats(1_ik)%mv_communicator=MV_COMM_ISEND_GROUPS
      case(MV_COMM_IBCAST)
        a_sp_mats(1_ik)%mv_communicator=MV_COMM_IBCAST
      case(MV_COMM_ISEND_IDX)
        a_sp_mats(1_ik)%mv_communicator=MV_COMM_ISEND_IDX
      case(MV_COMM_ISEND_COMM_BUFFER)
       a_sp_mats(1_ik)%mv_communicator=MV_COMM_ISEND_COMM_BUFFER
      case(MV_COMM_ISEND_IDX_AND_COMM_BUFFER)
        call cel_error("Error mv_communicator not implemented ", 1_ik,&
                        cel_is_in_debug,.TRUE._lk)
      case default
        a_sp_mats(1_ik)%mv_communicator=MV_COMM_ISEND_COMM_BUFFER
    end select

    do ii=2_ik, size(a_sp_mats)
      select case (input_parameters%mv_algorithm_off_diag) 
         case(MV_CSR)
          a_sp_mats(ii)%mv_algorithm=MV_CSR
        case(MV_JAD)
          a_sp_mats(ii)%mv_algorithm=MV_JAD
        case(MV_COO)
          a_sp_mats(ii)%mv_algorithm=MV_COO
        case default
          a_sp_mats(ii)%mv_algorithm=MV_CSR
      end select
    end do
  end if
  !$omp barrier
  
end subroutine cel_sp_mat_distr_set_mv_parameters

end module cel_sp_mat_distr_gl_module
