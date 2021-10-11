!     
! File:   cel_vec_module.f90
! Project CRESTA Exascale library (see details on https://www.cresta-project.eu)
! Send you email about the bugs to
! The High Performance Computing Center Stuttgart (HLRS) of the University of Stuttgart
! Dmitry Khabi,  (email: khabi@hlrs.de)
! Created on May 29, 2013
!

module cel_vec_module
use cel_types_module
use cel_timer_interface_module
use cel_error_module
use cel_base_print_module
use, intrinsic :: ISO_C_BINDING
implicit none

interface cel_convert_c
  module procedure cel_vec_c_to_f
end interface cel_convert_c

interface new
  module procedure cel_vec_new
end interface new

interface del
  module procedure cel_vec_del
end interface del

interface print
  module procedure cel_vec_print
end interface print

type, bind(c) :: cel_c_vec_type
  integer(kind=c_ik) :: length! elements_num number of elements
  type(C_PTR) ::  values! values array of the vec
end type cel_c_vec_type

type cel_vec_type
  real(kind=rk), allocatable, dimension(:) ::  values! values array of the vec
  integer(kind=ik) :: offset_indx
end type cel_vec_type


contains

subroutine cel_vec_c_to_f(c_vec, vec, to_deallocate_c)
  type(cel_c_vec_type), intent(inout) :: c_vec
  type(cel_vec_type), intent(inout) :: vec
  logical(kind=lk), intent(in) :: to_deallocate_c
  real(kind=rk), pointer :: c_real_array(:)
  integer(kind=c_ik) :: ierr
  integer(kind=ik) :: sys_stat

  call c_f_pointer(c_vec%values,c_real_array,[c_vec%length])
  if(allocated(vec%values)) deallocate(vec%values)
  allocate(vec%values(c_vec%length), stat=sys_stat)
   call cel_error("cel_vec_c_to_f allocate vec%values", sys_stat, cel_is_in_debug,.TRUE._lk)
  vec%values = c_real_array
  if(to_deallocate_c) ierr = cel_free_double_c(c_real_array(1))
  
end subroutine cel_vec_c_to_f


subroutine cel_vec_print(string,cel_vec)
  character(len=*) :: string
  type(cel_vec_type), intent(in) :: cel_vec
  
  write(*,'(2a)') trim(string)
  write(*,'(a,i0)') "offset_indx: ", cel_vec%offset_indx
  call print("values",cel_vec%values)

end subroutine cel_vec_print

subroutine cel_vec_new(vec, length)
    type(cel_vec_type), intent(inout) :: vec
    integer(kind=ik) :: length
    integer(kind=ik) :: sys_stat

    call cel_vec_del(vec)
    allocate(vec%values(length), stat=sys_stat)
    call cel_error("cel_vec_new", sys_stat, cel_is_in_debug, .TRUE._lk)
    vec%offset_indx=0_ik
end subroutine cel_vec_new

subroutine cel_vec_del(vec)
    type(cel_vec_type), intent(inout) :: vec

    if( allocated(vec%values) ) then
      deallocate(vec%values)
    end if
end subroutine cel_vec_del

!Computes vec_y = vec_x + alpha * vec_y
subroutine cel_vec_AYPX(vec_y, alpha, vec_x)
  type(cel_vec_type), intent(inout) :: vec_y
  real(kind=rk), intent(in) :: alpha
  type(cel_vec_type), intent(in) :: vec_x
  
  vec_y%values = vec_x%values + alpha*vec_y%values
end subroutine cel_vec_AYPX


end module cel_vec_module
