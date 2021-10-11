!     
! File:   cel_domain_3d_module.f90
! Project CRESTA (see details on https://cresta-project.eu) Exascale library
! Send you email about the bugs to
! The High Performance Computing Center Stuttgart (HLRS) of the University of Stuttgart
! Dmitry Khabi,  (email: khabi@hlrs.de)
! Created on May 29, 2013
!

module cel_domain_3d_module
use cel_types_module
use, intrinsic :: ISO_C_BINDING
implicit none

! 3D description of regular grid on domain [0;1]x[0;1]x[0;1] for c intefrace
type, bind(c) :: cel_c_domain_3D_type
  integer(kind=c_ik) :: block_indx ! block_indx index of sub-domain
  integer(kind=c_ik) ::  nx! nx number of nodes in X-direction
  integer(kind=c_ik) ::  ny! ny number of nodes in Y-direction
  integer(kind=c_ik) ::  nz! nz number of nodes in Z-direction
  integer(kind=c_ik) ::  dx! dx number of sub-domains in X-direction
  integer(kind=c_ik) ::  dy! dy number of sub-domains in Y-direction
  integer(kind=c_ik) ::  dz! dz number of sub-domains in Z-direction
  integer(kind=c_ik) ::  blocks! blocks number blocks (sub-domains) in 3D domain
  integer(kind=c_ik) ::  bx! bx number nodes in X-direction in sub-domain
  integer(kind=c_ik) ::  by! by number nodes in Y-direction in sub-domain
  integer(kind=c_ik) ::  bz! bz number nodes in Z-direction in sub-domain
  integer(kind=c_ik) ::  bxyz! number nodes in sub-domain
end type cel_c_domain_3D_type

! 3D description of regular grid on domain [0;1]x[0;1]x[0;1]
type cel_domain_3D_type
  integer(kind=ik) :: block_indx ! block_indx index of sub-domain
  integer(kind=ik) ::  nx! nx number of nodes in X-direction
  integer(kind=ik) ::  ny! ny number of nodes in Y-direction
  integer(kind=ik) ::  nz! nz number of nodes in Z-direction
  integer(kind=ik) ::  dx! dx number of sub-domains in X-direction
  integer(kind=ik) ::  dy! dy number of sub-domains in Y-direction
  integer(kind=ik) ::  dz! dz number of sub-domains in Z-direction
  integer(kind=ik) ::  blocks! blocks number blocks (sub-domains) in 3D domain
  integer(kind=ik) ::  bx! bx number nodes in X-direction in sub-domain
  integer(kind=ik) ::  by! by number nodes in Y-direction in sub-domain
  integer(kind=ik) ::  bz! bz number nodes in Z-direction in sub-domain
  integer(kind=ik) ::  bxyz! number nodes in sub-domain
end type cel_domain_3D_type

interface cel_convert_c
  module procedure cel_domain_3D_c_to_f
end interface cel_convert_c

contains

subroutine cel_domain_3D_c_to_f(cel_c_domain_3D, cel_domain_3D, to_deallocate_c)
  type(cel_c_domain_3D_type), intent(inout) :: cel_c_domain_3D
  type(cel_domain_3D_type), intent(inout) :: cel_domain_3D
  logical(kind=lk), intent(in) :: to_deallocate_c
  
  cel_domain_3D%block_indx = cel_c_domain_3D%block_indx
  cel_domain_3D%nx = cel_c_domain_3D%nx
  cel_domain_3D%ny = cel_c_domain_3D%ny
  cel_domain_3D%nz = cel_c_domain_3D%nz
  cel_domain_3D%dx = cel_c_domain_3D%dx
  cel_domain_3D%dy = cel_c_domain_3D%dy
  cel_domain_3D%dz = cel_c_domain_3D%dz
  cel_domain_3D%blocks = cel_c_domain_3D%blocks
  cel_domain_3D%bx = cel_c_domain_3D%bx
  cel_domain_3D%by = cel_c_domain_3D%by
  cel_domain_3D%bz = cel_c_domain_3D%bz
  cel_domain_3D%bxyz = cel_c_domain_3D%bxyz
  
end subroutine cel_domain_3D_c_to_f

end module cel_domain_3d_module
