!     
! File:  cel_string_module.f90
! Project CRESTA (see details on https://cresta-project.eu) Exascale library
! Send you email about the bugs to
! The High Performance Computing Center Stuttgart (HLRS) of the University of Stuttgart
! Dmitry Khabi,  (email: khabi@hlrs.de)
! Created on May 29, 2013


module cel_string_module
use cel_types_module
implicit none

contains
!Convert Integer to string
!integer(int) - integer to convert
!return string of integer
function cel_string_of(integer)
integer(kind=ik)  :: integer ! delivers an integer as string of correct length

character(len=20)               :: cel_string_of
character(len=64)              :: temp

write(temp,*) integer
cel_string_of=trim(adjustl(temp))

end function cel_string_of

end module cel_string_module
