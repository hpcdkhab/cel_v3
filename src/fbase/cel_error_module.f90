!     
! File:  cel_error_module.f90
! Project CRESTA (see details on https://cresta-project.eu) Exascale library
! Send you email about the bugs to
! The High Performance Computing Center Stuttgart (HLRS) of the University of Stuttgart
! Dmitry Khabi,  (email: khabi@hlrs.de)
! Created on May 29, 2013
! Autor: Uwe Kuester (HLRS)

module cel_error_module
use cel_types_module

implicit none
contains


subroutine cel_error(string, return_code, output_on, to_stop)
  character(*), intent(in) ::string
  integer(kind=ik), intent(in) ::return_code
  logical(kind=lk), intent(in) ::output_on
  logical(kind=lk), intent(in) ::to_stop
  if(return_code .NE. 0) then
    if(output_on .and. cel_is_in_debug) then
      write(cel_error_unit, "(79('!'))")
      write(cel_error_unit, '(A)') TRIM(string)
      write(cel_error_unit,'(A,I0)') "ERROR::function returned not nil:", return_code
      write(cel_error_unit, "(79('!'))")
      write(cel_error_unit, '(A,L)') 'Exit ...:',to_stop
    end if
    if(to_stop) STOP
  end if
end subroutine cel_error

subroutine cel_warning(string, return_code, output_on)
  character(*), intent(in) ::string
  integer(kind=ik), intent(in) ::return_code
  logical(kind=lk), intent(in) ::output_on
  if(return_code .NE. 0) then
    if(output_on .and. cel_is_in_debug) then
      write(cel_error_unit, "(79('!'))")
      write(cel_error_unit,'(A,I0)') "WARNING:", return_code
      write(cel_error_unit, '(A)') TRIM(string)
      write(cel_error_unit, "(79('!'))")
    end if
  end if
end subroutine cel_warning

subroutine cel_error_proc(string, proc_id, return_code, output_on, to_stop)
  character(*), intent(in) ::string
  integer(kind=mpik), intent(in) ::proc_id
  integer(kind=mpik), intent(in) ::return_code
  logical(kind=lk), intent(in) ::output_on
  logical(kind=lk), intent(in) ::to_stop
  if(return_code .NE. 0) then
    if(output_on .and. cel_is_in_debug) then
      write(cel_error_unit, "(79('!'))")
      write(cel_error_unit, '(I0,A,A)') proc_id,": ",TRIM(string)
      write(cel_error_unit,'(I0,A,A,I0)') proc_id,": ","ERROR::function returned not nil:", return_code
      write(cel_error_unit, '(I0,A,A)')  proc_id,": ",'Exit ...'
    end if
    if(to_stop) STOP
  end if
end subroutine cel_error_proc

subroutine cel_error_petsc(petsc_ierr)

  !-Dummy parameters
  integer(kind = mpik), intent(in) :: petsc_ierr

  if (petsc_ierr .NE. 0) then
    if(cel_is_in_debug) then
     write(cel_error_unit, "(79('!'))")
     write(cel_error_unit, '(A,I0,A)') 'PETSC ERROR :', petsc_ierr, ";"
     write(cel_error_unit, "(79('!'))")
     write(cel_error_unit, '(A)') 'Exit ...'
    end if
     stop
  end if

end subroutine cel_error_petsc

!=========
! Subroutine to evaluate file I/O errors
!=========
SUBROUTINE cel_file_err(filename, io_stat)

  !-Dummy parameters
  CHARACTER(Len = *), INTENT(in) :: filename
  INTEGER(kind = ik), INTENT(in) :: io_stat

  !==Code
  IF (io_stat /= 0) THEN
    if(cel_is_in_debug) then
      WRITE(cel_error_unit, "(79('!'))")
      WRITE(cel_error_unit, '(A,I0,A)') 'ERROR operating on file with io_stat: ', io_stat, ';'
      WRITE(cel_error_unit, '(A)') trim(filename)
      WRITE(cel_error_unit, "(79('!'))")
      WRITE(cel_error_unit, '(A)') 'Exit ...'
    end if
      STOP
  END IF
END SUBROUTINE cel_file_err

end module cel_error_module
