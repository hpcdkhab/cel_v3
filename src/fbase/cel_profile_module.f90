!     
! File:   cel_perf_module.f90
! Project CRESTA (see details on https://www.cresta-project.eu) Exascale library
! Send you email about the bugs to
! The High Performance Computing Center Stuttgart (HLRS) of the University of Stuttgart
! Autor: Uwe Küster
! Modified: Dmitry Khabi,  (email: khabi@hlrs.de)
! Created on Jul 10, 2013
!

module  cel_profile_module
!!use constants_module
use cel_interaction_module
use cel_omp_module
use cel_types_module
use cel_cpu_param_module
implicit none



type cel_profile_type
   integer(kind=ik)        :: id
   integer(kind=itk)        :: start_time_sec
   integer(kind=itk)        :: start_time_nanosec
   integer(kind=itk)        :: end_time_sec
   integer(kind=itk)        :: end_time_nanosec
end type cel_profile_type


contains
!***************************************************************************************************************/

subroutine cel_profile_write(cel_profiles,filename)
  type(cel_profile_type),dimension(:),intent(in)  :: cel_profiles
  character(len=*),intent(in)             :: filename
  logical(kind=lk) :: show_header
  integer :: unit,iostat
  integer(kind=ik) :: ii
  
  show_header = .false._lk
  unit = 3233
  open(unit=unit,iostat=iostat,file=filename,status='NEW',action='READWRITE')
  write(cel_output_unit,'(A,A,A,I0)')'check the profile file:"',&
        trim(filename),'" exists (>0-old file exists): ', iostat
  if(iostat .eq. 0) then
    show_header = .true._lk
  else
    open(unit=unit,iostat=iostat,file=filename,status='OLD',action='READWRITE',position='append')
    write(cel_output_unit,'(A,A,A,I0)')'reopen the profile file:"',&
        trim(filename),'" status (0-ok): ', iostat
    show_header = .false._lk
  end if
  if (iostat .gt. 0) then
    print*,'cel_profile_write: error, could not open the file:'
    print*,filename
  else
    if(show_header) then
      write(unit,'(A)',advance='no') "#id;"!1
      write(unit,'(A)',advance='no') "start sec;"!2
      write(unit,'(A)',advance='no') "start nanosec;"!3
      write(unit,'(A)',advance='no') "start sec;"!4
      write(unit,'(A)',advance='yes') "end nanosec;"!5
    end if
    do ii=1, size(cel_profiles)
      write(unit,'(I0,A)',advance='no') cel_profiles(ii)%id,";"
      write(unit,'(I0,A)',advance='no') cel_profiles(ii)%start_time_sec,";"
      write(unit,'(I0,A)',advance='no') cel_profiles(ii)%start_time_nanosec,";"
      write(unit,'(I0,A)',advance='no') cel_profiles(ii)%end_time_sec,";"
      write(unit,'(I0,A)',advance='yes') cel_profiles(ii)%end_time_nanosec,";"
    end do

    close(unit=unit)
  end if  

  
end subroutine cel_profile_write

!***********************************************************************
end module  cel_profile_module
!***********************************************************************

