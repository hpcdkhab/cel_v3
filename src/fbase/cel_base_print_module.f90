!     
! File:  cel_base_print_module.f90
! Project CRESTA (see details on https://cresta-project.eu) Exascale library
! Send you email about the bugs to
! The High Performance Computing Center Stuttgart (HLRS) of the University of Stuttgart
! Dmitry Khabi,  (email: khabi@hlrs.de)
! Created on May 29, 2013
! Autor: Uwe Kuester (HLRS)

module cel_base_print_module
use cel_types_module
use cel_string_module
implicit none
integer(kind=ik)  :: line_length=140

interface print
  module procedure cel_base_print_array_1D
  !module procedure cel_base_print_array_1D_mpik
  module procedure cel_base_print_array_2D
  module procedure cel_base_print_array_1D_real
  module procedure cel_base_print_array_2D_real
  module procedure cel_base_print_array_1D_perm
end interface print

contains
!print 1D integer array
!string (in) - label to print 
!array (in) - array to print
subroutine cel_base_print_array_1D(string,array)
  character(len=*)      :: string
  integer(kind=ik),dimension(:)  :: array
  
  integer(kind=ik)  :: nblock
  integer(kind=ik)  :: mm,nn
  integer(kind=ik)  :: size_array
  character(len=7)  :: fmt_int
  character(len=7)  :: fmt_chr
  integer(kind=ik)  :: field_length
  integer(kind=ik)  :: max_array,max_field_number

  size_array=size(array)
  write(*,'(2a)') trim(string)
  write(*,'(a,i0,a)') '',size_array,' elements'

  max_array=maxval(array)
  max_field_number=max(size_array,max_array)

  field_length=INT(log10(real(max_field_number)),kind=ik) + 1_ik + 1_ik
  if(minval(array) < 0 ) field_length=field_length+1

  nblock=line_length/field_length

  fmt_int='(i'//trim(cel_string_of(field_length))//')'
  fmt_chr='(a'//trim(cel_string_of(field_length))//')'



  do nn=1,size_array,nblock
    write(*,*)
    write(*,fmt_chr,advance='no') 'no:'
    do mm=nn,min(size_array,nn+nblock-1)
      write(*,fmt_int,advance='no') mm
    enddo
    write(*,'(a)',advance='yes') ''
    write(*,fmt_chr,advance='no') '   '
    do mm=nn,min(size_array,nn+nblock-1)
      write(*,fmt_int,advance='no') array(mm)
    enddo
    write(*,'(a)',advance='yes') ''
  enddo

  write(*,'(2a)') 'end: ',trim(string)
  write(*,'(a)',advance='yes') ''

end subroutine cel_base_print_array_1D

!print 1D integer array kind=mpik
!string (in) - label to print 
!array (in) - array to print
subroutine cel_base_print_array_1D_mpik(string,array)
  character(len=*)      :: string
  integer(kind=mpik),dimension(:)  :: array
  
  integer(kind=ik)  :: nblock
  integer(kind=ik)  :: mm,nn
  integer(kind=ik)  :: size_array
  character(len=7)  :: fmt_int
  character(len=7)  :: fmt_chr
  integer(kind=ik)  :: field_length
  integer(kind=ik)  :: max_array,max_field_number

  size_array=size(array)
  write(*,'(2a)') trim(string)
  write(*,'(a,i0,a)') '',size_array,' elements'

  max_array=maxval(array)
  max_field_number=max(size_array,max_array)

  field_length=INT(log10(real(max_field_number)),kind=ik) + 1_ik + 1_ik
  if(minval(array) < 0 ) field_length=field_length+1

  nblock=line_length/field_length

  fmt_int='(i'//trim(cel_string_of(field_length))//')'
  fmt_chr='(a'//trim(cel_string_of(field_length))//')'



  do nn=1,size_array,nblock
    write(*,*)
    write(*,fmt_chr,advance='no') 'no:'
    do mm=nn,min(size_array,nn+nblock-1)
      write(*,fmt_int,advance='no') mm
    enddo
    write(*,'(a)',advance='yes') ''
    write(*,fmt_chr,advance='no') '   '
    do mm=nn,min(size_array,nn+nblock-1)
      write(*,fmt_int,advance='no') array(mm)
    enddo
    write(*,'(a)',advance='yes') ''
  enddo

  write(*,'(2a)') 'end: ',trim(string)
  write(*,'(a)',advance='yes') ''

end subroutine cel_base_print_array_1D_mpik

!print 1D integer array
!string (in) - label to print 
!array (in) - array to print
subroutine cel_base_print_array_1D_perm(string,array,array_perm)
  character(len=*)      :: string
  integer(kind=ik),dimension(:)  :: array
  integer(kind=ik),dimension(:)  :: array_perm
  integer(kind=ik)  :: nblock
  integer(kind=ik)  :: mm,nn
  integer(kind=ik)  :: size_array
  character(len=7)  :: fmt_int
  character(len=7)  :: fmt_chr
  integer(kind=ik)  :: field_length
  integer(kind=ik)  :: max_array,max_field_number

  
  size_array=size(array_perm)
  write(*,'(2a)') trim(string)
  write(*,'(a,i0,a)') '',size_array,' elements'

  max_array=maxval(array)
  max_field_number=max(size_array,max_array)

  field_length=INT(log10(real(max_field_number)),kind=ik) + 1_ik + 1_ik
  if(minval(array) < 0 ) field_length=field_length+1

  nblock=line_length/field_length

  fmt_int='(i'//trim(cel_string_of(field_length))//')'
  fmt_chr='(a'//trim(cel_string_of(field_length))//')'



  do nn=1,size_array,nblock
    write(*,*)
    write(*,fmt_chr,advance='no') 'no:'
    do mm=nn,min(size_array,nn+nblock-1)
      write(*,fmt_int,advance='no') mm
    enddo
    write(*,'(a)',advance='yes') ''
    write(*,fmt_chr,advance='no') '   '
    do mm=nn,min(size_array,nn+nblock-1)
      write(*,fmt_int,advance='no') array(array_perm(mm))
    enddo
    write(*,'(a)',advance='yes') ''
  enddo

  write(*,'(2a)') 'end: ',trim(string)
  write(*,'(a)',advance='yes') ''

end subroutine cel_base_print_array_1D_perm

!print 2D integer array
!string (in) - label to print 
!array (in) - array to print
!mask (in) - array to print
subroutine cel_base_print_array_2D(string,array,mask)
  character(len=*)      :: string
  integer(kind=ik),dimension(:,:):: array
  logical(kind=lk),dimension(:),optional  :: mask

  integer(kind=ik)               :: nblock
  integer(kind=ik)               :: mm,nn,ll
  integer(kind=ik)               :: size_array_1,size_array_2
  character(len=7)      :: fmt_int
  character(len=7)      :: fmt_chr
  integer(kind=ik)              :: field_length
  integer(kind=ik)               :: max_array,max_field_number

  size_array_1=size(array,1)
  size_array_2=size(array,2)
  write(*,'(2a)') trim(string)
  write(*,'(2(a,i0),a)') 'array with ',size_array_1,'x',size_array_2,' elements'

  if(size_array_1<=0 .or. size_array_2 <=0) goto 1000

  max_array=maxval(array)
  max_field_number=max(size_array_1,max_array)

  field_length=INT(log10(real(max_field_number)),kind=ik) + 1_ik + 1_ik
  if(minval(array) < 0 ) field_length=field_length+1

  nblock=line_length/field_length

  fmt_int='(i'//trim(cel_string_of(field_length))//')'
  fmt_chr='(a'//trim(cel_string_of(field_length))//')'

  !!print*,' field_length=',field_length
  !!print*,' fmt_int=',fmt_int
  !!print*,' fmt_chr=',fmt_chr
  !!print*,' nblock=',nblock

  if(.not. present(mask)) then
  do nn=1,size_array_1,nblock
    write(*,*)
    write(*,fmt_chr,advance='no') 'no:'
    do mm=nn,min(size_array_1,nn+nblock-1)
      write(*,fmt_int,advance='no') mm
    enddo
      write(*,'(a)',advance='yes') ''
    do ll=1,size_array_2
      write(*,fmt_chr,advance='no') '   '
      do mm=nn,min(size_array_1,nn+nblock-1)
        write(*,fmt_int,advance='no') array(mm,ll)
      enddo
      write(*,'(a)',advance='yes') ''
    enddo
    write(*,'(a)',advance='yes') ''
  enddo

  else
    do nn=1,size_array_1,nblock
      write(*,fmt_chr,advance='no') 'no:'
      do mm=nn,min(size_array_1,nn+nblock-1)
        if(mask(mm)) write(*,fmt_int,advance='no') mm
      enddo
      write(*,'(a)',advance='yes') ''
      do ll=1,size_array_2
        write(*,fmt_chr,advance='no') '   '
        do mm=nn,min(size_array_1,nn+nblock-1)
          if(mask(mm))    write(*,fmt_int,advance='no') array(mm,ll)
        enddo
        write(*,'(a)',advance='yes') ''
      enddo
      write(*,'(a)',advance='yes') ''
    enddo
  endif

  1000 continue
       write(*,'(2a)') 'end: ',trim(string)
       write(*,'(a)',advance='yes') ''

end subroutine cel_base_print_array_2D

!print 1D real array
!string (in) - label to print 
!array (in) - array to print
subroutine cel_base_print_array_1D_real(string,array)
  character(len=*)      :: string
  real(kind=rk),dimension(:)  :: array
  
  integer(kind=ik)               :: nblock
  integer(kind=ik)               :: mm,nn
  integer(kind=ik)               :: size_array
  character(len=7)      :: fmt_int
  character(len=10)      :: fmt_fp
  character(len=7)      :: fmt_chr
  integer(kind=ik)               :: field_length

  size_array=size(array)
  write(*,'(2a)') trim(string)
  write(*,'(a,i0,a)') '',size_array,' elements'
  write(*,'(a,a,2(e10.3,a))') trim(string),': values in [',minval(array),',',maxval(array),']'

  field_length=10
  if(minval(array) < 0 ) field_length=field_length+1

  nblock=line_length/field_length

  fmt_int='(i'//trim(cel_string_of(field_length))//')'
  fmt_fp='(e'//trim(cel_string_of(field_length))//'.3)'
  fmt_chr='(a'//trim(cel_string_of(field_length))//')'

  !!print*,' fmt_int=',fmt_int
  !!print*,' fmt_fp=',fmt_fp
  !!print*,' fmt_chr=',fmt_chr

  do nn=1,size_array,nblock
    write(*,*)
    write(*,'(a4)',advance='no') 'no:'
    do mm=nn,min(size_array,nn+nblock-1)
      write(*,fmt_int,advance='no') mm
    enddo
    write(*,'(a)',advance='yes') ''
    write(*,'(a4)',advance='no') '   '
    do mm=nn,min(size_array,nn+nblock-1)
      if(abs(array(mm)) <= cel_real_small) then
        write(*,fmt_chr,advance='no') '-'
      else
        write(*,fmt_fp,advance='no') array(mm)
      endif
    enddo
    write(*,'(a)',advance='yes') ''
  enddo
  write(*,'(2a)') 'end: ',trim(string)
  write(*,'(a)',advance='yes') ''

end subroutine cel_base_print_array_1D_real


!print 2D real array
!string (in) - label to print 
!array (in) - array to print
!mask (in) - array to print
subroutine cel_base_print_array_2D_real(string,array,mask)
  character(len=*)      :: string
  real(kind=rk),dimension(:,:):: array
  logical(kind=lk),dimension(:),optional  :: mask

  integer(kind=ik)               :: nblock
  integer(kind=ik)               :: mm,nn,ll
  integer(kind=ik)               :: size_array_1,size_array_2
  character(len=7)      :: fmt_int
  character(len=7)      :: fmt_chr
  character(len=7)      :: fmt_fp
  integer(kind=ik)               :: field_length
  integer(kind=ik)               :: max_array,max_field_number

  size_array_1=size(array,1)
  size_array_2=size(array,2)
  write(*,'(2a)') trim(string)
  write(*,'(2(a,i0),a)') 'array with ',size_array_1,'x',size_array_2,' elements'

  if(size_array_1<=0 .or. size_array_2 <=0) goto 1000

  max_array=INT(maxval(array),kind=ik)
  max_field_number=max(size_array_1,max_array)

  field_length=INT(log10(real(max_field_number)),kind=ik) + 1_ik + 1_ik
  field_length=10
  if(minval(array) < 0 ) field_length=field_length+1

  nblock=line_length/field_length

  fmt_int='(i'//trim(cel_string_of(field_length))//')'
  fmt_chr='(a'//trim(cel_string_of(field_length))//')'
  fmt_fp='(e'//trim(cel_string_of(field_length))//'.3)'

  !!print*,' field_length=',field_length
  !!print*,' fmt_int=',fmt_int
  !!print*,' fmt_chr=',fmt_chr
  !!print*,' fmt_fp=',fmt_fp
  !!print*,' nblock=',nblock

  if(.not. present(mask)) then
    do nn=1,size_array_2,nblock
      write(*,*)
      write(*,fmt_chr,advance='no') 'no:'
      do mm=nn,min(size_array_2,nn+nblock-1)
        write(*,fmt_int,advance='no') mm
      enddo
      write(*,'(a)',advance='yes') ''
      do ll=1,size_array_1
        !!write(*,fmt_chr,advance='no') '   '
        write(*,fmt_int,advance='no') ll
        do mm=nn,min(size_array_2,nn+nblock-1)
          !!  write(*,fmt_fp,advance='no') array(ll,mm)
          if(abs(array(ll,mm)) <= cel_real_small) then
            write(*,fmt_chr,advance='no') '-    '
          else
            write(*,fmt_fp,advance='no') array(ll,mm)
          endif
        enddo
        write(*,'(a)',advance='yes') ''
      enddo
      write(*,'(a)',advance='yes') ''
    enddo
  else
    do nn=1,size_array_2,nblock
      write(*,fmt_chr,advance='no') 'no:'
      do mm=nn,min(size_array_2,nn+nblock-1)
        if(mask(mm)) write(*,fmt_int,advance='no') mm
      enddo
      write(*,'(a)',advance='yes') ''
      do ll=1,size_array_1
        write(*,fmt_chr,advance='no') '   '
        do mm=nn,min(size_array_2,nn+nblock-1)
          if(mask(mm))    write(*,fmt_int,advance='no') array(mm,ll)
        enddo
      write(*,'(a)',advance='yes') ''
      enddo
      write(*,'(a)',advance='yes') ''
    enddo
  endif

  1000 continue
       write(*,'(2a)') 'end: ',trim(string)
       write(*,'(a)',advance='yes') ''

end subroutine cel_base_print_array_2D_real

end module cel_base_print_module
