!     
! File:   cel_interaction_module.f90
! Project CRESTA (see details on https://cresta-project.eu) Exascale library
! Send you email about the bugs to
! The High Performance Computing Center Stuttgart (HLRS) of the University of Stuttgart
! Autor: Uwe Küster
! Modified: Dmitry Khabi,  (email: khabi@hlrs.de)
! Created on Jul 10, 2013
!

module cel_interaction_module

use cel_types_module

implicit none
!*************************************************************************************
     integer(kind=ik),parameter :: log_file_unit=99
     integer(kind=ik)           :: err_unit   

     interface input
          module procedure input_logical
          module procedure input_integer
          module procedure input_integer_array
          module procedure input_real
          module procedure input_char_string
     end interface

     interface print
          module procedure print_short_integer
          module procedure print_short_real
          module procedure print_short_complex
          module procedure print_short_string
     end interface

     interface print_err
          module procedure print_err_logical
          module procedure print_err_integer
          module procedure print_err_integer_array
          module procedure print_err_real4
          module procedure print_err_real8
          module procedure print_err_complex
          module procedure print_err_string
          module procedure print_err_string_array
     end interface


     interface string_of
          module procedure string_of_logical
          module procedure string_of_integer
          module procedure string_of_real
          module procedure string_of_double
     end interface

contains
!******************************************************************************
function parameter_of_keyword(line,keyword,content) result(result_value)

! takes a set of lines and searches for keyword
! if found, result_value will be the respective line number
! content is the rest of the line without the keyword
! if there is no rest content will be the next line;
! in this case result_value will be the number of this line.
! result_value is the line number where content has been found.
! If keyword is appearing more than once, the first
! occurence will be taken.

integer(kind=ik)   :: number_of_arguments
integer(kind=ik)   :: n_arg
integer(kind=ik)   :: arg_no
integer(kind=ik)   :: ind,ind_no
integer(kind=ik)   :: result_value
character(len=*),dimension(:)  :: line
character(len=80)   :: string
character(len=*),optional    :: content
character(len=*)    :: keyword

number_of_arguments=size(line)

arg_no=0_ik
ind_no=0_ik
do n_arg=1,number_of_arguments
   string=adjustl(line(n_arg))
   ind=index(string,trim(keyword))
   if( ind > 0 ) then
     arg_no=n_arg
     ind_no=ind
     exit
   endif
enddo
if(present(content)) then
   if(arg_no > 1 ) then
         content=string(ind_no+len_trim(keyword):)
      if(len_trim(content) == 0 ) then
         if( (arg_no+1) <= number_of_arguments ) then
            arg_no=arg_no+1
            content=line(arg_no)
         endif
      endif
   endif
endif

result_value=arg_no

if(arg_no > 0 ) then
!!print'(a,a,a,i0,a,a,a)','parameter_of_keyword: for keyword "',trim(keyword),'" found parameter ',arg_no,' "',trim(content),'"'
print'(a,a,a,i5,a,a,a)','parameter_of_keyword: for keyword "',trim(keyword),'" found parameter ',arg_no,' "',trim(content),'"'
else
print'(a,a,a)','parameter_of_keyword: no parameter found for keyword "',trim(keyword),'"'
endif

end function parameter_of_keyword
!******************************************************************************
function allocation_status(string,status,resume) result(result_value)
character(len=*) :: string
integer(kind=ik)          :: status
logical(kind=lk)          :: result_value
logical(kind=lk),optional :: resume
logical(kind=lk)          :: loc_resume


loc_resume=.false.
if(present(resume)) loc_resume=resume

  result_value=.false.
if(status /= 0 ) then
  print*,'allocation for ',trim(string),' not successful'
  print*,'status=',status
  if(.not. loc_resume) stop
  result_value=.true.
endif

end function allocation_status
!******************************************************************************
function are_identical_tokens(string1,string2) result(result_value)
character(len=*)  :: string1
character(len=*)  :: string2
integer(kind=ik)           :: ind
integer(kind=ik)           :: len1
integer(kind=ik)           :: len2
logical(kind=lk)           :: result_value

ind=index(trim(adjustl(string1)),trim(adjustl(string2)))

len1=len_trim(adjustl(string1))
len2=len_trim(adjustl(string2))

result_value = ind>0 .and. len1==len2

end function are_identical_tokens
!******************************************************************************
function upper_to_lower(string) result(result_string)
character(len=*)    :: string
!PGI!character(len=len(string))    :: result_string
character(len=str_length)    :: result_string
integer(kind=ik)                 :: mm

result_string=''
do mm=1,len(string)
   if( 'A' <= string(mm:mm) .and. string(mm:mm) <= 'Z' ) then
      result_string(mm:mm)=char(iachar(string(mm:mm))+32)
   else
      result_string(mm:mm)=string(mm:mm)
   endif
enddo

end function upper_to_lower
!******************************************************************************
function index_array(string,substring) result(result_array)
! delivers array of all indices in string where substring begins
! if substring is not present, result_array has size 0!

! vielleicht aendern in: result_array ist um 1 groesser als die Anzahl und dort 0!
! entspricht dann eher dem Verhalten von index()!
character(len=*)                    :: string
character(len=*)                    :: substring
integer(kind=ik),dimension(len_trim(string)) :: ind
integer(kind=ik),dimension(:),pointer        :: result_array
integer(kind=ik)                             :: indx
integer(kind=ik)                             :: nmax
integer(kind=ik)                             :: length
integer(kind=ik)                             :: nn


indx=index(string,substring)
length=len_trim(string)

      nmax=0
if(indx > 0 ) then
      nmax=nmax+1
      ind(nmax)=indx

   do nn=2,length
   indx=index(string(min(ind(nn-1)+1,length):),substring)
   if(indx > 0 ) then
      nmax=nmax+1
      ind(nmax)=indx
   else
      exit
   endif
   enddo
endif


allocate(result_array(nmax))

result_array=ind(1:nmax)

end function index_array
!******************************************************************************
function space_for_integer(longest_integer) result(result_value)
! calculates the amount for the decimal representation of longest_integer
! to be used for format definition 
integer(kind=ik)              :: longest_integer
integer(kind=ik)              :: result_value
character(len=60)    :: string


write(string,'(i60)') longest_integer
result_value=len_trim(adjustl(string))

end function space_for_integer
!******************************************************************************
function first_larger_decimal(value) result(result_value)
real(kind=rk)     :: value
real(kind=rk)     :: result_value
real(kind=rk)     :: abs_value
real(kind=rk)     :: signum_value
real(kind=rk)     :: dec_exponent
real(kind=rk)     :: factor_10
real(kind=rk)     :: mantissa
real(kind=rk)     :: int_dec_exponent
real(kind=rk)     :: ref_value
real(kind=rk)     :: adjust_array(8)
integer(kind=ik)             :: ii

! the function gives back the first value adjust_array(?)*10**? >=value
! in a way that decimals of value are respected
! example  -0.41 --> -0.40
! example   4.1 -->  5.0

adjust_array=(/1._real_kind,1.5_real_kind,2._real_kind,3._real_kind,4._real_kind,5._real_kind,8._real_kind,10._real_kind/)

abs_value=abs(value)
signum_value=sign(1._real_kind,value)

if( abs_value > 0. ) then
   dec_exponent=log(abs_value)/log(10._real_kind)
   
   int_dec_exponent=floor(dec_exponent)
   
   factor_10=10._real_kind**int_dec_exponent
   
   ! 1. =< mantissa < 10.
   mantissa=signum_value*abs_value/factor_10
   
   
   if(signum_value > 0._real_kind) then
      do ii=1,size(adjust_array)
      ref_value=signum_value*adjust_array(ii)
      if(ref_value >=mantissa) exit
      enddo
   else
      do ii=size(adjust_array),1,-1
      ref_value=signum_value*adjust_array(ii)
      if(ref_value >=mantissa) exit
      enddo
   endif
   
   result_value=ref_value*factor_10
   
else
   result_value=0._real_kind
endif

!!print*,'value=',value,' dec_exponent=',dec_exponent,' factor_10=',factor_10,' mantissa=',mantissa,' result_value=',result_value
!!print*,' signum_value=',signum_value
!!print*,' ref_value=',ref_value
!!print*,' result_value=',result_value

end function first_larger_decimal
!******************************************************************************
subroutine mark_string(position)
integer(kind=ik),dimension(:)                :: position
integer(kind=ik)                             :: kk,jj
integer(kind=ik)                             :: last  
integer(kind=ik)                             :: length
character(len=256)                  :: string

length=size(position)
if(length > 122-96) then
   print*,'mark_string: error'
   print*,'size(position)=',length
   print*,'is too large'
   print*,'It must not be larger than 122-96'
   call stop('in mark_string')
endif


string=' '
last=0_ik
do jj=1,length
   last=position(jj)
   string(last:last)='^'
enddo
write(*,'(a)') string(1:last)  

do jj=1,length
   last=position(jj)
   string(last:last)='|'
enddo
write(*,'(a)') string(1:last)  

do kk=1,last  
                     string(kk:kk)='.'
   if(mod(kk, 5)==0) string(kk:kk)='5'
   if(mod(kk,10)==0) string(kk:kk)='0'
enddo
!!write(*,'(a)') string(1:last)  
!!
!!string=' '
do jj=1,length
   last=position(jj)
   string(last:last)=char(96+jj)    ! 96+jj<=122
enddo
write(*,'(a)') string(1:last)  

end subroutine mark_string
!******************************************************************************
function string_of_integer(integer,format) result(result_string)
integer(kind=ik)                         :: integer
character(len=64)       :: result_string
character(len=*),optional       :: format

if(present(format)) then
   write(result_string,format) integer
else
   write(result_string,'(i64)') integer
endif

result_string=adjustl(result_string)

end function string_of_integer
!******************************************************************************
function string_of_real(real,format) result(result_string)
real(kind=4)    :: real
character(len=64)       :: result_string
character(len=*),optional       :: format

if(present(format)) then
   write(result_string,format) real
else
   write(result_string,'(e20.4)') real
endif

result_string=adjustl(result_string)

end function string_of_real
!******************************************************************************
function string_of_double(double,format) result(result_string)
real(kind=8)    :: double
character(len=64)       :: result_string
character(len=*),optional       :: format

if(present(format)) then
   write(result_string,format) double
else
   write(result_string,'(e20.4)') double
endif

result_string=adjustl(result_string)

end function string_of_double
!******************************************************************************
function string_of_logical(logical,format) result(result_string)
logical(kind=lk)                         :: logical
character(len=1)                :: result_string
character(len=*),optional       :: format

if(present(format)) then
   write(result_string,format) logical
else
   write(result_string,'(l1)') logical
endif

result_string=adjustl(result_string)

end function string_of_logical
!******************************************************************************
subroutine beep(number)
!! beeps number of times. If number not present only one time
integer(kind=ik),optional            :: number         
integer(kind=ik)                     :: no
integer(kind=ik)                     :: nn

if(present(number) ) then
   no=number
else
   no=1
endif

do nn=1,no
   write(*,*) char(7)
enddo

end subroutine beep
!******************************************************************************
!******************************************************************************
subroutine halt(string)
character(len=*),optional   :: string
character(len=str_length)      :: local_string

local_string=''
if(present(string) ) then
   local_string=adjustl(string)
else
   local_string='Halt: '
endif
   write(*,*) 'Halt: ',trim(local_string)
   write(0,*) 'Halt: ',trim(local_string)

write(*,*) 'to continue give <cr> '
write(0,*) 'to continue give <cr> '
read(*,*) 


  write(log_file_unit,*) '        ',' :: halt ',trim(local_string)


end subroutine halt
!******************************************************************************
!******************************************************************************
subroutine close(string)
character(len=*)            :: string
integer(kind=ik)                     :: iunit

iunit=give_unit_of(trim(string))

print*,trim(string),' with unit number ',iunit,' is closed'

close(iunit)

end subroutine close
!******************************************************************************
!******************************************************************************
subroutine open_logfile(string)
character(len=*),optional                   :: string
character(len=str_length)      ::filename 
if(present(string) ) then
   filename=trim(adjustl(string))
else
   filename='log_file'   
endif
open(unit=log_file_unit,file=trim(filename))

end subroutine open_logfile
!******************************************************************************
!******************************************************************************
function ask_for(string) result(result_logical)
character(len=*),optional   :: string
character(len=str_length)   :: local_string
logical(kind=lk)                     :: result_logical
character(len=1)            :: chara

  write(*,*)
if(present(string) ) then
   local_string=adjustl(string)
else
   local_string='please give logical(kind=lk): '
endif
   write(*,*) trim(local_string),' (t f)'

  read(*,*) chara

  result_logical=.false.
  if(chara(1:1) == 't') result_logical=.true.
  if(chara(1:1) == 'y') result_logical=.true.
  if(chara(1:1) == 'j') result_logical=.true.
  if(chara(1:1) == 'T') result_logical=.true.
  if(chara(1:1) == 'Y') result_logical=.true.
  if(chara(1:1) == 'J') result_logical=.true.

  print*,trim(local_string),'=',result_logical
  write(*,*)


  write(log_file_unit,*) result_logical,' :: ',trim(local_string)
  write(0,*) result_logical,' :: ',trim(local_string)

end function ask_for
!******************************************************************************
!******************************************************************************
subroutine input_logical(string,logical)
character(len=*),optional   :: string
character(len=str_length)      :: local_string
logical(kind=lk)                     :: logical
character(len=1)            :: chara

  write(*,*)
if(present(string) ) then
   local_string=adjustl(string)
else
   local_string= 'please give logical(kind=lk): '
endif
   write(*,*) trim(local_string),' (t f)'

  read(*,*) chara

  logical=.false.
  if(chara(1:1) == 't') logical=.true.
  if(chara(1:1) == 'y') logical=.true.
  if(chara(1:1) == 'j') logical=.true.
  if(chara(1:1) == 'T') logical=.true.
  if(chara(1:1) == 'Y') logical=.true.
  if(chara(1:1) == 'J') logical=.true.

  print*,trim(local_string),'=',logical
  write(*,*)


  write(log_file_unit,*) logical,' :: ',trim(local_string)
  write(0,*) logical,' :: ',trim(local_string)

end subroutine input_logical
!******************************************************************************
!******************************************************************************
subroutine input_integer(string,integer)
character(len=*),optional   :: string
character(len=str_length)      :: local_string
integer(kind=ik)                     :: integer

  write(*,*)
if(present(string) ) then
   local_string=adjustl(string)
else
   local_string= 'please give integer(kind=ik): '
endif
   write(*,*) trim(local_string)

  read(*,*) integer
  print*,trim(local_string),'=',integer
  write(*,*)


  write(log_file_unit,*) integer,' :: ',trim(local_string)
  write(0,*) integer,' :: ',trim(local_string)

end subroutine input_integer
!******************************************************************************
!******************************************************************************
subroutine input_integer_array(string,integer_array)
character(len=*),optional   :: string
character(len=str_length)   :: local_string
integer(kind=ik),dimension(:)        :: integer_array

  write(*,*)
if(present(string) ) then
   local_string=adjustl(string)//' with '//trim(string_of(INT(size(integer_array),kind=ik)))//' numbers'
else
   local_string= 'please give integer_array with'//trim(string_of(INT(size(integer_array),kind=ik)))//' numbers'
endif
   write(*,*) trim(local_string)

  read(*,*) integer_array
  print*,trim(local_string),'=',integer_array
  write(*,*)


  write(log_file_unit,*) integer_array,' :: ',trim(local_string)
  write(0,*) integer_array,' :: ',trim(local_string)

end subroutine input_integer_array
!******************************************************************************
!******************************************************************************
subroutine input_real(string,real)
character(len=*),optional   :: string
character(len=str_length)      :: local_string
real(kind=rk)        :: real

  write(*,*)
if(present(string) ) then
   local_string=adjustl(string)
else
   local_string= 'please give real: '
endif
   write(*,*) trim(local_string)

  read(*,*) real
  print*,trim(local_string),'=',real
  write(*,*)


  write(log_file_unit,*) real,' :: ',trim(local_string)
  write(0,*) real,' :: ',trim(local_string)

end subroutine input_real
!******************************************************************************
!******************************************************************************
subroutine input_char_string(string,char_string)
character(len=*),optional   :: string
character(len=str_length)      :: local_string
character(len=*)            :: char_string

  write(*,*)
if(present(string) ) then
   local_string=adjustl(string)
else
   local_string= 'please give char_string: '
endif
   write(*,*) trim(local_string)

  read(*,*) char_string
  print*,trim(local_string),'=',trim(char_string)
  write(*,*)


  write(log_file_unit,*) trim(char_string),' :: ',trim(local_string)
  write(0,*) trim(char_string),' :: ',trim(local_string)

end subroutine input_char_string
!******************************************************************************
!*****************************************************************************
subroutine print_short_integer(string,integer)
character(len=*)                      :: string
integer(kind=ik)                               :: integer

  print'(a,i5)', string,integer

end subroutine print_short_integer
!*****************************************************************************
!*****************************************************************************
subroutine print_short_real(string,real)
character(len=*)                      :: string
real(kind=rk)                  :: real

  print'(a,e15.5)', string,real

end subroutine print_short_real
!*****************************************************************************
!*****************************************************************************
subroutine print_short_complex(string,complex)
character(len=*)                      :: string
complex(kind=complex_kind)                  :: complex

  print'(a,2e15.5)', string,complex

end subroutine print_short_complex
!*****************************************************************************
!*****************************************************************************
subroutine print_short_string(string)
character(len=*)                      :: string

  print'(a)', string

end subroutine print_short_string
!*****************************************************************************
subroutine initialize_general_units()
   err_unit=-1
end subroutine initialize_general_units
!*****************************************************************************
subroutine define_error_file()

if(err_unit==-1) then
   err_unit=give_unit_of('program_ending_status.dat')
   write(err_unit,'(a)') ''
   write(err_unit,'(a)') ''
   write(err_unit,'(a)') 'here are potential errors!'
   write(err_unit,'(a)') ''
endif

end subroutine define_error_file
!*****************************************************************************
subroutine print_err_logical(string,logical,string2)
character(len=*)                      :: string
logical(kind=lk)                               :: logical
character(len=*),optional             :: string2

call define_error_file()

  write(err_unit,'(a,i5)',advance='no') string,logical
  if(present(string2))  then
     write(err_unit,'(1x,a)',advance='no') string2
  endif
  write(err_unit,'(a)') 

end subroutine print_err_logical
!*****************************************************************************
!*****************************************************************************
subroutine print_err_integer(string,integer,string2,integer2)
character(len=*)                      :: string
integer(kind=ik)                               :: integer
integer(kind=ik),optional                      :: integer2
character(len=*),optional             :: string2

call define_error_file()

  write(err_unit,'(a,i5)',advance='no') string,integer
  if(present(string2))  then
     write(err_unit,'(1x,a)',advance='no') string2
  endif
  if(present(integer2))  then
     write(err_unit,'(1x,i5)',advance='no') integer2
  endif
  write(err_unit,'(a)') 

end subroutine print_err_integer
!*****************************************************************************
subroutine print_err_integer_array(string,integer_array)
character(len=*)                      :: string
integer(kind=ik)                               :: no
integer(kind=ik),dimension(:)                  :: integer_array

call define_error_file()

  write(err_unit,'(a)',advance='no') string
  do no=1,size(integer_array)
     write(err_unit,'(1x,i5)',advance='no') integer_array(no)
  enddo
  write(err_unit,'(a)') 

end subroutine print_err_integer_array
!*****************************************************************************
!*****************************************************************************
subroutine print_err_real4(string,real,string2)
character(len=*)                      :: string
real(kind=4)                          :: real
character(len=*),optional             :: string2

call define_error_file()

  write(err_unit,'(a,e15.5)',advance='no') string,real
  if(present(string2))  then
     write(err_unit,'(1x,a)',advance='no') string2
  endif
  write(err_unit,'(a)') 

end subroutine print_err_real4
!*****************************************************************************
subroutine print_err_real8(string,real,string2)
character(len=*)                      :: string
real(kind=8)                          :: real
character(len=*),optional             :: string2

call define_error_file()

  write(err_unit,'(a,e15.5)',advance='no') string,real
  if(present(string2))  then
     write(err_unit,'(1x,a)',advance='no') string2
  endif
  write(err_unit,'(a)') 

end subroutine print_err_real8
!*****************************************************************************
!*****************************************************************************
subroutine print_err_complex(string,complex,string2)
character(len=*)                      :: string
complex(kind=complex_kind)                  :: complex
character(len=*),optional             :: string2

call define_error_file()

  write(err_unit,'(a,2e15.5)',advance='no') string,complex
  if(present(string2))  then
     write(err_unit,'(1x,a)',advance='no') string2
  endif
  write(err_unit,'(a)') 

end subroutine print_err_complex
!*****************************************************************************
!*****************************************************************************
subroutine print_err_string(string,string2,string3)
character(len=*)                      :: string
character(len=*),optional             :: string2
character(len=*),optional             :: string3

call define_error_file()

  write(err_unit,'(a)',advance='no') string
  if(present(string2))  then
     write(err_unit,'(1x,a)',advance='no') string2
  endif
  if(present(string3))  then
     write(err_unit,'(1x,a)',advance='no') string3
  endif
  write(err_unit,'(a)') 

end subroutine print_err_string
!*****************************************************************************
!*****************************************************************************
subroutine print_err_string_array(string,string_array)
character(len=*)                      :: string
character(len=*),dimension(:)         :: string_array
integer(kind=ik)                               :: no

call define_error_file()

  write(err_unit,'(a)',advance='no') string
  do no=1,size(string_array) 
     write(err_unit,'(a,a)',advance='no') ' ',trim(string_array(no))
  enddo
  write(err_unit,'(a)') 

end subroutine print_err_string_array
!*****************************************************************************
!*****************************************************************************
function strip_extension(ext_filename,extension) result(result_name)


character(len=*)                         :: ext_filename
character(len=*)                         :: extension
!PGI!character(len=len_trim(ext_filename))    :: string
!PGI!!!character(len=len_trim(ext_filename))    :: result_name
!PGI!character(len=len_trim(ext_filename)-len_trim(extension)-1)    :: result_name
character(len=str_length )    :: string
character(len=str_length )    :: result_name
integer(kind=ik)                                  :: ende
integer(kind=ik)                                  :: len_extension
integer(kind=ik)                                  :: len_string


! leading blanks stripped
string=adjustl(ext_filename)

len_extension=len_trim(extension)+1
len_string=len_trim(string)

!!print*,'filename=',trim(ext_filename)
!!print*,'string=',trim(string)
!!print*,'extension=',trim(extension)

ende=index(string,'.'//trim(extension),back=.true.)-1

!!print*,'len_extension=',len_extension
!!print*,'len_string=',len_string
!!print*,'ende=',ende

! still possible that extension ist innerpart of the ext_filename!
if(ende < 0  .or. (ende + len_extension /= len_string) ) then
print*,'strip_extension: ',trim(ext_filename),&
        &' has no extension ',trim(extension),' !'
call stop
endif

result_name=string(1:ende)

end  function strip_extension
!******************************************************************************
!******************************************************************************
subroutine stop(string)
character(len=*),optional   :: string

call define_error_file()

   write(*,*)
   write(0,*)
if( present(string) ) then
   write(*,*) 'stop: ',trim(string)
   write(0,*) 'stop: ',trim(string)
   write(err_unit,*) 'stop: ',trim(string)
else
   write(*,*) 'stop'
   write(0,*) 'stop'
   write(err_unit,*) 'stop'
endif
   write(*,*)
   write(*,*)
!! call stack_print(routine_name)
!! call stack_print(routine_name,err_unit)
stop
end subroutine stop
!******************************************************************************
!******************************************************************************
!******************************************************************************
subroutine crash(error)
integer(kind=ik)                 :: feld(10)
real(kind=rk)    :: error
integer(kind=ik)                 :: nn


print*,'crashed on artificial error to enable core file'
error=abs(error)+1.
print*,'crash error=',error
feld(1)=0_ik
do nn=1,1000
error=exp(error)+feld(int(error))
enddo

end subroutine crash
!******************************************************************************
!******************************************************************************
function give_unit_of(filename,form,action,status,iostat) result(result_unit)
! in the case 
!     that 'filename' is open the connected unit is returned
!     that 'filename' is not open the file is opened and a new unit is returned
!     that all unit numbers are connected the  program stops
!     
!     form has the same meaning as form= in open
!     action has the same meaning as action= in open

integer(kind=ik),optional                           :: iostat
integer(kind=ik)                                    :: loc_iostat

character(len=*),optional                  :: status
character(len=str_length)                  :: status_string

character(len=*),optional                  :: action
character(len=str_length)                  :: action_string

character(len=*),intent(in),optional       :: form
character(len=str_length)                  :: form_string

character(len=*)                           :: filename
character(len=str_length)                :: filename_string
!!character(len=len_trim(adjustl(filename))) :: filename_string

logical(kind=lk)                                    :: file_is_open
logical(kind=lk)                                    :: unit_is_open

integer(kind=ik)                                    :: result_unit
integer(kind=ik)                                    :: existing_file_number
integer(kind=ik)                                    :: test_unit

    

      if( .not. present(status)) then
         status_string='UNKNOWN'
      else
         status_string=adjustl(status)
      endif

      if( .not. present(form)) then
         form_string='FORMATTED'
      else
         form_string=adjustl(form)
      endif

      if( .not. present(action)) then
         action_string='READWRITE'
          else
         action_string=adjustl(action)
      endif

      filename_string=adjustl(filename)

       loc_iostat=-4711
      inquire(   file=trim(filename_string)&
             &,opened=file_is_open&
             &,number=existing_file_number&
             &)

      if( file_is_open ) then
         result_unit=existing_file_number
               print*,'file ',trim(filename_string)&
                    &,' has been already opened with unit=',result_unit
         loc_iostat=0
      else
         !! find a non occupied unit
         !! 10000 is a very large number
         result_unit=-1
         do test_unit=8,10000
            inquire(unit=test_unit,opened=unit_is_open,err=1000)
            if(.not.unit_is_open) then
               result_unit=test_unit
               goto 2000
            endif
         enddo
      
      print*
      print*,'error in give_unit_of(',trim(filename_string),'):'
      print*,'the file was not openend but no new unit has been found'
      print*,'unit=',result_unit
      print*,'form=',trim(form_string)
      print*,'action=',trim(action_string)
      print*,'status=',trim(status_string)
      print*,'iostat=',loc_iostat
      call stop('give_unit_of in interaction_module')
 2000 continue


               loc_iostat=0
         open(  unit=result_unit &
        &    ,  file=trim(filename_string) &
        &    ,  form=trim(form_string) &
        &    ,action=trim(action_string) &
        &    ,status=trim(status_string) &
        &    ,iostat=loc_iostat &
        &     )
               print'(a,a,a,i3,a,i10)','file ',trim(filename_string)&
                    &,' is initially opened with unit=',result_unit,' and iostat=',loc_iostat
      endif

if(present(iostat))  then
   iostat=loc_iostat
else
   if(loc_iostat /= 0 ) then
      print*
      print*,'error in give_unit_of(',trim(filename_string),'):'
      print*,'open was not successful'
      print*,'unit=',result_unit
      print*,'form=',trim(form_string)
      print*,'action=',trim(action_string)
      print*,'status=',trim(status_string)
      print*,'iostat=',loc_iostat
      call stop('give_unit_of in interaction_module')
   endif
endif

      return
1000  continue
      print*
      print*,'give_unit_of: error'
      print*,'no further unit can be opened'
      print*,'for file ',trim(filename_string)
      print*,'last tested unit_number=',test_unit
      call stop('give_unit_of')

end function give_unit_of
!******************************************************************************
function give_filename_of(unit) result(result_name)
integer(kind=ik)                   :: unit
logical(kind=lk)                   :: has_a_name
character(len=str_length) :: result_name
character(len=str_length) :: filename

      inquire(unit=unit,name=filename,named=has_a_name)
      if( has_a_name ) then 
         result_name=trim(filename)
      else
         result_name='unknow'//trim(string_of(unit))//'.dat'
      endif

end function give_filename_of
!******************************************************************************
function next_unused_filename(filename,extension) result(result_name)
character(len=*)      :: filename
character(len=*)      :: extension
character(len=256)    :: result_name
character(len=256)    :: temp_name
character(len=256)    :: num_string
integer(kind=ik)               :: nn
logical(kind=lk)               :: file_exists

   nn=0
do
   nn=nn+1
   write(num_string,'(i0)') nn
   result_name=trim(filename)//'_'//trim(adjustl(num_string))
   temp_name=trim(result_name)//'.'//trim(extension)
   inquire(file=trim(temp_name),exist=file_exists)
  !! print*,'file ',trim(temp_name)
   if(.not. file_exists) then
      print*,'file ',trim(temp_name), ' does not yet exist and will be taken '
      exit
   else
      print*,'file ',trim(temp_name), ' exists already'
   endif
enddo

end function next_unused_filename
!!!******************************************************************************
!!subroutine all_command_arguments(number_of_arguments,arguments) 
!!character(len=*),pointer,dimension(:)    :: arguments
!!integer(kind=ik)                                  :: number_of_arguments
!!integer(kind=ik)                                  :: arg
!!integer(kind=ik)                                  :: length
!!
!!
!!number_of_arguments=command_argument_count()
!!
!!allocate(arguments(0:number_of_arguments))
!!
!!    print*,'interaction_module.all_command_arguments: number_of_arguments=',number_of_arguments
!!
!!do arg=0,number_of_arguments-1
!!    call get_command_arguments(arg,arguments(arg),length)
!!    if(length > len(arguments(arg)) ) then
!!    print*,'interaction_module.all_command_arguments: warning:'
!!    print*,'length ',length,' of argument ',arg,' exceeds ',len(arguments(arg))
!!    print*,'the truncated argument is'
!!    print*,arguments(arg)
!!    endif
!!enddo
!!
!!end subroutine all_command_arguments
!******************************************************************************

end module cel_interaction_module
!
!****************************************************************************** 
!******************************************************************************



