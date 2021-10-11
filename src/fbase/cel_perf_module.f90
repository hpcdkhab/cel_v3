!     
! File:   cel_perf_module.f90
! Project CRESTA (see details on https://www.cresta-project.eu) Exascale library
! Send you email about the bugs to
! The High Performance Computing Center Stuttgart (HLRS) of the University of Stuttgart
! Autor: Uwe Küster
! Modified: Dmitry Khabi,  (email: khabi@hlrs.de)
! Created on Jul 10, 2013
!

module  cel_perf_module
!!use constants_module
use cel_interaction_module
use cel_omp_module
use cel_types_module
use cel_cpu_param_module
implicit none


integer(kind=ik),parameter :: perf_ik=8
integer(kind=ik),parameter :: perf_rk=rk

integer(kind=perf_ik),parameter :: i_add=1
integer(kind=perf_ik),parameter :: i_mult=2
integer(kind=perf_ik),parameter :: i_divide=3
integer(kind=perf_ik),parameter :: i_sqr=4
integer(kind=perf_ik),parameter :: i_abso=5
integer(kind=perf_ik),parameter :: i_load=6
integer(kind=perf_ik),parameter :: i_store=7
!integer,parameter :: i_oper=7






type cel_perf_type
   character(len=80),pointer,dimension(:)         :: string
   character(len=80),pointer,dimension(:)         :: properties
   integer(kind=perf_ik),pointer,dimension(:,:)   :: add
   integer(kind=perf_ik),pointer,dimension(:,:)   :: mult
   integer(kind=perf_ik),pointer,dimension(:,:)   :: divide
   integer(kind=perf_ik),pointer,dimension(:,:)   :: sqr
   integer(kind=perf_ik),pointer,dimension(:,:)   :: abso
   integer(kind=perf_ik),pointer,dimension(:,:)   :: load
   integer(kind=perf_ik),pointer,dimension(:,:)   :: store
   integer(kind=perf_ik),pointer,dimension(:,:)   :: omp_locks
   integer(kind=perf_ik),pointer,dimension(:,:)   :: number_of_repetitions
   integer(kind=perf_ik),pointer,dimension(:,:)   :: size
   integer(kind=perf_ik),pointer,dimension(:,:)   :: mat_blocks_num
   real(kind=8),pointer,dimension(:,:)            :: time
   integer(kind=perf_ik)                          :: max_type
   integer(kind=perf_ik)                          :: max_defined_type
   integer(kind=perf_ik)                          :: max_test
end type cel_perf_type

type cel_perf_counter_type
   integer(kind=perf_ik)   :: add
   integer(kind=perf_ik)   :: mult
   integer(kind=perf_ik)   :: divide
   integer(kind=perf_ik)   :: sqr
   integer(kind=perf_ik)   :: abso
   integer(kind=perf_ik)   :: load
   integer(kind=perf_ik)   :: store
   integer(kind=perf_ik)   :: omp_locks
   integer(kind=perf_ik)   :: number_of_repetitions
   integer(kind=perf_ik)   :: size
   integer(kind=perf_ik)   :: num_elements
   integer(kind=perf_ik)   :: send_bytes
   integer(kind=perf_ik)   :: get_bytes
   integer(kind=perf_ik)   :: send_messages_num
   integer(kind=perf_ik)   :: get_messages_num
   integer(kind=perf_ik)   :: num_allreduces
   real(kind=8),dimension(5)            :: time!local sum, global sum, global min, global max, global average
   real(kind=8),dimension(5)            :: time_comm!local sum, global sum, global min, global max, global average
   real(kind=8),dimension(5)            :: time_local_calc!local sum, global sum, global min, global max, global average
   real(kind=8),dimension(5)            :: time_halo_calc!local sum, global sum, global min, global max, global average
   real(kind=8),dimension(5)            :: time_jad_mv_phase_1!local sum, global sum, global min, global max, global average
   real(kind=8),dimension(5)            :: time_jad_mv_phase_2!local sum, global sum, global min, global max, global average
   real(kind=8),dimension(5)            :: time_jad_mv_phase_3!local sum, global sum, global min, global max, global average
   integer(kind=ik)                     :: cpu_freq_idx
end type cel_perf_counter_type

interface initialize
   module procedure cel_perf_initialize
end interface

interface reset
   module procedure cel_perf_counter_reset
end interface

interface print
   module procedure cel_perf_counter_print
end interface
!interface set_string_perf
!   module procedure cel_perf_set_string_perf_array
!   module procedure cel_set_string_perf_simple
!end interface

contains
!***************************************************************************************************************/
subroutine cel_perf_increment_perf_counter(perf,test,type,problem_size,add,mult,divide,sqr,&
  abs,load,store,omp_locks,repetitions,time, mat_blocks_num)
  type(cel_perf_type)  :: perf
  integer(kind=perf_ik)          :: test,type
  integer(kind=perf_ik)          :: problem_size, mat_blocks_num
  integer(kind=perf_ik) :: add, mult, divide, sqr, abs, load, store, omp_locks
  integer(kind=perf_ik)          :: repetitions
  real(kind=8)     :: time

  if(test < 1 .or. test > size(perf%size,1)) then
     print*,'cel_perf_module.cel_perf_increment_perf_counter: error'
     print*,'test=',test
     print*,'must be in intervall ',1,size(perf%size,1)
     stop
  endif
  if(type < 1 .or. type > size(perf%size,2)) then
     print*,'cel_perf_module.cel_perf_increment_perf_counter: error'
     print*,'type=',type
     print*,'must be in intervall ',1,size(perf%size,2)
     stop
  endif

  !print*,'incrementing counters for ',trim(perf%string(type))

  perf%max_defined_type=max(perf%max_defined_type,type)
  perf%size(test,type)                 =problem_size
  perf%add(test,type)                  =perf%add(test,type)   +add
  perf%mult(test,type)                 =perf%mult(test,type)  +mult
  perf%divide(test,type)               =perf%divide(test,type)+divide
  perf%sqr(test,type)                  =perf%sqr(test,type)   +sqr
  perf%abso(test,type)                 =perf%abso(test,type)  +abs
  perf%load(test,type)                 =perf%load(test,type)  +load
  perf%store(test,type)                =perf%store(test,type) +store
  perf%omp_locks(test,type)            =perf%omp_locks(test,type) +omp_locks
  perf%number_of_repetitions(test,type)=perf%number_of_repetitions(test,type) +repetitions
  perf%time(test,type)                 =perf%time(test,type)  +time
  perf%mat_blocks_num(test,type)       =mat_blocks_num

end subroutine cel_perf_increment_perf_counter
!***********************************************************************
subroutine cel_perf_initialize(perf,max_test,max_type)
  type(cel_perf_type)               :: perf
  integer(kind=perf_ik)                       :: max_type,max_test 

  perf%max_type=max_type
  perf%max_test=max_test


  !print*,'max_test,max_type=',max_test,max_type

  allocate(perf%string(max_type))
  allocate(perf%properties(max_type))
  allocate(perf%add(max_test,max_type))
  allocate(perf%mult(max_test,max_type))
  allocate(perf%divide(max_test,max_type))
  allocate(perf%sqr(max_test,max_type))
  allocate(perf%abso(max_test,max_type))
  allocate(perf%load(max_test,max_type))
  allocate(perf%store(max_test,max_type))
  allocate(perf%omp_locks(max_test,max_type))
  allocate(perf%number_of_repetitions(max_test,max_type))
  allocate(perf%time(max_test,max_type))
  allocate(perf%size(max_test,max_type))
  allocate(perf%mat_blocks_num(max_test,max_type))

  !print*,'perf%size=',size(perf%size,1),size(perf%size,2)

  perf%add   = 0
  perf%mult  = 0
  perf%divide= 0
  perf%sqr   = 0
  perf%abso  = 0
  perf%load  = 0
  perf%store = 0
  perf%omp_locks = 0
  perf%number_of_repetitions=0 
  perf%time  = 0.
  perf%max_defined_type=0
  perf%mat_blocks_num=0
end subroutine cel_perf_initialize
!***********************************************************************
subroutine cel_perf_detruct(perf)
  type(cel_perf_type)               :: perf

  perf%max_type=0
  perf%max_test=0

  !print*,'max_test,max_type=',max_test,max_type

  if(associated(perf%string))deallocate(perf%string)
  if(associated(perf%properties))deallocate(perf%properties)
  if(associated(perf%add))deallocate(perf%add)
  if(associated(perf%mult))deallocate(perf%mult)
  if(associated(perf%divide))deallocate(perf%divide)
  if(associated(perf%sqr))deallocate(perf%sqr)
  if(associated(perf%abso))deallocate(perf%abso)
  if(associated(perf%load))deallocate(perf%load)
  if(associated(perf%store))deallocate(perf%store)
  if(associated(perf%omp_locks))deallocate(perf%omp_locks)
  if(associated(perf%number_of_repetitions))deallocate(perf%number_of_repetitions)
  if(associated(perf%time))deallocate(perf%time)
  if(associated(perf%size))deallocate(perf%size)
  if(associated(perf%mat_blocks_num))deallocate(perf%mat_blocks_num)

  !print*,'perf%size=',size(perf%size,1),size(perf%size,2)

  perf%add   = 0
  perf%mult  = 0
  perf%divide= 0
  perf%sqr   = 0
  perf%abso  = 0
  perf%load  = 0
  perf%store = 0
  perf%omp_locks = 0
  perf%number_of_repetitions=0 
  perf%time  = 0.
  perf%max_defined_type=0
  perf%mat_blocks_num=0
end subroutine cel_perf_detruct

!***********************************************************************
subroutine cel_perf_counter_reset(cel_perf_counter)
  type(cel_perf_counter_type), intent(inout) :: cel_perf_counter
  
  cel_perf_counter%add = 0_perf_ik
  cel_perf_counter%mult = 0_perf_ik
  cel_perf_counter%divide = 0_perf_ik
  cel_perf_counter%sqr = 0_perf_ik
  cel_perf_counter%abso = 0_perf_ik
  cel_perf_counter%load = 0_perf_ik
  cel_perf_counter%store = 0_perf_ik
  cel_perf_counter%omp_locks = 0_perf_ik
  cel_perf_counter%number_of_repetitions = 1_perf_ik
  cel_perf_counter%size = 0_perf_ik
  cel_perf_counter%num_elements = 0_perf_ik
  cel_perf_counter%send_bytes = 0_perf_ik
  cel_perf_counter%get_bytes = 0_perf_ik
  cel_perf_counter%send_messages_num = 0_perf_ik
  cel_perf_counter%get_messages_num = 0_perf_ik
  cel_perf_counter%time = 0.0_rk
  cel_perf_counter%time_comm = 0.0_rk
  cel_perf_counter%time_local_calc = 0.0_rk
  cel_perf_counter%time_halo_calc = 0.0_rk
  cel_perf_counter%time_jad_mv_phase_1 = 0.0_rk
  cel_perf_counter%time_jad_mv_phase_2 = 0.0_rk
  cel_perf_counter%time_jad_mv_phase_3 = 0.0_rk
  cel_perf_counter%num_allreduces = 0_perf_ik
   
end subroutine cel_perf_counter_reset

subroutine cel_perf_counter_print(cel_perf_counter)
  type(cel_perf_counter_type), intent(inout) :: cel_perf_counter
  
  write(cel_output_unit,'(A)') "Performance counter data:"
  write(cel_output_unit,'(A,E13.6,A)') "add: ",real(cel_perf_counter%add,kind=rk)&
        &/(1.e6_rk), " x 10^6"
  write(cel_output_unit,'(A,E13.6,A)') "mult: ",real(cel_perf_counter%mult,kind=rk)&
        &/(1.e6_rk), " x 10^6"
  write(cel_output_unit,'(A,E13.6,A)') "divide: ",real(cel_perf_counter%divide,kind=rk)&
        &/(1.e6_rk), " x 10^6"
  write(cel_output_unit,'(A,E13.6,A)') "sqr: ",real(cel_perf_counter%sqr,kind=rk)&
        &/(1.e6_rk), " x 10^6"
  write(cel_output_unit,'(A,E13.6,A)') "abso: ",real(cel_perf_counter%abso,kind=rk)&
        &/(1.e6_rk), " x 10^6"
  write(cel_output_unit,'(A,E13.6,A)') "load: ",real(cel_perf_counter%load,kind=rk)&
        &/(1073741824.0_rk), " GiBytes"
  write(cel_output_unit,'(A,E13.6,A)') "store: ",real(cel_perf_counter%store,kind=rk)&
        &/(1073741824.0_rk), " GiBytes"
  write(cel_output_unit,'(A,E13.6,A)') "omp_locks: ",real(cel_perf_counter%omp_locks,kind=rk)&
        &, " "
  write(cel_output_unit,'(A,E13.6,A)') "number_of_repetitions: ",&
        real(cel_perf_counter%number_of_repetitions,kind=rk), ""
  write(cel_output_unit,'(A,E13.6,A)') "size: ",real(cel_perf_counter%size,kind=rk)&
        &, ""
  write(cel_output_unit,'(A,E13.6,A)') "num_alements: ",real(cel_perf_counter%num_elements,kind=rk)&
        &, ""
  write(cel_output_unit,'(A,E13.6,A)') "send_bytes: ",real(cel_perf_counter%send_bytes,kind=rk)&
        &/(1073741824.0_rk), " GiBytes"
  write(cel_output_unit,'(A,E13.6,A)') "get_bytes: ",real(cel_perf_counter%get_bytes,kind=rk)&
        &/(1073741824.0_rk), " GiBytes"
  write(cel_output_unit,'(A,E13.6,A)') "send_messages_num: ",&
        real(cel_perf_counter%send_messages_num,kind=rk)/(1.e6_rk), " x 10^6"
  write(cel_output_unit,'(A,E13.6,A)') "get_messages_num: ",&
        real(cel_perf_counter%get_messages_num,kind=rk)/(1.e6_rk), " x 10^6"
  write(cel_output_unit,'(A,E13.6,A)') "num_allreduces: ",&
        real(cel_perf_counter%num_allreduces,kind=rk), " "
  write(cel_output_unit,'(A,5(E13.6),A)') "time: ",&
        cel_perf_counter%time(1_ik),cel_perf_counter%time(2_ik),cel_perf_counter%time(3_ik),&
        cel_perf_counter%time(4_ik),cel_perf_counter%time(5_ik), " seconds"
  write(cel_output_unit,'(A,5(E13.6),A)') "time_comm: ",&
        cel_perf_counter%time_comm(1_ik),cel_perf_counter%time_comm(2_ik),cel_perf_counter%time_comm(3_ik),&
        cel_perf_counter%time_comm(4_ik),cel_perf_counter%time_comm(5_ik), " seconds"
  write(cel_output_unit,'(A,5(E13.6),A)') "time_local_calc: ",&
        cel_perf_counter%time_local_calc(1_ik), cel_perf_counter%time_local_calc(2_ik), &
        cel_perf_counter%time_local_calc(3_ik),&
        cel_perf_counter%time_local_calc(4_ik), cel_perf_counter%time_local_calc(5_ik), " seconds"
  write(cel_output_unit,'(A,5(E13.6),A)') "time_halo_calc: ",&
        cel_perf_counter%time_halo_calc(1_ik),cel_perf_counter%time_halo_calc(2_ik),&
        cel_perf_counter%time_halo_calc(3_ik),&
        cel_perf_counter%time_halo_calc(4_ik),cel_perf_counter%time_halo_calc(5_ik), " seconds"
  write(cel_output_unit,'(A,5(E13.6),A)') "time_jad_mv_phase_1: ",&
        cel_perf_counter%time_jad_mv_phase_1(1_ik),cel_perf_counter%time_jad_mv_phase_1(2_ik),&
        cel_perf_counter%time_jad_mv_phase_1(3_ik),&
        cel_perf_counter%time_jad_mv_phase_1(4_ik),cel_perf_counter%time_jad_mv_phase_1(5_ik), " seconds"
  write(cel_output_unit,'(A,5(E13.6),A)') "time_jad_mv_phase_2: ",&
        cel_perf_counter%time_jad_mv_phase_2(1_ik),cel_perf_counter%time_jad_mv_phase_2(2_ik),&
        cel_perf_counter%time_jad_mv_phase_2(3_ik),&
        cel_perf_counter%time_jad_mv_phase_2(4_ik),cel_perf_counter%time_jad_mv_phase_2(5_ik), " seconds"
  write(cel_output_unit,'(A,5(E13.6),A)') "time_jad_mv_phase_3: ",&
        cel_perf_counter%time_jad_mv_phase_3(1_ik),cel_perf_counter%time_jad_mv_phase_3(2_ik),&
        cel_perf_counter%time_jad_mv_phase_3(3_ik),&
        cel_perf_counter%time_jad_mv_phase_3(4_ik),cel_perf_counter%time_jad_mv_phase_3(5_ik), " seconds"
  
end subroutine cel_perf_counter_print

!***********************************************************************
subroutine cel_perf_set_string_perf_array(perf,string)
type(cel_perf_type)               :: perf
character(len=*),dimension(:) :: string
integer(kind=perf_ik)                       :: mm


if(size(perf%string) < size(string)) then
   print*,'error in set_string_perf:'
   print*,'size(perf%string)=',size(perf%string)
   print*,'is too small!'
   print*,'size(string)=',size(string)
   print*,'is larger'
   stop
endif

do mm=1,size(string)
  perf%string(mm)=adjustl(trim(string(mm)))
enddo


end subroutine cel_perf_set_string_perf_array
!***********************************************************************
subroutine cel_perf_set_string_perf_simple(perf,type,string)
  !subroutine cel_set_string_perf_simple(perf,type)
  type(cel_perf_type)               :: perf
  character(len=*)              :: string
  integer(kind=perf_ik)                       :: type



  !!  p!rinf*,'setting ',trim(string)
  perf%string(type)=adjustl(trim(string))



end subroutine cel_perf_set_string_perf_simple
!***********************************************************************
subroutine cel_perf_set_properties(perf,type,string)
  type(cel_perf_type)               :: perf
  character(len=*)              :: string
  integer (kind=perf_ik)                      :: type



  !!  print*,'setting ',trim(string)
  perf%properties(type)=adjustl(trim(string))



end subroutine cel_perf_set_properties

!***********************************************************************
subroutine cel_perf_incr_sum_perf(perf,type,perf_counter)
  type(cel_perf_type),intent(in)    :: perf
  integer(kind=perf_ik),intent(in)   :: type
  type(cel_perf_counter_type),intent(inout)    :: perf_counter
  integer(kind=perf_ik) ::test

  do test=1,perf%max_test
    perf_counter%add = perf_counter%add + perf%add(test,type)
    perf_counter%mult = perf_counter%mult + perf%mult(test,type)
    perf_counter%divide = perf_counter%divide + perf%divide(test,type)
    perf_counter%sqr = perf_counter%sqr + perf%sqr(test,type)
    perf_counter%abso = perf_counter%abso + perf%abso(test,type)
    perf_counter%load = perf_counter%load + perf%load(test,type)
    perf_counter%store = perf_counter%store + perf%store(test,type)
    perf_counter%omp_locks = perf_counter%omp_locks + perf%omp_locks(test,type)
    perf_counter%number_of_repetitions = perf_counter%number_of_repetitions + &
      &perf%number_of_repetitions(test,type)
    perf_counter%size = perf_counter%size + perf%size(test,type)
    perf_counter%time = perf_counter%time + perf%time(test,type)

  end do

end subroutine 
!***********************************************************************

subroutine cel_perf_write_performance_counters(cel_perf_counter,cel_omp,mem_level,&
   filename)
  type(cel_perf_counter_type),intent(in)  :: cel_perf_counter
  type(cel_omp_type),intent(in)           :: cel_omp
  integer(kind=ik),intent(in)             :: mem_level
  character(len=*),intent(in)             :: filename
  
  logical(kind=lk) :: show_header
  integer :: unit,iostat

  if(cel_omp%is_master) then
    show_header = .false._lk
    unit = 3232
    open(unit=unit,iostat=iostat,file=filename,status='NEW',action='READWRITE')
    write(cel_output_unit,'(A,A,A,I0)')'check the performance file:"',&
          trim(filename),'" exists (>0-old file exists): ', iostat
    if(iostat .eq. 0) then
      show_header = .true._lk
    else
      open(unit=unit,iostat=iostat,file=filename,status='OLD',action='READWRITE',position='append')
      write(cel_output_unit,'(A,A,A,I0)')'reopen the performance file:"',&
          trim(filename),'" status (0-ok): ', iostat
      show_header = .false._lk
    end if
    if (iostat .gt. 0) then
      print*,'performance: error, could not open the file:'
      print*,filename
    else
      if(show_header) then
        write(unit,'(A)',advance='no') "#procs;"!1
        write(unit,'(A)',advance='no') "threads;"!2
        write(unit,'(A)',advance='no') "size;"!3
        write(unit,'(A)',advance='no') "num_elements;"!4
        write(unit,'(A)',advance='no') "add;"!5
        write(unit,'(A)',advance='no') "mult;"!6
        write(unit,'(A)',advance='no') "divide;"!7
        write(unit,'(A)',advance='no') "sqr;"!8
        write(unit,'(A)',advance='no') "abso;"!9
        write(unit,'(A)',advance='no') "load;"!10
        write(unit,'(A)',advance='no') "store;"!11
        write(unit,'(A)',advance='no') "omp_locks;"!12
        write(unit,'(A)',advance='no') "number_of_repetitions;"!13
        write(unit,'(A)',advance='no') "send_bytes;"!14
        write(unit,'(A)',advance='no') "get_bytes;"!15
        write(unit,'(A)',advance='no') "send_messages_num;"!16
        write(unit,'(A)',advance='no') "get_messages_num;"!17
        write(unit,'(A)',advance='no') "loc_time_comm;"!18
        write(unit,'(A)',advance='no') "gl_time_comm;"!19
        write(unit,'(A)',advance='no') "gl_min time_comm;"!20
        write(unit,'(A)',advance='no') "gl_max_time_comm;"!21
        write(unit,'(A)',advance='no') "gl_average_time_comm;"!22
        write(unit,'(A)',advance='no') "loc_time;"!23
        write(unit,'(A)',advance='no') "gl_time;"!24
        write(unit,'(A)',advance='no') "gl_min_time;"!25
        write(unit,'(A)',advance='no') "gl_max_time;"!26
        write(unit,'(A)',advance='no') "gl_average_time;"!27
        write(unit,'(A)',advance='no') "loc_time_local_calc;"!28
        write(unit,'(A)',advance='no') "gl_time_local_calc;"!29
        write(unit,'(A)',advance='no') "gl_min_time_local_calc;"!30
        write(unit,'(A)',advance='no') "gl_max_time_local_calc;"!31
        write(unit,'(A)',advance='no') "gl_average time_local_calc;"!32
        write(unit,'(A)',advance='no') "loc_time_halo_calc;"!33
        write(unit,'(A)',advance='no') "gl_time_halo_calc;"!34
        write(unit,'(A)',advance='no') "gl_min_time_halo_calc;"!35
        write(unit,'(A)',advance='no') "gl_max_time_halo_calc;"!36
        write(unit,'(A)',advance='no') "gl_average time_halo_calc;"!37
        write(unit,'(A)',advance='no') "cpu_frequency;"!38
        write(unit,'(A)',advance='no') "loc_time_jad_mv_phase_1;"!39
        write(unit,'(A)',advance='no') "gl_time_jad_mv_phase_1;"!40
        write(unit,'(A)',advance='no') "gl_min_time_jad_mv_phase_1;"!41
        write(unit,'(A)',advance='no') "gl_max_time_jad_mv_phase_1;"!42
        write(unit,'(A)',advance='no') "gl_average time_jad_mv_phase_1;"!43
        write(unit,'(A)',advance='no') "loc_time_jad_mv_phase_2;"!44
        write(unit,'(A)',advance='no') "gl_time_jad_mv_phase_2;"!45
        write(unit,'(A)',advance='no') "gl_min_time_jad_mv_phase_2;"!46
        write(unit,'(A)',advance='no') "gl_miax_time_jad_mv_phase_2;"!47
        write(unit,'(A)',advance='no') "gl_average time_jad_mv_phase_2;"!48
        write(unit,'(A)',advance='no') "loc_time_jad_mv_phase_3;"!49
        write(unit,'(A)',advance='no') "gl_time_jad_mv_phase_3;"!50
        write(unit,'(A)',advance='no') "gl_min_time_jad_mv_phase_3;"!51
        write(unit,'(A)',advance='no') "gl_miax_time_jad_mv_phase_3;"!52
        write(unit,'(A)',advance='no') "gl_average_time_jad_mv_phase_3;"!53
        write(unit,'(A)',advance='no') "memory_level;"!54
        write(unit,'(A)',advance='no') "num_allreduces;"!55
        write(unit,'(A)',advance='yes') ""
      end if
      write(unit,'(I0,A)',advance='no') cel_omp%num_masters,";"
      write(unit,'(I0,A)',advance='no') cel_omp%num_threads,";"
      write(unit,'(I0,A)',advance='no') cel_perf_counter%size,";"
      write(unit,'(I0,A)',advance='no') cel_perf_counter%num_elements,";"
      write(unit,'(I0,A)',advance='no') cel_perf_counter%add,";"
      write(unit,'(I0,A)',advance='no') cel_perf_counter%mult,";"
      write(unit,'(I0,A)',advance='no') cel_perf_counter%divide,";"
      write(unit,'(I0,A)',advance='no') cel_perf_counter%sqr,";"
      write(unit,'(I0,A)',advance='no') cel_perf_counter%abso,";"
      write(unit,'(I0,A)',advance='no') cel_perf_counter%load,";"
      write(unit,'(I0,A)',advance='no') cel_perf_counter%store,";"
      write(unit,'(I0,A)',advance='no') cel_perf_counter%omp_locks,";"
      write(unit,'(I0,A)',advance='no') cel_perf_counter%number_of_repetitions,";"
      write(unit,'(I0,A)',advance='no') cel_perf_counter%send_bytes,";"
      write(unit,'(I0,A)',advance='no') cel_perf_counter%get_bytes,";"
      write(unit,'(I0,A)',advance='no') cel_perf_counter%send_messages_num,";"
      write(unit,'(I0,A)',advance='no') cel_perf_counter%get_messages_num,";"
      write(unit,'(E13.6,A)',advance='no') real(cel_perf_counter%time_comm(1_ik),kind=8),";"
      write(unit,'(E13.6,A)',advance='no') real(cel_perf_counter%time_comm(2_ik),kind=8),";"
      write(unit,'(E13.6,A)',advance='no') real(cel_perf_counter%time_comm(3_ik),kind=8),";"
      write(unit,'(E13.6,A)',advance='no') real(cel_perf_counter%time_comm(4_ik),kind=8),";"
      write(unit,'(E13.6,A)',advance='no') real(cel_perf_counter%time_comm(5_ik),kind=8),";"
      write(unit,'(E13.6,A)',advance='no') real(cel_perf_counter%time(1_ik),kind=8),";"
      write(unit,'(E13.6,A)',advance='no') real(cel_perf_counter%time(2_ik),kind=8),";"
      write(unit,'(E13.6,A)',advance='no') real(cel_perf_counter%time(3_ik),kind=8),";"
      write(unit,'(E13.6,A)',advance='no') real(cel_perf_counter%time(4_ik),kind=8),";"
      write(unit,'(E13.6,A)',advance='no') real(cel_perf_counter%time(5_ik),kind=8),";"
      write(unit,'(E13.6,A)',advance='no') real(cel_perf_counter%time_local_calc(1_ik),kind=8),";"
      write(unit,'(E13.6,A)',advance='no') real(cel_perf_counter%time_local_calc(2_ik),kind=8),";"
      write(unit,'(E13.6,A)',advance='no') real(cel_perf_counter%time_local_calc(3_ik),kind=8),";"
      write(unit,'(E13.6,A)',advance='no') real(cel_perf_counter%time_local_calc(4_ik),kind=8),";"
      write(unit,'(E13.6,A)',advance='no') real(cel_perf_counter%time_local_calc(5_ik),kind=8),";"
      write(unit,'(E13.6,A)',advance='no') real(cel_perf_counter%time_halo_calc(1_ik),kind=8),";"
      write(unit,'(E13.6,A)',advance='no') real(cel_perf_counter%time_halo_calc(2_ik),kind=8),";"
      write(unit,'(E13.6,A)',advance='no') real(cel_perf_counter%time_halo_calc(3_ik),kind=8),";"
      write(unit,'(E13.6,A)',advance='no') real(cel_perf_counter%time_halo_calc(4_ik),kind=8),";"
      write(unit,'(E13.6,A)',advance='no') real(cel_perf_counter%time_halo_calc(5_ik),kind=8),";"
      write(unit,'(E13.6,A)',advance='no') real(cel_perf_counter%cpu_freq_idx,kind=8),";"!real(cel_cpu_frequencies(cel_perf_counter%cpu_freq_idx+1),kind=8),";"
      write(unit,'(E13.6,A)',advance='no') real(cel_perf_counter%time_jad_mv_phase_1(1_ik),kind=8),";"
      write(unit,'(E13.6,A)',advance='no') real(cel_perf_counter%time_jad_mv_phase_1(2_ik),kind=8),";"
      write(unit,'(E13.6,A)',advance='no') real(cel_perf_counter%time_jad_mv_phase_1(3_ik),kind=8),";"
      write(unit,'(E13.6,A)',advance='no') real(cel_perf_counter%time_jad_mv_phase_1(4_ik),kind=8),";"
      write(unit,'(E13.6,A)',advance='no') real(cel_perf_counter%time_jad_mv_phase_1(5_ik),kind=8),";"
      write(unit,'(E13.6,A)',advance='no') real(cel_perf_counter%time_jad_mv_phase_2(1_ik),kind=8),";"
      write(unit,'(E13.6,A)',advance='no') real(cel_perf_counter%time_jad_mv_phase_2(2_ik),kind=8),";"
      write(unit,'(E13.6,A)',advance='no') real(cel_perf_counter%time_jad_mv_phase_2(3_ik),kind=8),";"
      write(unit,'(E13.6,A)',advance='no') real(cel_perf_counter%time_jad_mv_phase_2(4_ik),kind=8),";"
      write(unit,'(E13.6,A)',advance='no') real(cel_perf_counter%time_jad_mv_phase_2(5_ik),kind=8),";"
      write(unit,'(E13.6,A)',advance='no') real(cel_perf_counter%time_jad_mv_phase_3(1_ik),kind=8),";"
      write(unit,'(E13.6,A)',advance='no') real(cel_perf_counter%time_jad_mv_phase_3(2_ik),kind=8),";"
      write(unit,'(E13.6,A)',advance='no') real(cel_perf_counter%time_jad_mv_phase_3(3_ik),kind=8),";"
      write(unit,'(E13.6,A)',advance='no') real(cel_perf_counter%time_jad_mv_phase_3(4_ik),kind=8),";"
      write(unit,'(E13.6,A)',advance='no') real(cel_perf_counter%time_jad_mv_phase_3(5_ik),kind=8),";"
      write(unit,'(I0,A)',advance='no') mem_level,";"
      write(unit,'(I0,A)',advance='no') cel_perf_counter%num_allreduces,";"
      write(unit,'(A)',advance='yes') ""
      close(unit=unit)
    end if  
  end if  
  
end subroutine cel_perf_write_performance_counters
!***********************************************************************

subroutine cel_perf_write_performance(perf,filename)
  type(cel_perf_type) :: perf
  real(kind=8)     :: performance
  real(kind=8)     :: bandwidth
  real(kind=8)     :: operations_per_load
  real(kind=8)     :: operations_per_store, operations
  real(kind=8)     :: time_per_call
  integer(kind=8)     :: omp_locks, num_repetitions
  real(kind=8)     :: tiny_time
  integer(kind=perf_ik)                       :: max_test 
  integer(kind=perf_ik)                       :: max_defined_type
  integer(kind=perf_ik)                       :: type,test 
  integer(kind=perf_ik)                       :: problem_size, mat_blocks_num

  integer(kind=perf_ik)                       :: num
  logical(kind=lk)                       :: first_time
  character(len=256)            :: str
  character(len=*)              :: filename
  character(len=256)            :: loc_filename
  character(len=256)            :: test_filename 
  integer(kind=perf_ik)                       :: gnu_unit=1002
  integer(kind=perf_ik)                       :: test_unit=1010


  max_test=perf%max_test
  max_defined_type=perf%max_defined_type

  tiny_time=1.e-30

  loc_filename=next_unused_filename(filename,'gnu')
  !!unit=give_unit_of(trim(loc_filename)//'.dat')
  gnu_unit=give_unit_of(trim(loc_filename)//'.gnu')


  call system('mkdir '// trim(loc_filename))

  !!print*,'max_test=',max_test
  print*,'max_defined_type=',max_defined_type

                        num=0
                        first_time=.true.
  write(gnu_unit,'(a)') 'set style data points'
  write(gnu_unit,'(a)') 'set key right bottom'
  write(gnu_unit,'(a)') 'fun1(x)= log(x)/log(10)'
  write(gnu_unit,'(a)') 'fun2(x)= log(x)/log(10)'
  write(gnu_unit,'(a)') 'fun1(x)= x'
  write(gnu_unit,'(a)') 'fun2(x)= x'
  write(gnu_unit,'(a)') 'set logscale x'
  write(gnu_unit,'(a)') 'col1= 1'
  write(gnu_unit,'(a)') 'col2= 3'
  write(gnu_unit,'(a)') 'col2y2= 4'
  write(gnu_unit,'(3a)') 'set title "',trim(perf%properties(1)),'"'
  write(gnu_unit,'(3a)') 'set xlabel "matrix size"'
  write(gnu_unit,'(3a)') 'set ylabel "MFLOPs"'
  write(gnu_unit,'(3a)') 'set y2label "GB/s"'
  write(gnu_unit,'(a)') 'set yrange [0:]'
  write(gnu_unit,'(a)') 'set y2range [0:]'
  !write(gnu_unit,'(a)') 'set ytics nomirror
  write(gnu_unit,'(a)') 'set y2tics'
  write(gnu_unit,'(a)') 'set tics out'
  write(gnu_unit,'(a)') 'set autoscale  y'
  write(gnu_unit,'(a)') 'set autoscale y2'
  write(gnu_unit,'(a,a)') 'plot ',char(92)

  ! for plotting all tests together
     do type=1,max_defined_type
                    str= trim(perf%string(type))
  if(perf%number_of_repetitions(1,type) == 0 ) then
     print*,'no repetitions for ',trim(str)
     cycle
  endif
                    str= trim(perf%string(type))
                    str=str(1:index(str,' '))  ! only the part until the first blank
     test_filename= trim(perf%string(type))
     if(index(test_filename,' ')>0) test_filename=test_filename(1:index(test_filename,' '))
     test_filename=trim(loc_filename)//'/'//trim(test_filename)//'.dat'

     if(.not. first_time ) then
                        write(gnu_unit,'(a)',advance='yes') char(92)
                        write(gnu_unit,'(a)',advance='no') ','
     endif

                        num=num+1
                        first_time=.false.
                        !!write(gnu_unit,'(3a,i0)',advance='no') "'",trim(loc_filename)//'.dat',"' i ",num-1
                        write(gnu_unit,'(3a,i0)',advance='no') "'",trim(test_filename),"' i ",0
                        write(gnu_unit,'(a,a,a,a,a)',advance='no') ' u (','fun1(column(col1))','):(','fun2(column(col2))',')'
                        write(gnu_unit,'(a,a,a,a,a)',advance='no') ' t "',trim(str),' MFLOPS"'
                        write(gnu_unit,'(a,a,a,a,a)',advance='no') ' axes x1y1 '
                        
                        write(gnu_unit,'(a)',advance='yes') char(92)
                        write(gnu_unit,'(a)',advance='no') ','
                        write(gnu_unit,'(3a,i0)',advance='no') "'",trim(test_filename),"' i ",0
                        write(gnu_unit,'(a,a,a,a,a)',advance='no') ' u (','fun1(column(col1))','):(','fun2(column(col2y2))',')'
                        write(gnu_unit,'(a,a,a,a,a)',advance='no') ' t "',trim(str),' GB/s"'
                        write(gnu_unit,'(a,a,a,a,a)',advance='no') ' axes x2y2 '
                        
                        
                        
  enddo
                        write(gnu_unit,'(a)',advance='yes') ''
     write(gnu_unit,'(a)') 'pause -1'                   
  if (.false.) then
     


  ! for plotting all tests seperately
                        num=0
     do type=1,max_defined_type
     write(0,'(a,i0)') ' type=',type
  if(perf%number_of_repetitions(1,type) == 0 ) then
     cycle
  endif
  write(gnu_unit,'(3a)') 'set title "',trim(perf%properties(type)),'"'
     write(gnu_unit,'(a,a)') 'plot ',char(92)
                    str= trim(perf%string(type))
                    str=str(1:index(str,' '))  ! only the part until the first blank

     test_filename=trim(loc_filename//'.dat')
     test_filename= trim(perf%string(type))
     if(index(test_filename,' ')>0) test_filename=test_filename(1:index(test_filename,' '))
     test_filename=trim(loc_filename)//'/'//trim(test_filename)//'.dat'

                        num=num+1
                        !write(gnu_unit,'(3a,i0)',advance='no') "'",trim(test_filename),"' i ",num-1
                        write(gnu_unit,'(3a,i0)',advance='no') "'",trim(test_filename),"' i ",0
                        write(gnu_unit,'(a,a,a,a,a)',advance='no') ' u (','fun1(column(col1))','):(','fun2(column(col2))',')'
                        write(gnu_unit,'(a,a,a,a,a)',advance='no') ' t "',trim(str),' MFLOPS"'
                        write(gnu_unit,'(a,a,a,a,a)',advance='no') ' axes x1y1 '
                        
                        write(gnu_unit,'(a)',advance='yes') char(92)
                        write(gnu_unit,'(a)',advance='no') ','
                        write(gnu_unit,'(3a,i0)',advance='no') "'",trim(test_filename),"' i ",0
                        write(gnu_unit,'(a,a,a,a,a)',advance='no') ' u (','fun1(column(col1))','):(','fun2(column(col2y2))',')'
                        write(gnu_unit,'(a,a,a,a,a)',advance='no') ' t "',trim(str),' GB/s"'
                        write(gnu_unit,'(a,a,a,a,a)',advance='no') ' axes x2y2 '
                        write(gnu_unit,'(a)',advance='yes') ''
     write(gnu_unit,'(a)') 'pause -1'
  enddo

  endif

  write(0,*) ' nach gnu-file'







  do type=1,max_defined_type 


  test_filename= trim(perf%string(type))
  if(index(test_filename,' ')>0) test_filename=test_filename(1:index(test_filename,' '))  ! only the part until the first blank
  test_filename=trim(loc_filename)//'/'//trim(test_filename)//'.dat'
  test_unit=give_unit_of(test_filename)

  open(newunit=test_unit,file=test_filename)

  write(test_unit,'(a,a)') ''
  write(test_unit,'(a,a)') ''
  write(test_unit,'(a,a)') '# ',trim(perf%string(type)) 
  write(test_unit,'(a,10a    )')                   '# size   '  &
                                           & ,' time_per_call(sec) ' &
                                           & ,' perf(MFLOPS) ' &
                                           & ,' bandwidth(GB) ' &
                                           & ,' operations_per_load  ' &
                                           & ,' operations_per_store ' &
                                           & ,' operations  ' &
                                           & ,' mat_blocks_num  ' &
                                           & ,' omp_locks      ' &
                                           & ,' number of rep. ' 


  !!write(unit,'(a,a)') ''
  !!write(unit,'(a,a)') ''
  !!write(unit,'(a,a)') '# ',trim(perf%string(type)) 
  !!write(unit,'(a,5a    )')                   '# size   '  &
  !!                                         & ,' time_per_call ' &
  !!                                         & ,' performance   ' &
  !!                                         & ,' bandwidth     ' &
  !!                                         & ,' operations_per_load  ' &
  !!                                         & ,' operations_per_store '



  do test=1,max_test
     num_repetitions = perf%number_of_repetitions(test,type)
  if(num_repetitions == 0 ) then
     cycle
  endif
     operations=real(perf%add(test,type),kind=perf_rk) &
                       & +real(perf%mult(test,type),kind=perf_rk) &
                       & +real(perf%divide(test,type),kind=perf_rk) &
                       & +real(perf%sqr(test,type),kind=perf_rk)


  if(operations < 0) then
  print*
  print*,'cel_perf_write_performance: error'
  print*,'operations=',operations
  print*,'perf%add=',perf%add
  print*,'perf%mult=',perf%mult
  print*,'perf%divide=',perf%divide
  print*,'perf%sqr=',perf%sqr
  print*
  endif

     performance=operations/max(tiny_time,perf%time(test,type))*1.d-6

     bandwidth=real(perf%load(test,type)+perf%store(test,type),kind=8)/max(tiny_time,perf%time(test,type))*1.d-9
     if(perf%load(test,type) .gt. 0) then
        operations_per_load=operations/real(perf%load(test,type),kind=8)
     else
        operations_per_load=0.0
     endif
     if(perf%load(test,type) .gt. 0) then
        operations_per_store=operations/real(perf%store(test,type),kind=8)
     else
        operations_per_store=0.0
     endif
     time_per_call=perf%time(test,type)/real(num_repetitions,kind=8)
     problem_size=perf%size(test,type)
     mat_blocks_num = perf%mat_blocks_num(test,type)
     omp_locks = perf%omp_locks(test,type)

  !!write(unit,'(i10, 5e12.4)')  &
  !!                                         &  problem_size &
  !!                                         & ,time_per_call &
  !!                                         & ,performance &
  !!                                         & ,bandwidth &
  !!                                         & ,operations_per_load &
  !!                                         & ,operations_per_store
     write(test_unit,'(i12, 6e15.6, 3i12)')  &
                                           &  problem_size &
                                           & ,time_per_call &
                                           & ,performance &
                                           & ,bandwidth &
                                           & ,operations_per_load &
                                           & ,operations_per_store &
                                           & ,operations &
                                           & ,mat_blocks_num &
                                           & ,omp_locks &
                                           & ,num_repetitions
  enddo
  close(test_unit)
  enddo

  write(*,*)
  write(*,*)
  write(*,*)
  write(*,'(*(a))') 'to show the results on Unix/Linux, please run gnuplot ',trim(loc_filename)//'.gnu'
  write(*,'(*(a))') 'to show the results on Windows, please run wgnuplot ',trim(loc_filename)//'.gnu'
  write(*,*)
end subroutine cel_perf_write_performance
!***************************************************************************************************************/
subroutine gen_test_size(increase_factor,min_test_size,max_test_size,limit_for_single_step,test_size,max_test)
  real(kind=8)                             :: increase_factor
  integer(kind=4)                          :: min_test_size,max_test_size
  integer(kind=4)                          :: limit_for_single_step
  integer,allocatable,dimension(:),intent(out) :: test_size
  integer(kind=4)                          :: max_test

  integer,parameter                        :: maximal_number_of_tests=100000
  !integer,dimension(maximal_number_of_tests)     :: loc_test_size
  integer(kind=4)                          :: it
  real(kind=8)                             :: real_test_size


          write(*,*) 'gen_test_size: start'

          allocate(test_size(maximal_number_of_tests))
          real_test_size=1.01
          real_test_size=min_test_size+.01

          do it=1,maximal_number_of_tests

          test_size(it) = int(real_test_size)
          test_size(it) = max(1,int(real_test_size))
      !!    print'(a,i5,a,i8,a,e12.3)','it=',it,' test_size=',test_size(it),' real_test_size=',real_test_size

          if (test_size(it).gt.max_test_size) then
           print*,' largest possible test_size=',test_size(it)
           print*,' array dimension max_test_size=',max_test_size
           exit
          endif
          



             if( real_test_size*increase_factor <= limit_for_single_step) then
                ! The first results should be tested for small steps
                real_test_size=real_test_size+1.
             else
                real_test_size=max(real_test_size+1.,real_test_size*increase_factor)
             endif

          max_test=it

          enddo

  !!allocate(test_size(max_test))

  !!test_size(1:max_test)=loc_test_size(1:max_test)

  test_size=test_size(1:max_test)

      write(*,*) '        test_size  '
  do it=1,max_test
      write(*,'(i6,i7)') it,test_size(it)
  enddo
  write(*,*) 'gen_test_size: end'

end subroutine gen_test_size
!***********************************************************************
end module  cel_perf_module
!***********************************************************************

