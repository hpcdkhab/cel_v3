!==============================================================================
!> \file mod_vtkio.f90
!> Module for vtk output in Fortran
!>
!> The file holds a module with routines for vtk output
!>
!> \author Ralf Schneider
!> \date 11.06.2012

!==============================================================================
!> VTK output module
!> 
!> Module for VTK output of different data structures. 
!> The routines perform binary output of datasets in BIG ENDIAN see 
!>  http://vtk.1045678.n5.nabble.com/VTK-0012819-Wrong-byte-order-for-long-
!>         64bit-in-vtkDataReader-vtkDataWriter-td5092339.html
!> for details
Module cel_vtkio_module

  use cel_types_module

  Implicit None

  Private give_new_unit
  Private file_err

  interface cel_write_vtk

     Module Procedure write_vtk_structured_points_int4
     Module Procedure write_vtk_structured_points_int8
     Module Procedure write_vtk_structured_points_real8
     Module Procedure write_vtk_structured_points_real8_tensor
     Module Procedure write_vtk_unstructured_grid
     Module Procedure write_vtk_polydata_int4_1D
     Module Procedure write_vtk_polydata_int4_3D

  End interface cel_write_vtk

contains

  !============================================================================
  !> Subroutine which writes a vtk structured_points dataset
  subroutine write_vtk_unstructured_grid(nodes, elems, no_nodes, no_elems, &
       filename)

    Real(Kind=rk)    , Dimension(:,:), intent(in) :: nodes
    Integer(Kind=ik) , Dimension(:,:), intent(in) :: elems

    Integer(Kind=ik)                 , intent(in) :: no_nodes, no_elems
    Character(Len=*)                 , Intent(in) :: filename

    character(len=mcl)                :: tmp_line
    integer(kind=4)                   :: un_out, nnpe, cell_type
    integer(kind=ik)                  :: ii

    un_out = give_new_unit()
    
    if (size(elems(:,1)) == 8) then
       nnpe  = 8
       cell_type = 12
    Else if (size(elems(:,1)) == 3) then
       nnpe  = 3
       cell_type = 5
    Else
       Write(*,*)"Output of elements with ",size(elems(:,1))," nodes is not supported"
       write(*,*)"Program halted"
       Stop
    End if

    call Open_as_big_endian_stream(unit=un_out, file=trim(filename), &
         action='write', status='replace')

    write(un_out) '# vtk DataFile Version 3.0',achar(10)
    write(un_out) 'vtk output',achar(10)
    write(un_out) 'BINARY',achar(10)
    write(un_out) 'DATASET UNSTRUCTURED_GRID',achar(10)

    tmp_line=""
    write(tmp_line,"(A,1X,I0,1X,A,A)")'POINTS',no_nodes,'double',achar(10)
    write(un_out)trim(tmp_line)

    Do ii = 1, no_nodes
       Write(un_out)nodes(:,ii)
    End Do

    tmp_line=""
    write(tmp_line,"(A,1X,I0,1X,I0,A)")'CELLS',no_elems,no_elems*(nnpe+1),achar(10)
    write(un_out)trim(tmp_line)

    Do ii = 1, no_elems
       Write(un_out)nnpe,Int(elems(:,ii)-1,4)
    End Do

    tmp_line=""
    write(tmp_line,"(A,1X,I0,A)")'CELL_TYPES',no_elems,achar(10)
    write(un_out)trim(tmp_line)

    Do ii = 1, no_elems
       Write(un_out)cell_type
    End Do

    call close_big_endian_stream(un_out)

  end subroutine write_vtk_unstructured_grid

  !============================================================================
  !> Subroutine which writes a vtk structured_points dataset
  subroutine write_vtk_unstructured_grid_nodelist(nodes, no_nodes, &
       filename)

    Real(Kind=rk)    , Dimension(:,:), intent(in) :: nodes

    Integer(Kind=ik)                 , intent(in) :: no_nodes
    Character(Len=*)                 , Intent(in) :: filename

    character(len=mcl)                :: tmp_line
    integer(kind=4)                   :: un_out
    integer(kind=ik)                  :: ii

    un_out = give_new_unit()
    
    call Open_as_big_endian_stream(unit=un_out, file=trim(filename), &
         action='write', status='replace')

    write(un_out) '# vtk DataFile Version 3.0',achar(10)
    write(un_out) 'vtk output',achar(10)
    write(un_out) 'BINARY',achar(10)
    write(un_out) 'DATASET UNSTRUCTURED_GRID',achar(10)

    tmp_line=""
    write(tmp_line,"(A,1X,I0,1X,A,A)")'POINTS',no_nodes,'double',achar(10)
    write(un_out)trim(tmp_line)

    Do ii = 1, no_nodes
       Write(un_out)nodes(:,ii)
    End Do

    call close_big_endian_stream(un_out)

  end subroutine write_vtk_unstructured_grid_nodelist

  !============================================================================
  !> Subroutine which writes a vtk structured_points dataset
  subroutine write_vtk_structured_points_int4(matrix, filename, extend, &
                                              spacing, origin, desc, &
                                              head, init)

    integer(kind=4) , Dimension(:,:,:), intent(in) :: matrix
    character(len=*)                  , intent(in) :: filename
    Real(kind=rk)   , Dimension(3)    , intent(in) :: spacing,origin
    integer(kind=ik), Dimension(3)    , intent(in) :: extend
    character(len=*), optional        , intent(in) :: desc
    Logical         , optional        , intent(in) :: head, init

    character(len=mcl)               :: tmp_line, tmp_origin,tmp_real
    character(len=mcl)               :: tmp_pointdata
    character(len=12)                :: loc_desc
    integer(kind=4)                  :: un_out
    logical                          :: loc_head, loc_init

    integer(kind=ik)                 :: pointdata

    if (present(init)) then
       loc_init = init
    Else
       loc_init = .FALSE.
    End if

    if (present(head)) then
       loc_head = head
    Else
       loc_head = .FALSE.
    End if

    if (present(desc)) then
       loc_desc = desc(1:min(12,len_trim(desc)))
    Else
       loc_desc = "scalars_i4"
    End if

    un_out = give_new_unit()
    if (loc_init) then
       call open_as_big_endian_stream(unit=un_out, file=trim(filename), &
            action='write', status='replace')
       call write_vtk_head(un_out)
    Else
       call open_as_big_endian_stream(unit=un_out, file=trim(filename), &
            action='write', status='old', position='append')
    End if


    If (loc_head) then    
       !# -------
       !# Header 
       call write_vtk_structured_points_head(un_out,  extend, &
                                             spacing, origin   )    
    End If
    !# -------
    !# Scalars
    write(un_out)char(10)
    write(un_out) 'SCALARS '//trim(loc_desc)//' int',achar(10)

    !# -------
    !# Lockup_table
    write(un_out) 'LOOKUP_TABLE default',achar(10)


    write(un_out) matrix
    call close_big_endian_stream(un_out)

  end subroutine write_vtk_structured_points_int4

  !============================================================================
  !> Subroutine which writes a vtk structured_points dataset
  subroutine write_vtk_structured_points_int8(matrix, filename, extend, &
                                              spacing, origin, desc, &
                                              head, init)

    integer(kind=8) , Dimension(:,:,:), intent(in) :: matrix
    character(len=*)                  , intent(in) :: filename
    Real(kind=rk)   , Dimension(3)    , intent(in) :: spacing,origin
    integer(kind=ik), Dimension(3)    , intent(in) :: extend
    character(len=*), optional        , intent(in) :: desc
    Logical         , optional        , intent(in) :: head, init

    character(len=mcl)               :: tmp_line, tmp_origin,tmp_real
    character(len=mcl)               :: tmp_pointdata
    character(len=12)                :: loc_desc
    integer(kind=4)                  :: un_out
    logical                          :: loc_head, loc_init

    integer(kind=ik)                 :: pointdata

    if (present(init)) then
       loc_init = init
    Else
       loc_init = .FALSE.
    End if

    if (present(head)) then
       loc_head = head
    Else
       loc_head = .FALSE.
    End if
    
    if (present(desc)) then
       loc_desc = desc(1:min(12,len_trim(desc)))
    Else
       loc_desc = "scalars_i8"
    End if

    un_out = give_new_unit()
    if (loc_init) then
       call open_as_big_endian_stream(unit=un_out, file=trim(filename), &
            action='write', status='replace')
       call write_vtk_head(un_out)
    Else
       call open_as_big_endian_stream(unit=un_out, file=trim(filename), &
            action='write', status='old', position='append')
    End if

    If (loc_head) then
       !# -------
       !# Header 
       call write_vtk_structured_points_head(un_out,  extend, &
                                             spacing, origin   )    
    End If

    !# -------
    !# Scalars
    write(un_out)char(10)
    write(un_out) 'SCALARS '//trim(loc_desc)//' long',achar(10)
    !# -------
    !# Lockup_table
    write(un_out) 'LOOKUP_TABLE default',achar(10)


    write(un_out) matrix
    call close_big_endian_stream(un_out)

  end subroutine write_vtk_structured_points_int8

  !============================================================================
  !> Subroutine which writes a vtk structured_points dataset
  subroutine write_vtk_structured_points_real8(matrix,  filename, extend, &
                                               spacing, origin, desc, &
                                               head, init)

    Real(kind=8)    , Dimension(:,:,:), intent(in) :: matrix
    character(len=*)                  , intent(in) :: filename
    Real(kind=rk)   , Dimension(3)    , intent(in) :: spacing,origin
    integer(kind=ik), Dimension(3)    , intent(in) :: extend
    character(len=*), optional        , intent(in) :: desc
    Logical         , optional        , intent(in) :: head, init

    character(len=mcl)               :: tmp_line, tmp_origin,tmp_real
    character(len=mcl)               :: tmp_pointdata
    character(len=12)                :: loc_desc
    integer(kind=4)                  :: un_out
    logical                          :: loc_head, loc_init

    integer(kind=ik)                 :: pointdata

    if (present(init)) then
       loc_init = init
    Else
       loc_init = .FALSE.
    End if

    if (present(head)) then
       loc_head = head
    Else
       loc_head = .FALSE.
    End if

    if (present(desc)) then
       loc_desc = desc(1:min(12,len_trim(desc)))
    Else
       loc_desc = "scalars_r8"
    End if

    un_out = give_new_unit()
    if (loc_init) then
       call open_as_big_endian_stream(unit=un_out, file=trim(filename), &
            action='write', status='replace')
       call write_vtk_head(un_out)
    Else
       call open_as_big_endian_stream(unit=un_out, file=trim(filename), &
            action='write', status='old', position='append')
    End if

    If (loc_head) then
       !# -------
       !# Header 
       call write_vtk_structured_points_head(un_out,  extend, &
                                             spacing, origin   )
    End If

    !# -------
    !# Scalars
    write(un_out)char(10)
    write(un_out) 'SCALARS '//trim(loc_desc)//' double',achar(10)
    !# -------
    !# Lockup_table
    write(un_out) 'LOOKUP_TABLE default',achar(10)


    write(un_out) matrix
    call close_big_endian_stream(un_out)

  end subroutine write_vtk_structured_points_real8

  !============================================================================
  !> Subroutine which writes a vtk structured_points dataset
  subroutine write_vtk_structured_points_real8_3D1D(matrix,  filename, extend, &
                                               spacing, origin, desc, &
                                               head, init)

    Real(kind=8)    , Dimension(:)          , intent(in) :: matrix
    character(len=*)                        , intent(in) :: filename
    Real(kind=rk)   , Dimension(3), optional, intent(in) :: spacing,origin
    integer(kind=ik), Dimension(3), optional, intent(in) :: extend
    character(len=*), optional              , intent(in) :: desc
    Logical         , optional              , intent(in) :: head, init

    character(len=mcl)               :: tmp_line, tmp_origin,tmp_real
    character(len=mcl)               :: tmp_pointdata
    character(len=12)                :: loc_desc
    integer(kind=4)                  :: un_out
    logical                          :: loc_head, loc_init
    Real(kind=rk)   , Dimension(3)   :: loc_spacing, loc_origin
    integer(kind=ik), Dimension(3)   :: loc_extend

    integer(kind=ik)                 :: pointdata

    !==========================================================================

    !** Check for optional parameters *****************************************
    if (present(init)) then
       loc_init = init
    Else
       loc_init = .FALSE.
    End if

    if (present(head)) then
       loc_head = head
    Else
       loc_head = .FALSE.
    End if

    if (present(desc)) then
       loc_desc = desc(1:min(12,len_trim(desc)))
    Else
       loc_desc = "scalars_r8"
    End if

    if (present(extend)) then
       loc_extend = extend
    Else
       loc_extend(1) = size(matrix)
       loc_extend(2) = 1_ik
       loc_extend(3) = 1_ik
    End if

    if (present(spacing)) then
       loc_spacing = spacing
    Else
       loc_spacing = [1._rk,1._rk,1._rk]
    End if
    
    if (present(origin)) then
       loc_origin = origin
    Else
       loc_origin = [1._rk,1._rk,1._rk]
    End if

    !** Check on integrity of optional parameters *****************************
    If( (.not.loc_init .AND. loc_head) .OR. (loc_init .AND. .not.loc_head) ) then

       Write(*,*)"There is one of the parameters head or init present !!"
       Write(*,*)"with the other being absent                         !!"
       write(*,*)"This case is not meaningfull                        !!"
       Write(*,*)"Please fix your call to this routine.               !!"
       write(*,*)
       write(*,*)"PROGRAM HALTED"
       stop

    End If

    If( loc_head .AND. ( .NOT.present(origin)  .or. &
                         .NOT.present(spacing) .or. &
                         .NOT.present(extend)        ) ) then

       Write(*,*)"There is one of the parameters origin, spacing, extend !!"
       Write(*,*)"missing with the parameter head being present.         !!"
       write(*,*)"This case is not meaningfull                           !!"
       Write(*,*)"Please fix your call to this routine.                  !!"
       write(*,*)
       write(*,*)"PROGRAM HALTED"
       stop
       
    End If
       
    un_out = give_new_unit()
    if (loc_init) then
       call open_as_big_endian_stream(unit=un_out, file=trim(filename), &
            action='write', status='replace')
       call write_vtk_head(un_out)
    Else
       call open_as_big_endian_stream(unit=un_out, file=trim(filename), &
            action='write', status='old', position='append')
    End if

    If (loc_head) then
       !# -------
       !# Header 
       call write_vtk_structured_points_head(un_out,  extend, &
                                             spacing, origin   )
    End If

    !# -------
    !# Scalars
    write(un_out)char(10)
    write(un_out) 'SCALARS '//trim(loc_desc)//' double',achar(10)
    !# -------
    !# Lockup_table
    write(un_out) 'LOOKUP_TABLE default',achar(10)

    write(un_out) matrix
    call close_big_endian_stream(un_out)

  end subroutine write_vtk_structured_points_real8_3D1D

  !============================================================================
  !> Subroutine which writes a vtk structured_points dataset
  subroutine write_vtk_structured_points_real8_vector (&
       matrix,  filename, extend, spacing, origin, desc, head, init)

    Real(kind=8)    , Dimension(:,:,:,:)    , intent(in) :: matrix
    character(len=*)                        , intent(in) :: filename
    Real(kind=rk)   , Dimension(3), optional, intent(in) :: spacing,origin
    integer(kind=ik), Dimension(4), optional, intent(in) :: extend
    character(len=*), optional              , intent(in) :: desc
    Logical         , optional              , intent(in) :: head, init

    character(len=mcl)               :: tmp_line, tmp_origin,tmp_real
    character(len=mcl)               :: tmp_pointdata
    character(len=12)                :: loc_desc
    integer(kind=4)                  :: un_out
    logical                          :: loc_head, loc_init
    Real(kind=rk)   , Dimension(3)   :: loc_spacing, loc_origin
    integer(kind=ik), Dimension(4)   :: loc_extend
    integer(kind=ik)                 :: pointdata

    !** Check for optional parameters *****************************************
    if (present(init)) then
       loc_init = init
    Else
       loc_init = .FALSE.
    End if

    if (present(head)) then
       loc_head = head
    Else
       loc_head = .FALSE.
    End if

    if (present(desc)) then
       loc_desc = desc(1:min(12,len_trim(desc)))
    Else
       loc_desc = "vectors_r8"
    End if

    if (present(extend)) then
       loc_extend = extend
    Else
       loc_extend(1) = size(matrix)
       loc_extend(2) = 1_ik
       loc_extend(3) = 1_ik
       loc_extend(4) = 1_ik
    End if

    if (present(spacing)) then
       loc_spacing = spacing
    Else
       loc_spacing = [1._rk,1._rk,1._rk]
    End if
    
    if (present(origin)) then
       loc_origin = origin
    Else
       loc_origin = [1._rk,1._rk,1._rk]
    End if

    !** Check on integrity of optional parameters *****************************
    If( (.not.loc_init .AND. loc_head) .OR. (loc_init .AND. .not.loc_head) ) then

       Write(*,*)"There is one of the parameters head or init present !!"
       Write(*,*)"with the other being absent                         !!"
       write(*,*)"This case is not meaningfull                        !!"
       Write(*,*)"Please fix your call to this routine.               !!"
       write(*,*)
       write(*,*)"PROGRAM HALTED"
       stop

    End If

    If( loc_head .AND. ( .NOT.present(origin)  .or. &
                         .NOT.present(spacing) .or. &
                         .NOT.present(extend)        ) ) then

       Write(*,*)"There is one of the parameters origin, spacing, extend !!"
       Write(*,*)"missing with the parameter head being present.         !!"
       write(*,*)"This case is not meaningfull                           !!"
       Write(*,*)"Please fix your call to this routine.                  !!"
       write(*,*)
       write(*,*)"PROGRAM HALTED"
       stop
       
    End If

    un_out = give_new_unit()
    if (loc_init) then
       call open_as_big_endian_stream(unit=un_out, file=trim(filename), &
            action='write', status='replace')
       call write_vtk_head(un_out)
    Else
       call open_as_big_endian_stream(unit=un_out, file=trim(filename), &
            action='write', status='old', position='append')
    End if

    If (loc_head) then
       !# ------- 
       !# Header 
       call write_vtk_structured_points_head(un_out,  extend(2:4), &
                                             spacing, origin    )
    End If

    !# -------
    !# Scalars
    write(un_out)char(10)
    write(un_out) 'VECTORS '//trim(loc_desc)//' double',achar(10)

    write(un_out) matrix
    call close_big_endian_stream(un_out)

  end subroutine write_vtk_structured_points_real8_vector

  !============================================================================
  !> Subroutine which writes a vtk structured_points dataset
  subroutine write_vtk_structured_points_real8_tensor (&
       matrix,  filename, extend, spacing, origin, desc, head, init)

    Real(kind=8)    , Dimension(:,:,:,:,:), intent(in) :: matrix
    character(len=*)                      , intent(in) :: filename
    Real(kind=rk)   , Dimension(3)        , intent(in) :: spacing,origin
    integer(kind=ik), Dimension(5)        , intent(in) :: extend
    character(len=*), optional            , intent(in) :: desc
    Logical         , optional            , intent(in) :: head, init

    character(len=mcl)               :: tmp_line, tmp_origin,tmp_real
    character(len=mcl)               :: tmp_pointdata
    character(len=12)                :: loc_desc
    integer(kind=4)                  :: un_out
    logical                          :: loc_head, loc_init

    integer(kind=ik)                 :: pointdata

    if (present(init)) then
       loc_init = init
    Else
       loc_init = .FALSE.
    End if

    if (present(head)) then
       loc_head = head
    Else
       loc_head = .FALSE.
    End if

    if (present(desc)) then
       loc_desc = desc(1:min(12,len_trim(desc)))
    Else
       loc_desc = "tensors_r8"
    End if

    un_out = give_new_unit()
    if (loc_init) then
       call open_as_big_endian_stream(unit=un_out, file=trim(filename), &
            action='write', status='replace')
       call write_vtk_head(un_out)
    Else
       call open_as_big_endian_stream(unit=un_out, file=trim(filename), &
            action='write', status='old', position='append')
    End if

    If (loc_head) then
       !# -------
       !# Header 
       call write_vtk_structured_points_head(un_out,  extend(3:5), &
                                             spacing, origin    )
    End If

    !# -------
    !# Scalars
    write(un_out)char(10)
    write(un_out) 'TENSORS '//trim(loc_desc)//' double',achar(10)

    write(un_out) matrix
    call close_big_endian_stream(un_out)

  end subroutine write_vtk_structured_points_real8_tensor

  !============================================================================
  !> Subroutine which writes a head for a vtk structured_points dataset
  subroutine write_vtk_structured_points_head(un_out,  extend, &
                                              spacing, origin   ) 

    integer(kind=4)                   , intent(in) :: un_out
    Real(kind=rk)   , Dimension(3)    , intent(in) :: spacing,origin
    integer(kind=ik), Dimension(3)    , intent(in) :: extend

    character(len=mcl)               :: tmp_line, tmp_origin,tmp_real
    character(len=mcl)               :: tmp_pointdata

    integer(kind=ik)                 :: pointdata

    write(un_out) 'DATASET STRUCTURED_POINTS',achar(10)

    !# ---------
    !# Dimension
    tmp_real=""
    write(tmp_real,'(A,I0,1X,I0,1X,I0,A)')'DIMENSIONS ',&
         extend(1),extend(2),extend(3),achar(10)
    write(un_out) trim(tmp_real)

    !# ---------
    !# spacing
    tmp_line=""
    Write(tmp_line,'(A,F0.2,1X,F0.2,1X,F0.2,A)')'SPACING ', &
         spacing(1),spacing(2),spacing(3),achar(10)
    write(un_out)trim(tmp_line)

    !# ---------
    !# origin
    tmp_origin=""
    write(tmp_origin,'(A,F0.2,1X,F0.2,1X,F0.2,A)') 'ORIGIN ',&
         origin(1),origin(2),origin(3),achar(10)
    write(un_out)trim(tmp_origin)
 
    !# -------
    !# pointdata
    pointdata = extend(1) * extend(2) * extend(3)
    tmp_pointdata=""
    write(tmp_pointdata,'(A,1X,I0,A)') 'POINT_DATA',pointdata
    write(un_out)trim(tmp_pointdata)

  End subroutine write_vtk_structured_points_head

  !============================================================================
  !> Subroutine which writes a vtk structured_points dataset
  subroutine write_vtk_polydata_int4_1D (matrix, filename, desc, head)

    Integer(kind=4) , Dimension(:)            , intent(in) :: matrix
    character(len=*)                          , intent(in) :: filename
    character(len=*), optional                , intent(in) :: desc
    Logical         , optional                , intent(in) :: head

    character(len=mcl)               :: tmp_pointdata
    character(len=12)                :: loc_desc
    integer(kind=4)                  :: un_out
    logical                          :: loc_head, loc_init

    integer(kind=ik)                 :: pointdata

    if (present(head)) then
       loc_head = head
    Else
       loc_head = .FALSE.
    End if

    if (present(desc)) then
       loc_desc = desc(1:min(12,len_trim(desc)))
    Else
       loc_desc = "scalars_i8"
    End if

    un_out = give_new_unit()
    call open_as_big_endian_stream(unit=un_out, file=trim(filename), &
                                   action='write', status='old', &
                                   position='append')
    If (loc_head) then
       !# -------
       !# Header 
       write(tmp_pointdata,'(A,1X,I0,A)') 'POINT_DATA', size(matrix)
       write(un_out)trim(tmp_pointdata)

    End If

    !# Scalars
    write(un_out)char(10)
    write(un_out) 'SCALARS '//trim(loc_desc)//' int',achar(10)
    !# -------
    !# Lockup_table
    write(un_out) 'LOOKUP_TABLE default',achar(10)

    write(un_out) matrix
    call close_big_endian_stream(un_out)

  end subroutine write_vtk_polydata_int4_1D

  !============================================================================
  !> Subroutine which writes a vtk structured_points dataset
  subroutine write_vtk_polydata_int8_1D (matrix, filename, desc, head)

    Integer(kind=8) , Dimension(:)            , intent(in) :: matrix
    character(len=*)                          , intent(in) :: filename
    character(len=*), optional                , intent(in) :: desc
    Logical         , optional                , intent(in) :: head

    character(len=mcl)               :: tmp_pointdata
    character(len=12)                :: loc_desc
    integer(kind=4)                  :: un_out
    logical                          :: loc_head, loc_init

    integer(kind=ik)                 :: pointdata

    if (present(head)) then
       loc_head = head
    Else
       loc_head = .FALSE.
    End if

    if (present(desc)) then
       loc_desc = desc(1:min(12,len_trim(desc)))
    Else
       loc_desc = "scalars_i8"
    End if
 
    un_out = give_new_unit()
    call open_as_big_endian_stream(unit=un_out, file=trim(filename), &
                                   action='write', status='old', &
                                   position='append')
    If (loc_head) then
       !# -------
       !# Header 
       write(tmp_pointdata,'(A,1X,I0,A)') 'POINT_DATA', size(matrix)
       write(un_out)trim(tmp_pointdata)

    End If

    !# Scalars
    write(un_out)char(10)
    write(un_out) 'SCALARS '//trim(loc_desc)//' int',achar(10)
    !# -------
    !# Lockup_table
    write(un_out) 'LOOKUP_TABLE default',achar(10)

    write(un_out) Int(matrix,4)
    call close_big_endian_stream(un_out)

  end subroutine write_vtk_polydata_int8_1D

  !============================================================================
  !> Subroutine which writes a vtk structured_points dataset
  subroutine write_vtk_polydata_int4_3D (matrix, filename, desc, head)

    Integer(kind=4) , Dimension(:,:,:)        , intent(in) :: matrix
    character(len=*)                          , intent(in) :: filename
    character(len=*), optional                , intent(in) :: desc
    Logical         , optional                , intent(in) :: head

    character(len=mcl)               :: tmp_pointdata
    character(len=12)                :: loc_desc
    integer(kind=4)                  :: un_out
    logical                          :: loc_head, loc_init

    integer(kind=ik)                 :: pointdata

    if (present(head)) then
       loc_head = head
    Else
       loc_head = .FALSE.
    End if

    if (present(desc)) then
       loc_desc = desc(1:min(12,len_trim(desc)))
    Else
       loc_desc = "scalars_i8"
    End if

    un_out = give_new_unit()
    call open_as_big_endian_stream(unit=un_out, file=trim(filename), &
                                   action='write', status='old', &
                                   position='append')
    If (loc_head) then
       !# -------
       !# Header 
       write(tmp_pointdata,'(A,1X,I0,A)') 'POINT_DATA', size(matrix)
       write(un_out)trim(tmp_pointdata)

    End If

    !# Scalars
    write(un_out)char(10)
    write(un_out) 'SCALARS '//trim(loc_desc)//' int',achar(10)
    !# -------
    !# Lockup_table
    write(un_out) 'LOOKUP_TABLE default',achar(10)

    write(un_out) matrix
    call close_big_endian_stream(un_out)

  end subroutine write_vtk_polydata_int4_3D

  !============================================================================
  !> Subroutine which writes a vtk Polydata Real-8 vector dataset
  subroutine write_vtk_polydata_Real8_1D (matrix, filename, desc, head)

    Real(kind=8)    , Dimension(:)            , intent(in) :: matrix
    character(len=*)                          , intent(in) :: filename
    character(len=*), optional                , intent(in) :: desc
    Logical         , optional                , intent(in) :: head

    character(len=mcl)               :: tmp_pointdata
    character(len=12)                :: loc_desc
    integer(kind=4)                  :: un_out
    logical                          :: loc_head, loc_init

    integer(kind=ik)                 :: pointdata

    if (present(head)) then
       loc_head = head
    Else
       loc_head = .FALSE.
    End if

    if (present(desc)) then
       loc_desc = desc(1:min(12,len_trim(desc)))
    Else
       loc_desc = "scalars_r8"
    End if

    un_out = give_new_unit()
    call open_as_big_endian_stream(unit=un_out, file=trim(filename), &
                                   action='write', status='old', &
                                   position='append')
    If (loc_head) then
       !# -------
       !# Header 
       write(tmp_pointdata,'(A,1X,I0,A)') 'POINT_DATA', size(matrix(:))
       write(un_out)trim(tmp_pointdata)

    End If

    !# Vectors
    write(un_out)char(10)
    write(un_out) 'SCALARS '//trim(loc_desc)//' double',achar(10)
    !# -------
    !# Lockup_table
    write(un_out) 'LOOKUP_TABLE default',achar(10)
    !# -------
    write(un_out) matrix
    call close_big_endian_stream(un_out)

  end subroutine write_vtk_polydata_Real8_1D

  !============================================================================
  !> Subroutine which writes a vtk Polydata Real-8 vector dataset
  subroutine write_vtk_polydata_Real8_vec (matrix, filename, desc, head)

    Real(kind=8)    , Dimension(:,:)          , intent(in) :: matrix
    character(len=*)                          , intent(in) :: filename
    character(len=*), optional                , intent(in) :: desc
    Logical         , optional                , intent(in) :: head

    character(len=mcl)               :: tmp_pointdata
    character(len=12)                :: loc_desc
    integer(kind=4)                  :: un_out
    logical                          :: loc_head, loc_init

    integer(kind=ik)                 :: pointdata

    if (present(head)) then
       loc_head = head
    Else
       loc_head = .FALSE.
    End if

    if (present(desc)) then
       loc_desc = desc(1:min(12,len_trim(desc)))
    Else
       loc_desc = "scalars_i8"
    End if

    un_out = give_new_unit()
    call open_as_big_endian_stream(unit=un_out, file=trim(filename), &
                                   action='write', status='old', &
                                   position='append')
    If (loc_head) then
       !# -------
       !# Header 
       write(tmp_pointdata,'(A,1X,I0,A)') 'POINT_DATA', size(matrix(1,:))
       write(un_out)trim(tmp_pointdata)

    End If

    !# Vectors
    write(un_out)char(10)
    write(un_out) 'VECTORS '//trim(loc_desc)//' double',achar(10)
    !# -------
    write(un_out) matrix
    call close_big_endian_stream(un_out)

  end subroutine write_vtk_polydata_Real8_vec

  !============================================================================
  !> Subroutine which writes a vtk structured_points dataset
  subroutine write_vtk_polydata_grid(grids,  filename, orig, init)

    Real(kind=rk)   , Dimension(:,:)          , intent(in) :: grids
    character(len=*)                          , intent(in) :: filename
    Real(kind=ik)   , Dimension(3), optional  , intent(in) :: orig
    Logical         , optional                , intent(in) :: init

    character(len=mcl)               :: tmp_char
    character(len=12)                :: loc_desc
    integer(kind=4)                  :: un_out
    logical                          :: loc_init

    integer(kind=ik)                 :: no_points,ii

    if (present(init)) then
       loc_init = init
    Else
       loc_init = .FALSE.
    End if

    no_points = size(grids(1,:))

    un_out = give_new_unit()
    if (loc_init) then

       call open_as_big_endian_stream(unit=un_out, file=trim(filename), &
            action='write', status='replace')
       call write_vtk_head(un_out)

    Else

       call open_as_big_endian_stream(unit=un_out, file=trim(filename), &
            action='write', status='old', position='append')

    End if

    write(un_out) 'DATASET POLYDATA',achar(10)

    write(tmp_char,'(A,I0,A)')"POINTS ",no_points," double"
    write(un_out)trim(tmp_char),achar(10)

    if (present(orig)) then
       Do ii = 1, no_points
          write(un_out)grids(:,ii)+orig
       End Do
    Else
       write(un_out)grids
    End if

    write(un_out)achar(10)

    call close_big_endian_stream(un_out)

  end subroutine write_vtk_polydata_grid

  !============================================================================
  !> Subroutine which writes a vtk structured_points dataset
  subroutine write_vtk_coo_matrix(cols, rows, values,  filename, init)

    integer(kind=ik)   , Dimension(:)          , intent(in) :: cols
    integer(kind=ik)   , Dimension(:)          , intent(in) :: rows
    Real(kind=rk)      , Dimension(:)          , intent(in) :: values
    character(len=*)                          , intent(in) :: filename
    Logical         , optional                , intent(in) :: init

    character(len=mcl)               :: tmp_char
    character(len=12)                :: loc_desc
    integer(kind=4)                  :: un_out
    logical                          :: loc_init

    integer(kind=ik)                 :: no_points,ii

    if (present(init)) then
       loc_init = init
    Else
       loc_init = .FALSE.
    End if

    no_points = size(cols)

    un_out = give_new_unit()
    if (loc_init) then

       call open_as_big_endian_stream(unit=un_out, file=trim(filename), &
            action='write', status='replace')
       call write_vtk_head(un_out)

    Else

       call open_as_big_endian_stream(unit=un_out, file=trim(filename), &
            action='write', status='old', position='append')

    End if

    write(un_out) 'DATASET POLYDATA',achar(10)

    write(tmp_char,'(A,I0,A,A)')"POINTS ",no_points," double",achar(10)
    write(un_out)trim(tmp_char),achar(10)

    Do ii = 1, no_points
        write(un_out) REAL(cols(ii),kind=rk),&
         REAL(rows(ii),kind=rk), 0.0_rk
    End Do

    write(un_out)achar(10)

    call close_big_endian_stream(un_out)
    call write_vtk_polydata_Real8_1D (values, filename, "values", .true.)

  end subroutine write_vtk_coo_matrix
  !============================================================================
  !> Subroutine which writes a head for a vtk binary file
  subroutine write_vtk_head(un_out)

    integer(kind=4)                           , intent(in) :: un_out

    write(un_out) '# vtk DataFile Version 3.0',achar(10)
    write(un_out) 'vtk output',achar(10)
    write(un_out) 'BINARY',achar(10)

  End subroutine write_vtk_head

  !############################################################################
  !############################################################################
  !> \name vtkio private auxilliary routines
  !> @{
  !> These routines exist in most codes. So they are declared private to aviod
  !> interference with other implementations

  !============================================================================
  !> Function which returns new free unit
  function give_new_unit() result(new_unit)

    Integer(kind=4) :: new_unit

    Integer(kind=4) :: ii

    Logical :: unit_is_open

    Do ii = 3000, huge(new_unit)-1

       inquire(unit=ii, opened=unit_is_open)

       if( .not.unit_is_open ) then
          new_unit = ii
          Exit
       end if

    End Do

    if ( unit_is_open ) then

       WRITE(*,*)
       WRITE(*,*)'Something bad and unexpected happened during search ',&
            'for free unit'
       WRITE(*,*)'Could not find a new unit between 3000 and huge(Int(kind=4))'
       WRITE(*,*)' '
       WRITE(*,*)'PROGRAM STOPPED'
       STOP
    END IF

  End function give_new_unit

  !============================================================================
  !> Subroutine for I/O error handling while operating on files
  SUBROUTINE file_err(in_file,io_stat)

    INTEGER             :: io_stat
    CHARACTER (LEN=*)   :: in_file

    IF (io_stat /= 0) Then
       WRITE(*,*)
       WRITE(*,"(80('='))")
       WRITE(*,"('EE ',A,T77,' EE')")   'Operation on file :'       
       WRITE(*,"('EE ',A          )")   in_file
       WRITE(*,"('EE ',A,T77,' EE')")   'faild !!'
       WRITE(*,"('EE ',A,I0,T77,' EE')")'With I/O Status ',io_stat
       WRITE(*,"('EE PROGRAM STOPPED ..... ',T77,' EE',/,'<',77('='),'>')")
       STOP
    End IF

  END SUBROUTINE file_err

  !============================================================================
  !> Subroutine for opening files with big endian encoding
  !> 
  Subroutine open_as_big_endian_stream(unit,file,action,status,position)

    integer(kind=4) ,intent(in)           :: unit
    character(len=*),intent(in)           :: file,action,status
    character(len=*),intent(in), optional :: position

    integer(kind=4)             :: ier
    character(len=mcl)          :: loc_pos

    If (present(position)) then
       loc_pos = position
    Else
       loc_pos='rewind'
    End If

    !**************************************************************************
    !** GFortran, Intel implementation
    Open(unit=unit, file=trim(file), action=trim(action), &
         status=trim(status), &
         access="stream", convert="big_endian", position=trim(loc_pos))

!!$    !**************************************************************************
!!$    !** CRAY-Fortran implementation
!!$    call asnunit (unit,"-N swap_endian",ier)
!!$    IF (ier /= 0) Then
!!$       WRITE(*,*)
!!$       WRITE(*,"(80('='))")
!!$       WRITE(*,"('EE ',A,T77,' EE')")  'Asign operation on file :'       
!!$       WRITE(*,"('EE ',A          )")  file
!!$       WRITE(*,"('EE ',A,I0,T77,' EE')")   'Connected to unit :',unit
!!$       WRITE(*,"('EE ',A,T77,' EE')")   'faild !!'
!!$       WRITE(*,"('EE ',A,I0,T77,' EE')")'With error flag ',ier
!!$       WRITE(*,"('EE PROGRAM STOPPED ..... ',T77,' EE',/,'<',77('='),'>')")
!!$       STOP
!!$    End IF

!!$    Open(unit=unit, file=trim(file), action=trim(action), &
!!$         status=trim(status), access="stream")

  End Subroutine open_as_big_endian_stream

  !============================================================================
  !> Subroutine for closing files with big endian encoding
  !> 
  Subroutine close_big_endian_stream(unit)

    integer(kind=4) ,intent(in) :: unit
    integer(kind=4)             :: ier

    !**************************************************************************
    !** GFortran, Intel implementation
    CLOSE(unit=unit)

!!$
!!$    !**************************************************************************
!!$    !** CRAY-Fortran implementation
!!$    call asnunit (unit,"-R",ier)
!!$   IF (ier /= 0) Then
!!$       WRITE(*,*)
!!$       WRITE(*,"(80('='))")
!!$       WRITE(*,"('EE ',A,T77,' EE')")  'Asign release on unit :',unit
!!$       WRITE(*,"('EE ',A,T77,' EE')")  'faild !!'
!!$       WRITE(*,"('EE ',A,I0,T77,' EE')")'With error flag ',ier
!!$       WRITE(*,"('EE PROGRAM STOPPED ..... ',T77,' EE',/,'<',77('='),'>')")
!!$       STOP
!!$    End IF

!!$    CLOSE(unit=unit)

  End Subroutine close_big_endian_stream
  !> @}
  !# End of memeber group "vtkio private auxilliary routines" #################

End Module cel_vtkio_module
