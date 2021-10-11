!
! File:   cel_omp_mv_generate_matrix_module.f90
! Project CRESTA (see details on https://cresta-project.eu) Exascale library
! Send you email about the bugs to
! The High Performance Computing Center Stuttgart (HLRS) of the University of Stuttgart
! Dmitry Khabi,  (email: khabi@hlrs.de)
! Created on May 29, 2013
!

module cel_sp_mat_generate_module
use cel_types_module
use cel_error_module
use cel_domain_3d_module
use cel_vec_module
use cel_sp_mat_module
use cel_c_sp_mat_interface_module
use cel_discretisation_task_module
use cel_input_parameters_module
use cel_base_print_module
use cel_mmio_module
!external libraries
use OMP_LIB
use, intrinsic :: iso_c_binding
implicit none


contains

subroutine cel_sp_mat_generate(sp_mat, &
  boundary_vec, th_solution_vec, &
  input_parameters, omp, test_idx, err_code, block_indx)
  type(cel_sp_mat_type), intent(inout) :: sp_mat
  type(cel_vec_type), intent(inout) :: boundary_vec
  type(cel_vec_type), intent(inout) :: th_solution_vec
  type(cel_input_parameters_type), intent(in) :: input_parameters
  type(cel_omp_type), intent(in) :: omp
  integer(kind=ik), intent(in) :: test_idx
  integer(kind=ik), intent(inout) :: err_code
  integer(kind=ik), optional, intent(in) :: block_indx
  type(cel_c_discretisation_task_type) :: discretisation_task_c
  type(cel_c_domain_3D_type) :: domain_3D_c
  type(cel_domain_3D_type) :: domain_3D
  type(cel_c_sp_mat_type) :: sp_mat_c
  type(cel_c_vec_type) :: boundary_vec_c
  type(cel_c_vec_type) :: th_solution_vec_c
  real(kind=rk) :: temp_real
  integer(kind=rk) :: temp_int
  
  if(present(block_indx)) then
    domain_3D_c%block_indx=block_indx
  else
    domain_3D_c%block_indx=omp%master_num
  endif
  temp_real = REAL(input_parameters%nx, kind=rk)*(input_parameters%increase_factor-1.0_rk)*(test_idx-1_ik)
  temp_int = INT(temp_real,kind=ik)
  domain_3D_c%nx=input_parameters%nx+temp_int
  domain_3D_c%nx=domain_3D_c%nx/input_parameters%dx
  domain_3D_c%nx=domain_3D_c%nx*input_parameters%dx
  temp_real = REAL(input_parameters%ny, kind=rk)*(input_parameters%increase_factor-1.0_rk)*(test_idx-1_ik)
  temp_int = INT(temp_real,kind=ik)
  domain_3D_c%ny=input_parameters%ny+temp_int
  domain_3D_c%ny=domain_3D_c%ny/input_parameters%dy
  domain_3D_c%ny=domain_3D_c%ny*input_parameters%dy
  temp_real = REAL(input_parameters%nz, kind=rk)*(input_parameters%increase_factor-1.0_rk)*(test_idx-1_ik)
  temp_int = INT(temp_real,kind=ik)
  domain_3D_c%nz=input_parameters%nz+temp_int
  domain_3D_c%nz=domain_3D_c%nz/input_parameters%dz
  domain_3D_c%nz=domain_3D_c%nz*input_parameters%dz
  domain_3D_c%dx=input_parameters%dx
  domain_3D_c%dy=input_parameters%dy
  domain_3D_c%dz=input_parameters%dz
  discretisation_task_c%problem = POISSON_3D
  if(domain_3D_c%ny .eq. domain_3D_c%nx .and. domain_3D_c%ny .eq. domain_3D_c%nz) then
    discretisation_task_c%stencil = STENCIL_27
  else
     discretisation_task_c%stencil = STENCIL_7
  end if
  discretisation_task_c%boundary_cond = X_2_PLUS_Y_2_PLUS_Z_2
  
  if(domain_3D_c%nx + domain_3D_c%ny + domain_3D_c%nz .le. 2_ik) then
    err_code = 1_ik
    call cel_error("Error in cel_sp_mat_generate:check the parameters: too many processe for the matrix ?", &
    err_code, .TRUE._lk,.TRUE._lk)
  end if
  
  err_code =  INT(cel_c_sp_mat_interface &
    (INT(domain_3D_c%block_indx,kind=c_int), discretisation_task_c, domain_3D_c, sp_mat_c, &
    boundary_vec_c, th_solution_vec_c,  INT(1, kind=ceik) ), kind=ik)

  call cel_error("Error in cel_sp_mat_generate:cel_c_sp_mat_interface", &
    err_code, .TRUE._lk,.TRUE._lk)
  if(err_code .EQ. 0_ik) then

    call cel_convert_c(domain_3D_c,domain_3D, .TRUE._lk)
    call cel_convert_c(sp_mat_c,sp_mat, .TRUE._lk)
    call cel_convert_c(boundary_vec_c,boundary_vec, .TRUE._lk)
    call cel_convert_c(th_solution_vec_c,th_solution_vec, .TRUE._lk)

    sp_mat%rows=sp_mat%rows+1_ik
    sp_mat%cols=sp_mat%cols+1_ik
    sp_mat%block_rows_num=1_ik
    sp_mat%block_cols_num=1_ik
  end if
  
  
end subroutine cel_sp_mat_generate



subroutine cel_sp_mat_generate_simple(sparse_matrix_coo, &
  boundary_vector, th_solution_vector, &
  input_parameters, omp, test_idx, err_code)
  type(cel_sp_mat_type), intent(inout) :: sparse_matrix_coo
  type(cel_vec_type), intent(inout) :: boundary_vector
  type(cel_vec_type), intent(inout) :: th_solution_vector
  type(cel_input_parameters_type), intent(in) :: input_parameters
  type(cel_omp_type), intent(in) :: omp
  integer(kind=ik), intent(in) :: test_idx
  integer(kind=ik), intent(inout) :: err_code
  integer(kind=ik) :: rows_num, mmax, nn, mm, idx
  logical (kind=lk) :: output_on
  real(kind=rk) :: temp_real
  integer(kind=rk) :: temp_int
  
  temp_real = REAL(input_parameters%nx, kind=rk)*(input_parameters%increase_factor-1.0_rk)*(test_idx-1_ik)
  temp_int = input_parameters%nx+INT(temp_real,kind=ik)
  temp_real = REAL(input_parameters%nx, kind=rk)*(input_parameters%increase_factor-1.0_rk)*(test_idx-1_ik)
  temp_int = temp_int+input_parameters%ny+INT(temp_real,kind=ik)
  temp_real = REAL(input_parameters%nx, kind=rk)*(input_parameters%increase_factor-1.0_rk)*(test_idx-1_ik)
  temp_int = temp_int+input_parameters%nz+INT(temp_real,kind=ik)
  rows_num = temp_int
  mmax = 27

  call cel_sp_mat_new(sparse_matrix_coo, rows_num, rows_num,&
  1_ik, 1_ik ,mmax*rows_num)
  sparse_matrix_coo%block_indx=omp%master_num
  
  idx = 1
  do mm=1,mmax
    do nn=1,rows_num
      sparse_matrix_coo%rows(idx) = nn
      sparse_matrix_coo%cols(idx) = min(max(1,mm+nn-(mmax-1)/2),rows_num)    ! not to exceed [1,nmax]
      sparse_matrix_coo%values(idx) = REAL(nn+mm,kind=rk)
      idx = idx + 1
    enddo
  enddo

  call cel_sp_mat_format_sparse_31_sort(sparse_matrix_coo, output_on, err_code)
  
  call new(boundary_vector,rows_num)
  call new(th_solution_vector,rows_num)
  boundary_vector%values=1.0
  th_solution_vector%values=1.0
  
end subroutine cel_sp_mat_generate_simple

subroutine cel_sp_mat_generate_simple_test(sparse_matrix_coo, &
  boundary_vector, th_solution_vector, &
  input_parameters, omp, test_idx, err_code)
  type(cel_sp_mat_type), intent(inout) :: sparse_matrix_coo
  type(cel_vec_type), intent(inout) :: boundary_vector
  type(cel_vec_type), intent(inout) :: th_solution_vector
  type(cel_input_parameters_type), intent(in) :: input_parameters
  type(cel_omp_type), intent(in) :: omp
  integer(kind=ik), intent(in) :: test_idx
  integer(kind=ik), intent(inout) :: err_code
  integer(kind=ik) :: ii,jj,kk, idx
  integer(kind=ik) :: num_local_rows, num_rows, block_size, num_local_blocks, num_local_raw_blocks
  integer(kind=ik) :: offset_local_raw_blocks, mult_factor
  
  err_code = 0_ik
  if(omp%is_master) then
    num_rows = input_parameters%nx*input_parameters%ny*input_parameters%nz
    num_local_rows = num_rows / omp%num_masters
    if(.not. (num_rows .eq. num_local_rows*omp%num_masters)) then
      call cel_error("cel_sp_mat_generate_simple_test Error check input parameter (nx,ny,nz)", &
      1_ik, .TRUE._lk,.TRUE._lk)
    end if
    block_size = input_parameters%num_chunks*input_parameters%chunk_size
    num_local_raw_blocks = num_local_rows / block_size
    if(.not. (num_local_rows .eq. num_local_raw_blocks*block_size)) then
      call cel_error("cel_sp_mat_generate_simple_test Error check input parameter (chunk_size,num_chunks)", &
      2_ik, .TRUE._lk,.TRUE._lk)
    end if
    num_local_blocks = num_local_raw_blocks*1_ik
    offset_local_raw_blocks = omp%master_num*num_local_raw_blocks
    
    mult_factor = 2_ik
    if(omp%master_num .EQ. 0) mult_factor = 1_ik
    if(omp%master_num .EQ. omp%num_masters-1_ik) mult_factor = 1_ik
    
    call cel_sp_mat_new(sparse_matrix_coo, num_local_rows, num_rows,&
    1_ik, 1_ik ,mult_factor*block_size*block_size*num_local_blocks)
    
    sparse_matrix_coo%block_indx=omp%master_num
    call new(boundary_vector,num_local_rows)
    call new(th_solution_vector,num_local_rows)
    th_solution_vector%values=1.0
    idx = 1_ik
   
    do ii=1,num_local_raw_blocks
      do jj=1,block_size
        do kk=1,block_size*mult_factor
        sparse_matrix_coo%rows(idx) = (ii+offset_local_raw_blocks-1_ik)*block_size+jj
        sparse_matrix_coo%cols(idx) = (ii+offset_local_raw_blocks-1_ik)*block_size+kk
        sparse_matrix_coo%values(idx) = REAL(ii*(omp%master_num+1),kind=rk)
        idx = idx + 1_ik
        enddo
        boundary_vector%values((ii-1)*block_size+jj)=REAL(ii*(omp%master_num+1)*block_size*mult_factor,kind=rk)
      enddo
    enddo
    !call cel_sp_mat_format_sparse_31_sort(sparse_matrix_coo, .FALSE._lk, err_code)
  end if

end subroutine cel_sp_mat_generate_simple_test

  

subroutine cel_sp_mat_read_matrix_market_file(sp_mat, matrix_dir,&
  filename,&
  block_rows_num,&
  block_cols_num,&
  omp, &
  err_code)
  type(cel_sp_mat_type), intent(inout)  :: sp_mat
  character(len=*), intent(in)         :: matrix_dir
  character(len=*), intent(in)         :: filename
  integer(kind=ik)                     :: block_rows_num
  integer(kind=ik)                     :: block_cols_num
  type(cel_omp_type), intent(in)       :: omp
  integer (kind=ik), intent(out)          :: err_code
  character(len=1024)                  :: path_to_matrix
  integer(kind=ik)                     :: iunit, io_stat
  integer(kind=ik)                     :: nrows, ncols, nnz
  
  character(len=7)                   :: field
  character(len=19)                  :: symm
  character(len=10)                  :: rep
  
  integer(kind=ik) :: ival(1)
  complex cval(1)
  
  err_code = 0_ik
      
  if(omp%is_master) then

    write(path_to_matrix,'(A,A,A)') trim(matrix_dir),'/',trim(filename)
    write(*,'(A)') trim(path_to_matrix)
    iunit = 1027
    open(unit=iunit,file=path_to_matrix,iostat = io_stat,access='sequential',&
     action='read')
    call cel_file_err(path_to_matrix, io_stat)
    
    call mminfo(iunit,rep,field,symm,nrows,ncols,nnz)
    if(cel_is_in_debug) then
      print *,' Matrix is type: ',rep,' ',field,' ',symm
      write(*,'(A,I0,A,I0,A,I0,A)')'  Matrix size: ',nrows,' by ',ncols,' with ', &
            nnz,' nonzeros.'
     !write(*,'(A,I0,A,I0)')'  block_rows_num: ',block_rows_num,'  block_cols_num ', block_cols_num
    end if
    
    call new(sp_mat, nrows, ncols,&
      block_rows_num, block_cols_num ,nnz)
    
    call mmread(iunit,rep,field,symm,nrows,ncols,nnz,huge(1_ik),&
                  sp_mat%rows,sp_mat%cols,ival,sp_mat%values,cval)  
    !call print("P_MATRIX",sp_mat)
    
    close(iunit,iostat = io_stat)
    call cel_file_err(path_to_matrix, io_stat)
  end if

end subroutine cel_sp_mat_read_matrix_market_file

end module  cel_sp_mat_generate_module
