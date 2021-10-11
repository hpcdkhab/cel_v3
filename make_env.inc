#     
# File:  Makefile
# Project CRESTA (see details on https://cresta-project.eu) Exascale library
# Send you email about the bugs to
# The High Performance Computing Center Stuttgart (HLRS) of the University of Stuttgart
# Dmitry Khabi,  (email: khabi@hlrs.de)
# Created on March 25, 2013
#

#YES; NO
export IS_DEBUG=NO
#GNU; INTEL; CRAY;
export COMPILER=CRAY
export MPI=
#YES;NO
export CRAY=YES
#NO; VAMPIR
PROFILER=NO

export COMPILER_DIR=$(shell pwd)
export DIRECTIVES=-D_CEL_PERF

export PERF_TOOL=
ifndef IS_DEBUG
$(error env. variable IS_DEBUG in make_env.inc is not set)
endif
ifndef COMPILER
$(error env. variable COMPILER in make_env.inc is not set )
endif

ifeq ($(COMPILER),INTEL)
  ifeq ($(CRAY),YES)
  # INTEL compiler
   export CC=cc
   export MPI_CC = cc
   export F90 = ftn
   export MPI_F90 = ftn
   export F90_MODULE_DIR = -module 
  else
  # INTEL compiler
   export CC = icc
   export F90 = ifort
   ifeq ($(MPI),OPENMPI)
     export MPI_CC = mpicc
     export MPI_F90 = mpifort
  else
     export MPI_CC = mpiicc
     export MPI_F90 = mpif90
   endif
   export F90_MODULE_DIR = -module
  endif
  #intel mkl includes flags  so as to satisfy the circular dependencies
  export INCLUDE_BEGIN_GROUP =-Wl,--start-group
  export INCLUDE_END_GROUP =-Wl,--end-group
  ifeq ($(IS_DEBUG),YES)
  # DEBUG
    export CFLAGS = -O0 -restrict -vec-report1   -g -traceback -C -debug full   $(DIRECTIVES)
    export FFLAGS =  -warn all -O0  -c  -traceback -C -debug full -assume realloc_lhs -heap-arrays 50 -openmp  -check bounds -g  $(DIRECTIVES) 
    export FLFLAGS  =  -warn all -O0  -openmp -check bounds  $(DIRECTIVES)
  else
  # RELEASE
    export CFLAGS = -O3 -openmp -restrict   -ansi-alias -unroll-aggressive -g  $(DIRECTIVES)
    export FFLAGS =    -O3 -openmp -warn  all -c -assume realloc_lhs -heap-arrays 50 -g  $(DIRECTIVES) 
    #export FFLAGS =    -O3 -openmp -warn  all - -g  $(DIRECTIVES)
    export FLFLAGS  =   -O3  -warn all -openmp  $(DIRECTIVES)

  #flags=             -c   -check:all  -traceback -C -debug full -assume realloc_lhs 
  #flags=             -c   -check:all  -traceback -C -debug full -assume realloc_lhs -heap-arrays 50
  #flags=   -c -assume realloc_lhs -O3 -heap-arrays 50 -openmp
  #flags= -warn all   -c   -check:all  -traceback -C -debug full -assume realloc_lhs -heap-arrays 50 -openmp
  #flags=   -c -openmp

  endif
  
endif

ifeq ($(COMPILER),GNU)

  ifeq ($(CRAY),YES)
   #GNU compiler
   export CC=cc
   export MPI_CC = cc
   export F90 = ftn
   export MPI_F90 = ftn
   export F90_MODULE_DIR = -J
  else
   #GNU compiler
   export CC=gcc
   export MPI_CC = mpicc
   export F90 = gfortran
   export MPI_F90 = mpif90
   export F90_MODULE_DIR = -J
  endif
  ifeq ($(IS_DEBUG),YES)
  # DEBUG
    export CFLAGS = -O0 -g -Wall -Wno-unused-variable -fstrict-aliasing  -mavx -fopenmp  $(DIRECTIVES)
    export FFLAGS =  -O0 -g -Wall -Wno-unused-variable -mavx  -fopenmp  -fbounds-check  $(DIRECTIVES)
    export FLFLAGS  =  -O0 -g -Wall -Wno-unused-variable -fopenmp -fbounds-check  $(DIRECTIVES)
  else
  # RELEASE
  
    export CFLAGS = -O3 -g -Wall -Wno-unused-variable -ffast-math -fstrict-aliasing -funroll-loops -mavx -fopenmp   $(DIRECTIVES)
    export FFLAGS =  -O3 -g -Wall -Wno-unused-variable -ffast-math -fstrict-aliasing -funroll-loops  -mavx  -fopenmp -fbounds-check $(DIRECTIVES)
    export FLFLAGS  =  -O3 -g -Wall -fopenmp  $(DIRECTIVES)

  endif
  
endif

ifeq ($(COMPILER),CRAY)
  export CC=cc
  export MPI_CC = cc
  export F90 = ftn
  export MPI_F90 = ftn
  export F90_MODULE_DIR = -J

  ifeq ($(IS_DEBUG),YES)
  # DEBUG
    export CFLAGS = -O0 -g -DCRAY_COMPILER -h omp   $(DIRECTIVES)  
    export FFLAGS = -O0 -g -em -h omp   $(DIRECTIVES) -R bcdps -rm
    export FLFLAGS  = -O0 -g   $(DIRECTIVES)
  else
  # RELEASE
    export CFLAGS =  -DCRAY_COMPILER -h display_opt   $(DIRECTIVES)
    #export FFLAGS = -O3 -DCRAY_COMPILER -h display_opt -em -h omp -hpl=$(COMPILER_DIR)/cray.pl -r m  $(DIRECTIVES) -h autoprefetch -h noautothread -h page_align_allocate -h scalar3 -h thread3 -h unroll2 -h vector3 -h wp  -h cache3 -h fusion2
    export FFLAGS = -O3 -DCRAY_COMPILER -h display_opt -em -hpl=$(COMPILER_DIR)/cray.pl -h omp -r m $(DIRECTIVES) -h autoprefetch -h noautothread  -h scalar3 -h thread3 -h unroll2 -h vector3  -h cache3 -h fusion2 -h wp
    #export FLFLAGS  =  -h omp -h wp -hpl=$(COMPILER_DIR)/cray.pl  $(DIRECTIVES) 
    export FLFLAGS  =   -h omp -h wp -hpl=$(COMPILER_DIR)/cray.pl
  endif

endif

ifeq ($(COMPILER),PGI)
  ifeq ($(CRAY),YES)
  # PGI compiler
   export CC=cc
   export MPI_CC = cc
   export F90 = ftn
   export MPI_F90 = ftn
   export F90_MODULE_DIR = -module 
  else
  # PGI compiler
     echo "PGI compiler is only defined for CRAY"
     exit
  endif
  ifeq ($(IS_DEBUG),YES)
  # DEBUG
    export CFLAGS = -O0 -g $(DIRECTIVES)
    export FFLAGS =  -O0 -g  $(DIRECTIVES) 
    export FLFLAGS  =  -O0   $(DIRECTIVES)
  else
  # RELEASE
    export CFLAGS = -O3  -g  $(DIRECTIVES)
    export FFLAGS =    -O3  -g  $(DIRECTIVES) 
    #export FFLAGS =    -O3  -g  $(DIRECTIVES)
    export FLFLAGS  =   -O3   $(DIRECTIVES)

  endif
  
endif


ifeq ($(PROFILER),VAMPIR)
  export CC=vtcc
  export MPI_CC =vtcc
  export F90 = vtf90
  export MPI_F90 = vtf90
  export F90_MODULE_DIR = -J
endif
