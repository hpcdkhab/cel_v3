#!/bin/bash

#PBS -l nodes=1:mem32gb
#PBS -l walltime=00:35:00
#PBS -M khabi@hlrs.de
#PBS -m abe
#PBS -N cel_mpi_omp_cg

export PROJ_DIR=/zhome/academic/HLRS/hlrs/hpcdkhab/myProjects/CEL_CRESTA/cresta/wp4/CEL_V2/
export EXEC=$PROJ_DIR/tests/cel_mpi_omp/cel_mpi_omp_cg
export IO_DIR=/lustre/ws1/ws/hpcdkhab-cel-0
module load compiler/intel/13.0.1
module load mpi/impi/4.1.0-intel-13.0.0
cd $IO_DIR
mpiexec -trace  -n 2  $EXEC  -attempts 1 -verbosity 1 -nx 32 -ny 32 -nz 32 -dx 2 -dy 1 -dz 1 -tests_num 1 -threads 8 -generatematrix 1 -num_chunks 2 -chunk_size 512 -iter_max 1000 -eps 0.1
