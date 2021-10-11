#!/bin/bash
 
#PBS -l nodes=32:mem32gb:ppn=8
#PBS -l walltime=00:60:00
#PBS -M khabi@hlrs.de
#PBS -m abe
#PBS -N CRESTA_BENCHMARK_PETSC_3D_256x256x256

PROJ_DIR=/zhome/academic/HLRS/hlrs/hpcdkhab/myProjects/cel
cd $PROJ_DIR
module load compiler/intel/13.1.0
module load mpi/impi/4.1.0-intel-13.0.0
module load numlib/intel/mkl/11.0

mpirun -n 8 $PROJ_DIR/tests/cel_test_petsc -x 256 -y 256 -z 256 -a 2 -b 2 -c 2 -v 0 -s 1 -o 0
mpirun -n 8 $PROJ_DIR/tests/cel_test_petsc -x 256 -y 256 -z 256 -a 2 -b 2 -c 2 -v 0 -s 2 -o 0
mpirun -n 8 $PROJ_DIR/tests/cel_test_petsc -x 256 -y 256 -z 256 -a 2 -b 2 -c 2 -v 0 -s 3 -o 0

mpirun -n 16 $PROJ_DIR/tests/cel_test_petsc -x 256 -y 256 -z 256 -a 4 -b 2 -c 2 -v 0 -s 1 -o 0
mpirun -n 16 $PROJ_DIR/tests/cel_test_petsc -x 256 -y 256 -z 256 -a 4 -b 2 -c 2 -v 0 -s 2 -o 0
mpirun -n 16 $PROJ_DIR/tests/cel_test_petsc -x 256 -y 256 -z 256 -a 4 -b 2 -c 2 -v 0 -s 3 -o 0

mpirun -n 32 $PROJ_DIR/tests/cel_test_petsc -x 256 -y 256 -z 256 -a 4 -b 4 -c 2 -v 0 -s 1 -o 0
mpirun -n 32 $PROJ_DIR/tests/cel_test_petsc -x 256 -y 256 -z 256 -a 4 -b 4 -c 2 -v 0 -s 2 -o 0
mpirun -n 32 $PROJ_DIR/tests/cel_test_petsc -x 256 -y 256 -z 256 -a 4 -b 4 -c 2 -v 0 -s 3 -o 0

mpirun -n 64 $PROJ_DIR/tests/cel_test_petsc -x 256 -y 256 -z 256 -a 4 -b 4 -c 4 -v 0 -s 1 -o 0
mpirun -n 64 $PROJ_DIR/tests/cel_test_petsc -x 256 -y 256 -z 256 -a 4 -b 4 -c 4 -v 0 -s 2 -o 0
mpirun -n 64 $PROJ_DIR/tests/cel_test_petsc -x 256 -y 256 -z 256 -a 4 -b 4 -c 4 -v 0 -s 3 -o 0

mpirun -n 128 $PROJ_DIR/tests/cel_test_petsc -x 256 -y 256 -z 256 -a 8 -b 4 -c 4 -v 0 -s 1 -o 0
mpirun -n 128 $PROJ_DIR/tests/cel_test_petsc -x 256 -y 256 -z 256 -a 8 -b 4 -c 4 -v 0 -s 2 -o 0
mpirun -n 128 $PROJ_DIR/tests/cel_test_petsc -x 256 -y 256 -z 256 -a 8 -b 4 -c 4 -v 0 -s 3 -o 0

mpirun -n 256 $PROJ_DIR/tests/cel_test_petsc -x 256 -y 256 -z 256 -a 8 -b 8 -c 4 -v 0 -s 1 -o 0
mpirun -n 256 $PROJ_DIR/tests/cel_test_petsc -x 256 -y 256 -z 256 -a 8 -b 8 -c 4 -v 0 -s 2 -o 0
mpirun -n 256 $PROJ_DIR/tests/cel_test_petsc -x 256 -y 256 -z 256 -a 8 -b 8 -c 4 -v 0 -s 3 -o 0

mpirun -n 512 $PROJ_DIR/tests/cel_test_petsc -x 256 -y 256 -z 256 -a 8 -b 8 -c 8 -v 0 -s 1 -o 0
mpirun -n 512 $PROJ_DIR/tests/cel_test_petsc -x 256 -y 256 -z 256 -a 8 -b 8 -c 8 -v 0 -s 2 -o 0
mpirun -n 512 $PROJ_DIR/tests/cel_test_petsc -x 256 -y 256 -z 256 -a 8 -b 8 -c 8 -v 0 -s 3 -o 0
