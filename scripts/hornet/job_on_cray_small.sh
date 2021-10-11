#!/bin/bash
 
#PBS -l walltime=00:15:00
#PBS -lmppwidth=2048
#PBS -lmppnppn=32
#PBS -M hpcdkhab@hlrs.de
#PBS -m n
#PBS -N CEL_POISSON_3D_1024
#PBS -o /snx1600/ws7/ws/hpcdkhab-cel-0/cel/cray_small.out
#PBS -e /snx1600/ws7/ws/hpcdkhab-cel-0/cel/cray_small.err

PROJ_DIR=/snx1600/ws7/ws/hpcdkhab-cel-0/cel/

export CRAY_OMP_CHECK_AFFINITY=TRUE
#module load perftools

bash -x ${PROJ_DIR}/run_experiment_cray_small.sh "IBCAST_CSR_COO_SMALL" 3 3 8 8 1 1 1 1

