#!/bin/bash
 
#!/bin/bash
 
#PBS -l walltime=00:180:00
#PBS -lmppwidth=1024
#PBS -lmppnppn=32
#PBS -M hpcdkhab@hlrs.de
#PBS -m n
#PBS -N CRESTA_CEL_POISSON_3D
#PBS -o /snx1600/ws7/ws/hpcdkhab-cel-0/cel/cray.err
#PBS -e /snx1600/ws7/ws/hpcdkhab-cel-0/cel/cray.out

PROJ_DIR=/snx1600/ws7/ws/hpcdkhab-cel-0/cel

export CRAY_OMP_CHECK_AFFINITY=TRUE

bash -x ${PROJ_DIR}/run_experiment_cray.sh 

