#!/bin/bash
 
#!/bin/bash
 
#PBS -l walltime=00:120:00
#PBS -lmppwidth=64
#PBS -lmppnppn=32
#PBS -M hpcdkhab@hlrs.de
#PBS -m n
#PBS -N CEL_POISSON_3D
#PBS -o /zhome/academic/HLRS/hlrs/hpcdkhab/myProjects/CEL_V3/CEL_V3/cray_small_1.err
#PBS -e /zhome/academic/HLRS/hlrs/hpcdkhab/myProjects/CEL_V3/CEL_V3/cray_small_1.out

PROJ_DIR=/zhome/academic/HLRS/hlrs/hpcdkhab/myProjects/CEL_V3/CEL_V3

export CRAY_OMP_CHECK_AFFINITY=TRUE

bash -x ./run_experiment_cray.sh

