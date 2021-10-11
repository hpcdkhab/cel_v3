#!/bin/bash
 
#PBS -l walltime=00:30:00
#PBS -lmppwidth=4096
#PBS -lmppnppn=32
#PBS -M hpcdkhab@hlrs.de
#PBS -m n
#PBS -N CEL
#PBS -o /snx1600/ws7/ws/hpcdkhab-cel-0/cel/cray.out
#PBS -e /snx1600/ws7/ws/hpcdkhab-cel-0/cel/cray.err



PROJ_DIR=/snx1600/ws7/ws/hpcdkhab-cel-0/cel/

export CRAY_OMP_CHECK_AFFINITY=TRUE
#module load perftools

#PREFIX;SIZE_START;SIZE_END;PROCS_START;PROCS_END;PERF_NEW;MV_ALGORITHM_ID;MV_COMMUNICATOR_ID;BENCHMARK_ID

#PROCS_ARRAY=("not used" "1" "2" "4" "8" "12" "16" "18" "24" "32"  "36" "48"  "64" "96" "128" "144" "192" "256" "384" "512" "576" "768" "1024" "1152" "1536" "1728" "2304" )
#XX=("not used"   "24" "48" "96" "192" "288" "384")

bash -x ${PROJ_DIR}/run_experiment_cray.sh  "1_512_192_HORNET_ISEND_COMM_CSR_COO" 4 4 1 17 1 1 1 1

bash -x ${PROJ_DIR}/run_experiment_cray.sh  "1_512_192_HORNET_ISEND_GROUPS_CSR_COO" 4 4 1 17 1 1 2 1

bash -x ${PROJ_DIR}/run_experiment_cray.sh "1_512_192_HORNET_IBCAST_GROUPS_CSR_COO" 4 4 1 17 1 1 3 1

