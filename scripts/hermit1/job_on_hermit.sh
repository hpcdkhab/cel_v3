#!/bin/bash

#PBS -l walltime=00:15:00
#PBS -lmppwidth=256
#PBS -lmppnppn=32
#PBS -M hpcdkhab@hlrs.de
#PBS -m abe
#PBS -N HERMIT_CEL
#PBS -o /univ_1/ws1/ws/hpcdkhab-cel-0/cel/
#PBS -e /univ_1/ws1/ws/hpcdkhab-cel-0/cel/


export PROJ_DIR=//univ_1/ws1/ws/hpcdkhab-cel-0/cel/

export CRAY_OMP_CHECK_AFFINITY=TRUE

#PREFIX;SIZE_START;SIZE_END;PROCS_START;PROCS_END;PERF_NEW;MV_ALGORITHM_ID;MV_COMMUNICATOR_ID;BENCHMARK_ID



#bash -x ${PROJ_DIR}/run_experiment_cray_hermit.sh "HERMIT_ISEND_COMM_CSR_COO" 1 4 1 19 1 1 1 1

#bash -x ${PROJ_DIR}/run_experiment_cray_hermit.sh "HERMIT_ISEND_GROUPS_CSR_COO" 1 4 1 19 1 1 2 1

bash -x ${PROJ_DIR}/run_experiment_cray_hermit.sh "3_HERMIT_IBCAST_GROUPS_CSR_COO" 1 1 1 6 1 1 3 1

