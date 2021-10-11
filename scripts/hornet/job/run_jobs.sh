#!/bin/sh

LUSTRE_DIR_SED="\/nas_home\/hpcdkhab\/job\/data\/"
PROJ_DIR_SED="\/nas_home\/hpcdkhab\/job\/"
LUSTRE_DIR="/nas_home/hpcdkhab/job/data/"
PROJ_DIR="/nas_home/hpcdkhab/job/"


SIZE_START=1
SIZE_END=1
PROC_START=4
PROC_END=4
CORE_START=2
CORE_END=2
PBS_NODES="2"
PREFIX="CG"
MV_COMM_START=5
MV_COMM_END=5
MV_ALGORITHM_START=1
MV_ALGORITHM_END=1
MV_ALGORITHM_OFF_DIAG_START=3
MV_ALGORITHM_OFF_DIAG_END=3
BENCHMARK_START=5
BENCHMARK_END=5
PRECONDITIONER_START=0
PRECONDITIONER_END=0
MATRIX_SCALER_START=1
MATRIX_SCALER_END=1
TOL_EP="1.0E-16"


RUN_ATTEMPTS=1
RUN_APRUN_PARM_1=" "
RUN_APRUN_PARM_N=" "
SUBMIT_COMM="source"
RUN_CHECK_RESULT=1
write_performance=1
PERF_NEW=1
XX=("not used" "12"  "24" "48" "96" "192" "288" "384" "576")
RUN_XX="(\"not used\" "12"   "24" "48" "96" "192" "288" "384" "576")"

PROCS_ARRAY=("not used" "1" "2" "4" "8" "12" "16" "18" "24" "32"  "36" "48"  "64" "96" "128" "144" "192" "256" "384" "512" "576" "768" "1024" "1152" "1536" "1728" "2304" )
COMM_ARRAY=("not_used" "ISEND" "ISEND_GROUPS" "IBCAST_GROUPS" "ISEND_IDX" "ISEND_COMM_BUFER" )
MV_ALGORITHM_ARRAY=("not used" "MV_CSR" "MV_JAD" "MV_COO" )
BENCHMARK_ARRAY=("not used" "BENCHMARK_MV" "BENCHMARK_R_D" "BENCHMARK_NORM2" "BENCHMARK_CG" "BENCHMARK_GMRES")
PRECONDITIONER_ARRAY=("null" "JACOBI")
MATRIX_SCALER_ARRAY=("MATRIX_SCALER_NULL" "MATRIX_SCALER_DIAGONAL")
QSUB_REJECT_SLEEP=1
MAIN_SLEEP=1




for bench_idx in `seq ${BENCHMARK_START} ${BENCHMARK_END}`
do
for prec_idx in `seq ${PRECONDITIONER_START} ${PRECONDITIONER_END}`
do
for matscal_idx in `seq ${MATRIX_SCALER_START} ${MATRIX_SCALER_END}`
do
for mvalg_idx in `seq ${MV_ALGORITHM_START} ${MV_ALGORITHM_END}`
do
for mvalg_offd_idx in `seq ${MV_ALGORITHM_OFF_DIAG_START} ${MV_ALGORITHM_OFF_DIAG_END}`
do
for k in `seq ${MV_COMM_START} ${MV_COMM_END}`
do
 for j in `seq ${SIZE_START} ${SIZE_END}`
 do
  file_name_part=${PREFIX}_${BENCHMARK_ARRAY[bench_idx]}_${PRECONDITIONER_ARRAY[prec_idx]}_${MATRIX_SCALER_ARRAY[matscal_idx]}_mv_comm_${k}_mv_alg_${MV_ALGORITHM_ARRAY[mvalg_idx]}_${MV_ALGORITHM_ARRAY[mvalg_offd_idx]}_size_${XX[$j]}_procs_${PROCS_ARRAY[${PROC_START}]}_${PROCS_ARRAY[${PROC_END}]}_core_${CORE_START}_${CORE_END}
  concatenate_file="${LUSTRE_DIR}concatenate_${file_name_part}.sh"
  result_file="${LUSTRE_DIR}result_${file_name_part}.dat"
  run_file=run_experiment_cray_${file_name_part}.sh
  cp run_experiment_cray.template ${run_file}

  sed_string="s/USER_RUN_ATTEMPTS/${RUN_ATTEMPTS}/g"
  sed -i "${sed_string}" ${run_file} 
  sed_string="s/USER_RUN_APRUN_PARM_1/${RUN_APRUN_PARM_1}/g"
  sed -i "${sed_string}" ${run_file} 
  sed_string="s/USER_RUN_APRUN_PARM_N/${RUN_APRUN_PARM_N}/g"
  sed -i "${sed_string}" ${run_file} 
  sed_string="s/USER_RUN_XX/${RUN_XX}/g"
  sed -i "${sed_string}" ${run_file} 
  sed_string="s/USER_RUN_CHECK_RESULT/${RUN_CHECK_RESULT}/g"
  sed -i "${sed_string}" ${run_file}
  sed_string="s/USER_TOL_EPS/${TOL_EP}/g"

  sed -i "${sed_string}" ${run_file}


  #sed_string="s/USER_MV_ALGORITHM/${mvalg_idx}/g"
  #sed -i "${sed_string}" ${run_file} 
  sed_string="s/user_write_performance/${write_performance}/g"
  sed -i "${sed_string}" ${run_file} 

  touch ${concatenate_file} 
  touch ${result_file}
  for p in `seq ${PROC_START} ${PROC_END}`
  do 
   for i in `seq ${CORE_START} ${CORE_END}`
   do
     prop=${COMM_ARRAY[k]}_size_${XX[j]}_proc_${PROCS_ARRAY[p]}_core_${i}_for_${file_name_part}

     perf_file="${prop}.dat"

     job_file_name=./job_on_cray_${prop}.sh
     cp job_on_cray.template  ${job_file_name} 
      
     echo "prepare job:"  ${job_file_name}

     let JOB_WIDTH=${PROCS_ARRAY[p]}*8

     sed_string="s/USER_RUN_SCRIPT_NAME/${run_file}/g"
     sed -i "${sed_string}" ${job_file_name} 

     sed_string="s/USER_JOB_WIDTH/${JOB_WIDTH}/g"
     sed -i "${sed_string}" ${job_file_name} 

     
     sed_string="s/USER_PROC_START/${p}/g"
     sed -i "${sed_string}" ${job_file_name} 

     sed_string="s/USER_PBS_NODES/${PBS_NODES}/g"
     sed -i "${sed_string}" ${job_file_name} 


     sed_string="s/USER_MV_ALGORITHM/${mvalg_idx}/g"
     sed -i "${sed_string}" ${job_file_name}
     
     sed_string="s/user_matrix_scaler/${matscal_idx}/g"
     sed -i "${sed_string}" ${job_file_name}
 
     sed_string="s/user_mv_algorithm_off_diag/${mvalg_offd_idx}/g"
     sed -i "${sed_string}" ${job_file_name} 

     sed_string="s/USER_PROC_END/${p}/g"
     sed -i "${sed_string}" ${job_file_name} 

     sed_string="s/USER_CORE_START/${i}/g"
     sed -i "${sed_string}" ${job_file_name} 

     sed_string="s/USER_CORE_END/${i}/g"
     sed -i "${sed_string}" ${job_file_name} 

     sed_string="s/USER_SIZE_START/${j}/g"
     sed -i "${sed_string}" ${job_file_name} 

     sed_string="s/USER_SIZE_END/${j}/g"
     sed -i "${sed_string}" ${job_file_name} 

     sed_string="s/USER_PREFIX/${PREFIX}/g"
     sed -i "${sed_string}" ${job_file_name} 
 
     sed_string="s/USER_PERF_NEW/${PERF_NEW}/g"
     sed -i "${sed_string}" ${job_file_name} 

     sed_string="s/USER_PROJ_DIR/${PROJ_DIR_SED}/g"
     sed -i "${sed_string}" ${job_file_name} 

     sed_string="s/USER_LUSTRE_DIR/${LUSTRE_DIR_SED}/g"
     sed -i "${sed_string}" ${job_file_name} 

     sed_string="s/USER_OUT/${prop}/g"
     sed -i "${sed_string}" ${job_file_name} 
     
     sed_string="s/USER_BENCHMARK/${bench_idx}/g"
     sed -i "${sed_string}" ${job_file_name}

     sed_string="s/user_preconditioner/${prec_idx}/g"
     sed -i "${sed_string}" ${job_file_name}

     sed_string="s/USER_ERR/${prop}/g"
     sed -i "${sed_string}" ${job_file_name} 
 
     sed_string="s/USER_JOB_NAME/${prop}/g"
     sed -i "${sed_string}" ${job_file_name} 
 
     sed_string="s/USER_MV_COMMUNICATOR/${k}/g"
     sed -i "${sed_string}" ${job_file_name} 
 
     sed_string="s/USER_PERF_FILE/${perf_file}/g"
     sed -i "${sed_string}" ${job_file_name} 
 

    ${SUBMIT_COMM} ${job_file_name} 
    while [ $? -ne 0 ]; do
       echo "could not submit the job" ${job_file_name}
       echo "will try in ${QSUB_REJECT_SLEEP} minutes"
       sleep ${QSUB_REJECT_SLEEP}
       ${SUBMIT_COMM} ${job_file_name} 
     done
    echo "job submitted: " ${job_file_name}
    echo "cat ${LUSTRE_DIR}${perf_file} >> ${result_file}" >> ${concatenate_file}  
    sleep 1
   done
   echo "sleep ${MAIN_SLEEP} seconds"
   sleep ${MAIN_SLEEP}
  done
 done

done
done
done
done
done
done
