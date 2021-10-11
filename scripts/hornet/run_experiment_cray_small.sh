#!/bin/sh

PROJ_DIR="/snx1600/ws7/ws/hpcdkhab-cel-0/cel/"
BIN="${PROJ_DIR}/cel_cg"
BIN_PAT="${PROJ_DIR}/cel_cg+pat"


APRUN_PARM_1="-ss -j 1 -N1 -cc cpu -d 8"
APRUN_PARM_N="-ss -j 1 -N2 -cc cpu -d 8"

XX=("not used"  "12" "24" "48" "96" "192" "288" "384")

PROCS_ARRAY=("not used" "1" "2" "4" "8" "12" "16" "24" "27" "36" "48"  "64" "96" "128" "144" "192" "256" "384" "512")
ATTEMPTS=1000
DX=("not used" "1" "2" "2" "2" "3"  "4"  "3"  "3"  "4"  "4"   "4"  "6"  "8"   "6"   "8"   "8"   "8"   "8"  )
DY=("not used" "1" "1" "2" "2" "2"  "2"  "3"  "3"  "3"  "4"   "4"  "4"  "4"   "6"   "6"   "8"   "8"   "8"  )
DZ=("not used" "1" "1" "1" "2" "2"  "2"  "2"  "3"  "3"  "3"   "4"  "4"  "4"   "4"   "4"   "4"   "6"   "8"  )

PREFIX=$1
SIZE_START=$2
SIZE_END=$3
PROCS_START=$4
PROCS_END=$5
PERF_NEW=$6


MV_ALGORITHM_ID=$7
MV_COMMUNICATOR_ID=$8
BENCHMARK_ID=$9

echo size: ${XX[${SIZE_START}]} - ${XX[${SIZE_END}]} 
echo procs: ${PROCS_ARRAY[${PROCS_START}]} - ${PROCS_ARRAY[${PROCS_END}]} 
cd ${PROJ_DIR}
export CRAY_OMP_CHECK_AFFINITY=TRUE
export MPICH_RANK_REORDER_METHOD=1
export PAT_RT_SUMMARY=1

perf_file="${PROJ_DIR}/${PREFIX}_poisson_3D_crayxc30_mv_size_${XX[${SIZE_START}]}-${XX[${SIZE_END}]}_procs_${PROCS_ARRAY[${PROCS_START}]}-${PROCS_ARRAY[${PROCS_END}]}.dat"
if [ ${PERF_NEW} -ge 1 ];  then
  rm ${perf_file} -rf
fi

for j in `seq ${SIZE_START} ${SIZE_END}`
do
   for p in `seq ${PROCS_START} ${PROCS_END}`
  do 
   for i in 8
   do
   echo "procs:" ${PROCS_ARRAY[p]}  "threads per proc:" $i "test id:" ${j}
   echo "DX:" ${DX[p]} "DY:"  ${DY[p]} "DZ:" ${DZ[p]} "XX:" ${XX[j]} "ATTEMPTS:" ${ATTEMPTS} 
  export MPICH_RANK_REORDER_METHOD=1
#  export PAT_RT_EXPFILE_DIR=${PROJ_DIR}/craypat_revine_study_${XX[j]}_${PROCS_ARRAY[p]}/
#  export PAT_RT_EXPFILE_NAME=pat_statistic_${XX[j]}_${PROCS_ARRAY[p]}
#  export PAT_RT_HWPC "PAPI_FP_OPS"
#  rm ${PAT_RT_EXPFILE_DIR} -rf
#  mkdir ${PAT_RT_EXPFILE_DIR} -p
# case ${PROCS_ARRAY[p]} in
#  1)
#    aprun -n ${PROCS_ARRAY[p]} ${APRUN_PARM_1} ${BIN_PAT}   -verbosity 1 -nx ${XX[j]} -ny ${XX[j]} -nz ${XX[j]} -dx ${DX[p]} -dy ${DY[p]} -dz ${DZ[p]} -prefetcher 0 -increase_factor 1.1 -tests_num 1 -generatematrix 1  -iter_max 1000  -eps 0.001 -attempts ${ATTEMPTS} -check_result 0 -write_matrix 0 -write_performance  1 -performance_filename temp.file  -threads ${i}
#  ;;
#  *) 
#    aprun -n ${PROCS_ARRAY[p]} ${APRUN_PARM_N} ${BIN_PAT}  -verbosity 1 -nx ${XX[j]} -ny ${XX[j]} -nz ${XX[j]} -dx ${DX[p]} -dy ${DY[p]} -dz ${DZ[p]} -prefetcher 0 -increase_factor 1.1 -tests_num 1 -generatematrix 1  -iter_max 1000  -eps 0.001 -attempts ${ATTEMPTS} -check_result 0 -write_matrix 0 -write_performance  0 -performance_filename temp.file  -threads ${i}
# ;;
# esac 
# sleep 1
# cd ${PAT_RT_EXPFILE_DIR}
#  pat_report  -o report.txt ${PAT_RT_EXPFILE_DIR}/${PAT_RT_EXPFILE_NAME}.xf
#  cd ${PROJ_DIR}
#  rm MPICH_RANK_ORDER -f
#  cp ${PAT_RT_EXPFILE_DIR}/MPICH_RANK_ORDER.*  ${PROJ_DIR}/MPICH_RANK_ORDER
#  sleep 1
#   export MPICH_RANK_REORDER_METHOD=1
#  if [ -e "${PROJ_DIR}/MPICH_RANK_ORDER" ] 
#  then
#   export MPICH_RANK_REORDER_METHOD=3
#  fi
   case ${PROCS_ARRAY[p]} in
    1)
      aprun -n ${PROCS_ARRAY[p]} ${APRUN_PARM_1} ${BIN}   -verbosity 1 -nx ${XX[j]} -ny ${XX[j]} -nz ${XX[j]} -dx ${DX[p]} -dy ${DY[p]} -dz ${DZ[p]} -prefetcher 0 -increase_factor 1.1 -tests_num 1 -generatematrix 1  -iter_max 1000  -eps 0.001 -attempts ${ATTEMPTS} -check_result 0 -write_matrix 0 -write_performance  1 -performance_filename ${perf_file} -threads ${i} -write_profile 0 -mv_algorithm ${MV_ALGORITHM_ID} -mv_communicator ${MV_COMMUNICATOR_ID} -benchmark_id ${BENCHMARK_ID}
    ;;
    *) 
      aprun -n ${PROCS_ARRAY[p]} ${APRUN_PARM_N} ${BIN}  -verbosity 1 -nx ${XX[j]} -ny ${XX[j]} -nz ${XX[j]} -dx ${DX[p]} -dy ${DY[p]} -dz ${DZ[p]} -prefetcher 0 -increase_factor 1.1 -tests_num 1 -generatematrix 1  -iter_max 1000  -eps 0.001 -attempts ${ATTEMPTS} -check_result 0 -write_matrix 0 -write_performance  1 -performance_filename ${perf_file} -threads ${i} -write_profile 0 -mv_algorithm ${MV_ALGORITHM_ID} -mv_communicator ${MV_COMMUNICATOR_ID} -benchmark_id ${BENCHMARK_ID}
   ;;
   esac 
#   mv *.vtk ${PAT_RT_EXPFILE_DIR}/
  sleep 5
  done
 done
done

echo "#" ${APRUN_PARM} >> ${perf_file}

