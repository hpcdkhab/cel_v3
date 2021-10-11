#!/bin/sh

perf_file="poisson_3d_64_x_64_star_27_procs_per_node_2_pinned.dat"

rm ${perf_file} -rf

PROCS_ARRAY=("not used" "1" "2" "4" "8" "16")
DX=("not used" "1" "2" "2" "2" "4"  "4")
DY=("not used" "1" "1" "2" "2" "2"  "4")    
DZ=("not used" "1" "1" "1" "2" "2"  "2")

#export I_MPI_PIN_PROCESSOR_LIST=p0,p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,p16
for p in {1..2}
do 
 for i in {2..8}
 do
  for j in {1..5} 
  do
   echo "procs:" ${PROCS_ARRAY[p]}  "threads per proc:" $i "test id:" ${j}
   echo "DX:" ${DX[p]} "DY:"  ${DY[p]} "DZ:" ${DZ[p]}
   mpirun -n ${p} ./tests/cel_cg/cel_cg   -verbosity 1 -nx 16 -ny 16 -nz 16 -dx ${DX[p]} -dy ${DY[p]} -dz ${DZ[p]} -prefetcher 0 -increase_factor 1.1 -tests_num 1 -generatematrix 1  -iter_max 1000  -eps 0.001 -attempts 100 -check_result 0 -write_matrix 0 -write_performance  1 -performance_filename ${perf_file} -threads ${i}
   sleep 2
  done
 done
done


