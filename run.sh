#!/bin/bash

# Limit ROOT and BLAS threading
#export ROOT_TTHREADS=1
#export OMP_NUM_THREADS=1
#export MKL_NUM_THREADS=1
#export OPENBLAS_NUM_THREADS=1
#export VECLIB_MAXIMUM_THREADS=1
#export NUMEXPR_NUM_THREADS=1

q_ranges=("q10to100" "q1to10")

for qrange in "${q_ranges[@]}"; do
  mkdir -p "./log/${qrange}"   # create log dir for this qrange
  mkdir -p "./results/${qrange}"
  chmod 777 "log"
  chmod 777 "results"
  for i in $(seq -w 0 99); do
    input_file="inputfiles/ep_10x100_${qrange}/split/part_${i}.list"
    res_dir="./results/${qrange}"
    output_file="${res_dir}/JetOut_${qrange}_${i}.histos.root"
    log_dir="./log/${qrange}"
    log_file="${log_dir}/out_${qrange}_${i}.log"
    
    echo "Running ROOT job for $qrange part $i ..."
    root -l -q "src/JetAnalyzer.C(\"$input_file\",\"$output_file\")" &> "$log_file"
  done
done

