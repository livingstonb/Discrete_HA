#!/bin/bash
# run all specifications
declare -a batches=("batch_all" "batch_Qa" "batch_Qb" "batch_Qc")

for i in ${batches[@]}
do
  sbatch "${i}.sbatch"
done
