#!/bin/bash
# run all specifications

cd ~/GitHub/MPCrecode/batch/sbatch
array=$(ls *.sbatch)

for i in ${array[@]}
do
  sbatch "$i"
done
