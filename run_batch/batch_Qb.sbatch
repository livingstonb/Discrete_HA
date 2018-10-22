#!/bin/bash

#SBATCH --job-name=MPCs_Qb
#SBATCH --output=/home/livingstonb/output/matlab_Qb.out
#SBATCH --error=/home/livingstonb/output/matlab_Qb.err
#SBATCH --partition=broadwl
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=3
#SBATCH --mem-per-cpu=4000

module load matlab

matlab -nodisplay < /home/livingstonb/GitHub/MPCrecode/run_batch/temp_master/master_Qb.m
