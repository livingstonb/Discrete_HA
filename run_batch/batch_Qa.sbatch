#!/bin/bash

#SBATCH --job-name=MPCs_Qa
#SBATCH --output=/home/livingstonb/output/matlab_Qa.out
#SBATCH --error=/home/livingstonb/output/matlab_Qa.err
#SBATCH --partition=broadwl
#SBATCH --time=03:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=3
#SBATCH --mem-per-cpu=4000

module load matlab

matlab -nodisplay < /home/livingstonb/GitHub/MPCrecode/run_batch/temp_master/master_Qa.m
