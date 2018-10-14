#!/bin/bash

#SBATCH --job-name=MPC_single_spec
#SBATCH --output=/home/livingstonb/output/matlab.out
#SBATCH --error=/home/livingstonb/output/matlab.err
#SBATCH --partition=broadwl
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=3
#SBATCH --mem-per-cpu=4000

# set variables for use in master.m
setenv Batch 0
setenv Server 1
setenv SmallGrid 0

module load matlab

matlab -nodisplay < master.m
