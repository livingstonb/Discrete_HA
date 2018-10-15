#!/bin/bash

#SBATCH --job-name=MPC_single_spec
#SBATCH --output=/home/livingstonb/output/matlab.out
#SBATCH --error=/home/livingstonb/output/matlab.err
#SBATCH --partition=broadwl
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=3
#SBATCH --mem-per-cpu=4000
#SBATCH --export=Batch=0,Server=1,SmallGrid=0

# set variables for use in master.m
export Batch=0
export Server=1
export SmallGrid=0

module load matlab

matlab -nodisplay < master.m
