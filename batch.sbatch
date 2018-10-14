#!/bin/bash

#SBATCH --job-name=MPCs
#SBATCH --output=/home/livingstonb/output/matlab.out
#SBATCH --error=/home/livingstonb/output/matlab.err
#SBATCH --partition=broadwl
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=3
#SBATCH --mem-per-cpu=4000

# set variables for use in master.m
setenv Batch 1
setenv Server 1
setenv SmallGrid 0

module load matlab

matlab -nodisplay < master.m
