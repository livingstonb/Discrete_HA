#!/bin/bash

#SBATCH --job-name=MPCs
#SBATCH --output==/home/livingstonb/output/matlab_simple.out
#SBATCH --error=/home/livingstonb/output/matlab_simple.err
#SBATCH --partition=broadwl
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --mem-per-cpu=2000

module load matlab

matlab -nodisplay < master.m
