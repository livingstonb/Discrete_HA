#!/bin/bash

#SBATCH --job-name=MPCs
#SBATCH --output=/home/livingstonb/output/matlab.out
#SBATCH --error=/home/livingstonb/output/matlab.err
#SBATCH --partition=broadwl
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=3
#SBATCH --mem-per-cpu=4000
#SBATCH --export=Batch=1,Server=1,SmallGrid=0

module load matlab

matlab -nodisplay < master.m
