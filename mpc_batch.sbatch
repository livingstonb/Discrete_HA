#!/bin/bash

#SBATCH --job-name=MPCs
#SBATCH --output=/home/livingstonb/output/matlab.out
#SBATCH --error=/home/livingstonb/output/matlab.err
#SBATCH --partition=broadwl
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=6
#SBATCH --mem-per-cpu=2000

module load matlab

matlab -nodisplay < master.m
