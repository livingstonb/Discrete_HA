#!/bin/bash

#SBATCH --job-name=MPCs_Qc
#SBATCH --output=/home/livingstonb/output/matlab_Qc.out
#SBATCH --error=/home/livingstonb/output/matlab_Qc.err
#SBATCH --partition=broadwl
#SBATCH --time=03:30:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=3
#SBATCH --mem-per-cpu=4000

module load matlab

matlab -nodisplay < master.m
