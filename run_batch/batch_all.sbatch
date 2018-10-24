#!/bin/bash

#SBATCH --job-name=MPCs_larger_grid_for_beta_het
#SBATCH --output=/home/livingstonb/output/matlab.out
#SBATCH --error=/home/livingstonb/output/matlab.err
#SBATCH --partition=broadwl
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=2
#SBATCH --mem-per-cpu=6000

module load matlab

matlab -nodisplay < /home/livingstonb/GitHub/MPCrecode/master.m
