#!/bin/bash
#SBATCH --job-name=disc
#SBATCH --output=/home/livingstonb/GitHub/Discrete_HA/Output/run%a.out
#SBATCH --error=/home/livingstonb/GitHub/Discrete_HA/Output/run%a.err
#SBATCH --partition=broadwl
#SBATCH --time=12:00:00
#SBATCH --array=1-39
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=10000

module load matlab

matlab -nodisplay < /home/livingstonb/GitHub/Discrete_HA/master.m
