#!/bin/bash
#SBATCH --job-name=disc_nb1
#SBATCH --output=/home/livingstonb/GitHub/MPCrecode/Output/run%a.out
#SBATCH --error=/home/livingstonb/GitHub/MPCrecode/Output/run%a.err
#SBATCH --partition=broadwl
#SBATCH --time=01:00:00
#SBATCH --array=1-64
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=4000

module load matlab

matlab -nodisplay < /home/livingstonb/GitHub/MPCrecode/master.m