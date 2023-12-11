#!/bin/bash

#SBATCH -J simulate_pPGM
#SBATCH -o ./%x.%j.%N.out
#SBATCH -e ./%x.%j.%N.err
#SBATCH -D ./
#SBATCH --get-user-env
#SBATCH --clusters=cm2_tiny
#SBATCH --partition=cm2_tiny
#SBATCH --nodes=1-1
#SBATCH --cpus-per-task=56
# 56 is the maximum reasonable value for CooLMUC-2
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --mail-user=josefine.stark@tum.de
#SBATCH --export=NONE
#SBATCH --time=12:00:00
#SBATCH --begin=now+0hours
#SBATCH --mem=32G 
module load slurm_setup

module load r

Rscript generate_pPGM_cluster.R
