#!/bin/bash

#SBATCH -J H0_NCT
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

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/lrz/sys/spack/release/22.2.1/opt/x86_64/libjpeg-turbo/2.1.0-gcc-urdhzdt/lib64/

module load r

Rscript test_H0_NCT.R
