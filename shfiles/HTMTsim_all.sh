#!/bin/bash
#SBATCH -J HTMTsim_all
#SBATCH -c 32
#SBATCH --mem=32G
#SBATCH -p small_cpu
#SBATCH --tmp=32G
#SBATCH --output=HTMTsim_all_%j.out
#SBATCH --error=HTMTsim_all_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jason.berger@uni-wuerzburg.de

mkdir -p results
srun Rscript HTMTsim_all.R
