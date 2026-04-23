#!/bin/bash
#SBATCH -J ArrayJob_HTMT
#SBATCH -p small_cpu
#SBATCH --array=1-96
#SBATCH --qos=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=6G
#SBATCH --output=outfiles/%x_%A_%a.out
#SBATCH --error=errorfiles/%x_%A_%a.err
#SBATCH --time=24:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jason.berger@uni-wuerzburg.de


echo "Running Task: ${SLURM_ARRAY_TASK_ID} on host $(hostname)"
module load R
srun Rscript code/rerunsim.R