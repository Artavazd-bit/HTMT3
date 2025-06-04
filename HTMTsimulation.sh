#!/bin/bash
#SBATCH -J HTMTsimulation          
#SBATCH -c 32                   
#SBATCH --mem=8G                
#SBATCH -p small_cpu
#SBATCH --tmp=5G                             
#SBATCH --mail-type=ALL          
#SBATCH --mail-user=jason.berger@uni-wuerzburg.de

cd ~/R-projects/HTMT/newrun
srun R --save < HTMTsimulation.R