#!/bin/bash
#SBATCH -J HTMTsim50800          
#SBATCH -c 32                   
#SBATCH --mem=8G                
#SBATCH -p small_cpu
#SBATCH --tmp=5G                             
#SBATCH --mail-type=ALL          
#SBATCH --mail-user=jason.berger@uni-wuerzburg.de

cd ~/confruns/
srun R --save < HTMTsim50800.R