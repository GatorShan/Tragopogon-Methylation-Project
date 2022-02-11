#!/bin/bash
#SBATCH --job-name=processBismarkAln_CG_V1
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=shan158538@ufl.edu	
#SBATCH --ntasks=1
#SBATCH --mem=60gb
#SBATCH --time=6-00:00:00
#SBATCH --output=processBismarkAln_CG_V1_%j.out 
#SBATCH --error=processBismarkAln_CG_V1_%j.error

date; hostname; pwd

module purge
module load R

#Run R script 
Rscript processBismarkAln_CG_V1.r

date
