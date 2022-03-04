#!/bin/bash
#SBATCH --job-name=processBismarkAln_subgenome_compare_CHH_V1
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=shan158538@ufl.edu	
#SBATCH --ntasks=1
#SBATCH --mem=90gb
#SBATCH --time=3-00:00:00
#SBATCH --output=processBismarkAln_subgenome_compare_CHH_V1_%j.out 
#SBATCH --error=processBismarkAln_subgenome_compare_CHH_V1_%j.error
#SBATCH --qos=soltis-b

date; hostname; pwd

module purge
module load R

#Run R script 
Rscript processBismarkAln_subgenome_compare_CHH_V1.r

date