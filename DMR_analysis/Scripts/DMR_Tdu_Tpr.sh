#!/bin/bash
#SBATCH --job-name=DMR_Tdu_Tpr
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=shan158538@ufl.edu	
#SBATCH --ntasks=1
#SBATCH --mem=60gb
#SBATCH --time=6-00:00:00
#SBATCH --output=DMR_Tdu_Tpr.%j.out 

date; hostname; pwd

module purge
module load R

#Run R script 
Rscript DMR_Tdu_Tpr.r

date