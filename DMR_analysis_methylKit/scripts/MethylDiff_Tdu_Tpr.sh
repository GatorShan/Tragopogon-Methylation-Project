#!/bin/bash
#SBATCH --job-name=MethylDiff_Tdu_Tpr
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=shan158538@ufl.edu	
#SBATCH --ntasks=1
#SBATCH --mem=60gb
#SBATCH --time=2-00:00:00
#SBATCH --output=MethylDiff_Tdu_Tpr_%j.out 
#SBATCH --error=MethylDiff_Tdu_Tpr_%j.error

date; hostname; pwd

module purge
module load R

#Run R script 
Rscript MethylDiff_Tdu_Tpr.r

date
