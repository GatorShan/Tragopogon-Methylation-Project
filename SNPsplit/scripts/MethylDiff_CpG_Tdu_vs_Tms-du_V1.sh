#!/bin/bash
#SBATCH --job-name=MethylDiff_CpG_Tdu_vs_Tms-du_V1
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=shan158538@ufl.edu	
#SBATCH --ntasks=1
#SBATCH --mem=30gb
#SBATCH --time=12:00:00
#SBATCH --output=MethylDiff_CpG_Tdu_vs_Tms-du_V1_%j.out 
#SBATCH --error=MethylDiff_CpG_Tdu_vs_Tms-du_V1_%j.error

date; hostname; pwd

module purge
module load R

#Run R script 
Rscript MethylDiff_CpG_Tdu_vs_Tms-du_V1.r

date
