#!/bin/bash
#SBATCH --job-name=MethylDiff_Tdu-Tpr_CpG_mincov3_V4
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=shan158538@ufl.edu	
#SBATCH --ntasks=1
#SBATCH --mem=5gb
#SBATCH --time=06:00:00
#SBATCH --output=MethylDiff_Tdu-Tpr_CpG_mincov3_V4_%j.out 
#SBATCH --error=MethylDiff_Tdu-Tpr_CpG_mincov3_V4_%j.error
#SBATCH --qos=soltis

date; hostname; pwd

module purge
module load R

#Run R script 
Rscript MethylDiff_Tdu-Tpr_CpG_mincov3_V4.r

date
