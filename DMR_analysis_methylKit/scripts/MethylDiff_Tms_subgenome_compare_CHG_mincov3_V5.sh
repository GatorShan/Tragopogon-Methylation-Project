#!/bin/bash
#SBATCH --job-name=MethylDiff_Tms_subgenome_compare_CHG_mincov3_V5
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=shan158538@ufl.edu	
#SBATCH --ntasks=1
#SBATCH --mem=5gb
#SBATCH --time=06:00:00
#SBATCH --output=MethylDiff_Tms_subgenome_compare_CHG_mincov3_V5_%j.out 
#SBATCH --error=MethylDiff_Tms_subgenome_compare_CHG_mincov3_V5_%j.error
#SBATCH --qos=soltis

date; hostname; pwd

module purge
module load R

#Run R script 
Rscript MethylDiff_Tms_subgenome_compare_CHG_mincov3_V5.r

date
