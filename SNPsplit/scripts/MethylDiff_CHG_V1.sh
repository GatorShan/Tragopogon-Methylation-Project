#!/bin/bash
#SBATCH --job-name=MethylDiff_CHG_V1
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=shan158538@ufl.edu	
#SBATCH --ntasks=1
#SBATCH --mem=30gb
#SBATCH --time=2-00:00:00
#SBATCH --output=MethylDiff_CHG_V1_%j.out 
#SBATCH --error=MethylDiff_CHG_V1_%j.error
#SBATCH --qos=soltis-b

date; hostname; pwd

module purge
module load R

#Run R script 
Rscript MethylDiff_CHG_V1.r

date