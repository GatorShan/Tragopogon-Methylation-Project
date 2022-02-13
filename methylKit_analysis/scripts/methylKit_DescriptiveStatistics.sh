#!/bin/bash
#SBATCH --job-name=methylKit_DescriptiveStatistics
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=shan158538@ufl.edu	
#SBATCH --ntasks=1
#SBATCH --mem=60gb
#SBATCH --time=2-00:00:00
#SBATCH --output=methylKit_DescriptiveStatistics_%j.out 
#SBATCH --error=methylKit_DescriptiveStatistics_%j.error

date; hostname; pwd

module purge
module load R

#Run R script 
Rscript methylKit_DescriptiveStatistics.r

date
