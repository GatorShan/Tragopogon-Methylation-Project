#!/bin/bash
#SBATCH --job-name=MethylDiff_CHH_V1
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=shan158538@ufl.edu	
#SBATCH --ntasks=1
#SBATCH --mem=150gb
#SBATCH --time=3-00:00:00
#SBATCH --output=MethylDiff_CHH_V1_%j.out 
#SBATCH --error=MethylDiff_CHH_V1_%j.error
#SBATCH --qos=soltis-b

date; hostname; pwd

module purge
module load R

#Run R script 
Rscript MethylDiff_CHH_V1.r

date
