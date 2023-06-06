#!/bin/bash
#SBATCH --job-name=processBismarkAln_CHH_minCOV-3_V1
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=shan158538@ufl.edu	
#SBATCH --ntasks=1
#SBATCH --mem=120gb
#SBATCH --time=1-00:00:00
#SBATCH --output=processBismarkAln_CHH_minCOV-3_V1_%j.out 
#SBATCH --error=processBismarkAln_CHH_minCOV-3_V1_%j.error
#SBATCH --qos=soltis-b

date; hostname; pwd

module purge
module load R

#Run R script
Rscript processBismarkAln_CHH_minCOV-3_V1.r

date
