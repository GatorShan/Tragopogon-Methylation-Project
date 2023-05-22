#!/bin/bash
#SBATCH --job-name=overlap_min-cov_CHG
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=shan158538@ufl.edu	
#SBATCH --ntasks=1
#SBATCH --mem=80gb
#SBATCH --time=00:30:00
#SBATCH --output=overlap_min-cov_CHG_%j.out 
#SBATCH --error=overlap_min-cov_CHG_%j.error
#SBATCH --qos=soltis-b

date; hostname; pwd

cd /orange/soltis/shan158538/Methylation_output/bismark_coverage_files

module purge
module load python3

overlap_min-cov_4.py DES1_CHG.gz.bismark.cov S1_CHG.gz.bismark.cov S2_CHG.gz.bismark.cov S3_CHG.gz.bismark.cov S4_CHG.gz.bismark.cov S5_CHG.gz.bismark.cov

date
