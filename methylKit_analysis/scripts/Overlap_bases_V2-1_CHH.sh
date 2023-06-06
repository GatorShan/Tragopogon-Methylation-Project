#!/bin/bash
#SBATCH --job-name=Overlap_bases_V2-1_CHH
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=shan158538@ufl.edu	
#SBATCH --ntasks=1
#SBATCH --mem=180gb
#SBATCH --time=01:00:00
#SBATCH --output=Overlap_bases_V2-1_CHH_%j.out 
#SBATCH --error=Overlap_bases_V2-1_CHH_%j.error
#SBATCH --qos=soltis-b

date; hostname; pwd

cd /blue/soltis/shan158538/Methylation/OutPut/Overlap_bases/CHH_mincov_3

module purge
module load python

Overlap_bases.V2.py \
    Tdu_1_cov3_CHH.txt \
    Tdu_2_cov3_CHH.txt \
    Tms_1_du_cov3_CHH.txt \
    Tms_1_pr_cov3_CHH.txt \
    Tms_2_du_cov3_CHH.txt \
    Tms_2_pr_cov3_CHH.txt \
    Tpr_1_cov3_CHH.txt \
    Tpr_2_cov3_CHH.txt

date
