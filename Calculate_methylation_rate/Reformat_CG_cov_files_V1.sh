#!/bin/bash

#SBATCH --job-name=Reformat_CG_cov_files_V1
#SBATCH --mail-user=shan158538@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --output=Reformat_CG_cov_files_V1_%j.out
#SBATCH --error=Reformat_CG_cov_files_V1_%j.error
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=10gb
#SBATCH --time=2-00:00:00
#SBATCH --qos=soltis

date;hostname;pwd

module purge
module load python

IN=/orange/soltis/shan158538/Methylation_output/bismark_coverage_files
CONTEXT=CG
OUT=/orange/soltis/shan158538/Methylation_output/feature_methylation/gene_CG

cd ${IN}

for SAMPLE in DES1 S1 S2 S3 S4 S5; do
	echo "Reformatting coverage file from sample ${SAMPLE} with ${CONTEXT} context"
	Reformat_bismark_cov_V3.py \
		${SAMPLE}_CpG.gz.bismark.cov \
		-o=${OUT}/${SAMPLE}/
done

date
