#!/bin/bash

#SBATCH --job-name=Reformat_CHG_cov_files_shared_V1
#SBATCH --mail-user=shan158538@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --output=Reformat_CHG_cov_files_shared_V1_%j.out
#SBATCH --error=Reformat_CHG_cov_files_shared_V1_%j.error
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=10gb
#SBATCH --time=1-00:00:00
#SBATCH --qos=soltis-b

date;hostname;pwd

module purge
module load python

IN=/orange/soltis/shan158538/Methylation_output/bismark_coverage_files
CONTEXT=CHG
OUT=/orange/soltis/shan158538/Methylation_output/feature_methylation/gene_${CONTEXT}_shared

cd ${IN}

for SAMPLE in DES1 S1 S2 S3 S4 S5; do
	echo "Reformatting coverage file from sample ${SAMPLE} with ${CONTEXT} context"
	Reformat_bismark_cov_V3.py \
		${SAMPLE}_${CONTEXT}.gz.bismark_shared_filtered.cov \
		-o=${OUT}/${SAMPLE}/
done

date
