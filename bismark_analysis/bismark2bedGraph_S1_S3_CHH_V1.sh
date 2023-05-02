#!/bin/bash

#SBATCH --job-name=bismark2bedGraph_S1_S3_CHH_V1
#SBATCH --mail-user=shan158538@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --output=bismark2bedGraph_S1_S3_CHH_V1_%j.out
#SBATCH --error=bismark2bedGraph_S1_S3_CHH_V1_%j.error
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8gb
#SBATCH --time=5-00:00:00
#SBATCH --qos=soltis

date;hostname;pwd

IN=/orange/soltis/shan158538/Methylation_output/bismark_methylation_extractor

module purge
module load bismark/0.22.3

for SAMPLE in S1 S3; do
	for CONTEXT in CHH; do
		echo "Extracting methylation status of ${CONTEXT} context from sample ${SAMPLE}"
		bismark2bedGraph \
			--output ${SAMPLE}_${CONTEXT} \
			--dir ${IN}/${SAMPLE} \
			--CX_context \
			--buffer_size 90% \
			${IN}/${SAMPLE}/${CONTEXT}_OB_${SAMPLE}_cat_R1_val_1_bismark_bt2_pe_nameSorted.deduplicated.txt.gz \
			${IN}/${SAMPLE}/${CONTEXT}_OT_${SAMPLE}_cat_R1_val_1_bismark_bt2_pe_nameSorted.deduplicated.txt.gz
	done
done

date
