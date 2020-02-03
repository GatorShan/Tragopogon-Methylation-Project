#!/bin/bash

#SBATCH --job-name=bismark_methylation_extractor_DES1_V1
#SBATCH --mail-user=shan158538@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --output=bismark_methylation_extractor_DES1_V1_%j.out
#SBATCH --error=bismark_methylation_extractor_DES1_V1_%j.error
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --mem=10gb
#SBATCH --time=4-00:00:00
#SBATCH --qos=soltis-b

date;hostname;pwd

IN=/orange/soltis/shan158538/Methylation_output/bismark_deduplicate
OUT=/orange/soltis/shan158538/Methylation_output/bismark_methylation_extractor/DES1

module purge
module load bismark/0.22.3

bismark_methylation_extractor \
	${IN}/HMCWKCCXY_s8_1_4981-LF_17_SL334590_val_1_bismark_bt2_pe_nameSorted.deduplicated.bam \
	--paired-end \
	--no_overlap \
	--output ${OUT} \
	--gzip \
	--multicore 3 \
	--bedGraph \
	--scaffolds

date
