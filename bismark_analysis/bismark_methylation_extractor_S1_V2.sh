#!/bin/bash

#SBATCH --job-name=bismark_methylation_extractor_S1_V2
#SBATCH --mail-user=shan158538@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --output=bismark_methylation_extractor_S1_V2_%j.out
#SBATCH --error=bismark_methylation_extractor_S1_V2_%j.error
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=10gb
#SBATCH --time=4-00:00:00
#SBATCH --qos=soltis

date;hostname;pwd

IN=/orange/soltis/shan158538/Methylation_output/bismark_deduplicate
OUT=/orange/soltis/shan158538/Methylation_output/bismark_methylation_extractor/S1

module purge
module load bismark/0.22.3

bismark_methylation_extractor \
	${IN}/S1_cat_R1_val_1_bismark_bt2_pe_nameSorted.deduplicated.bam \
	--paired-end \
	--no_overlap \
	--output ${OUT} \
	--gzip \
	--multicore 4 \
	--bedGraph \
	--scaffolds \
	--ignore_r2 2

date
