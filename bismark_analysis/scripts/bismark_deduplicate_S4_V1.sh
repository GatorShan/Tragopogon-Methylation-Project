#!/bin/bash

#SBATCH --job-name=bismark_deduplicate_S4_V1
#SBATCH --mail-user=shan158538@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --output=bismark_deduplicate_S4_V1_%j.out
#SBATCH --error=bismark_deduplicate_S4_V1_%j.error
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=40gb
#SBATCH --time=1-00:00:00
#SBATCH --qos=soltis-b

date;hostname;pwd

IN=/orange/soltis/shan158538/Methylation_output/bismark_alignment/S4
OUT=/orange/soltis/shan158538/Methylation_output/bismark_deduplicate

module purge
module load bismark/0.22.3

deduplicate_bismark \
	--output_dir ${OUT} \
	--paired \
	${IN}/S4_cat_R1_val_1_bismark_bt2_pe_nameSorted.bam

date
