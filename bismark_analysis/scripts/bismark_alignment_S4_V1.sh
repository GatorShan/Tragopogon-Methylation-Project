#!/bin/bash

#SBATCH --job-name=bismark_alignment_S4_V1
#SBATCH --mail-user=shan158538@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --output=bismark_alignment_S4_V1_%j.out
#SBATCH --error=bismark_alignment_S4_V1_%j.error
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=20gb
#SBATCH --time=10-00:00:00
#SBATCH --qos=soltis

date;hostname;pwd

REF=/orange/soltis/shan158538/Methylation_output/bismark_genome_prep
IN=/orange/soltis/shan158538/Methylation_output/trim_galore
OUT=/orange/soltis/shan158538/Methylation_output/bismark_alignment/S4
TEMP=/orange/soltis/shan158538/methylpy_temp

module purge
module load bismark/0.22.3

bismark \
	--genome ${REF} \
	-1 ${IN}/S4_cat_R1_val_1.fq.gz \
	-2 ${IN}/S4_cat_R2_val_2.fq.gz \
	--output_dir ${OUT} \
	--multicore 5 \
	--temp_dir ${TEMP}

date
