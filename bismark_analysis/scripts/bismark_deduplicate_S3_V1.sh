#!/bin/bash

#SBATCH --job-name=bismark_deduplicate_S3_V1
#SBATCH --mail-user=shan158538@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --output=bismark_deduplicate_S3_V1_%j.out
#SBATCH --error=bismark_deduplicate_S3_V1_%j.error
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=20gb
#SBATCH --time=10:00:00
#SBATCH --qos=soltis-b

date;hostname;pwd

IN1=/orange/soltis/shan158538/Methylation_output/bismark_alignment/DES1
IN2=/orange/soltis/shan158538/Methylation_output/bismark_alignment/S1
IN3=/orange/soltis/shan158538/Methylation_output/bismark_alignment/S3
OUT=/orange/soltis/shan158538/Methylation_output/bismark_deduplicate

module purge
##module load samtools/1.9
module load bismark/0.22.3

### Processing S3

deduplicate_bismark \
	--output_dir ${OUT} \
	--paired \
	${IN3}/S3_cat_R1_val_1_bismark_bt2_pe_nameSorted.bam

date
