#!/bin/bash

#SBATCH --job-name=bismark_deduplicate_DES1_S1_S3_V1
#SBATCH --mail-user=shan158538@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --output=bismark_deduplicate_DES1_S1_S3_V1_%j.out
#SBATCH --error=bismark_deduplicate_DES1_S1_S3_V1_%j.error
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=10gb
#SBATCH --time=4-00:00:00
#SBATCH --qos=soltis-b

date;hostname;pwd

IN1=/orange/soltis/shan158538/Methylation_output/bismark_alignment/DES1
IN2=/orange/soltis/shan158538/Methylation_output/bismark_alignment/S1
IN3=/orange/soltis/shan158538/Methylation_output/bismark_alignment/S3
OUT=/orange/soltis/shan158538/Methylation_output/bismark_deduplicate

module purge
module load samtools/1.9
module load bismark/0.22.3

### Processing DES1
samtools sort \
	${IN1}/HMCWKCCXY_s8_1_4981-LF_17_SL334590_val_1_bismark_bt2_pe.bam \
	-o ${IN1}/HMCWKCCXY_s8_1_4981-LF_17_SL334590_val_1_bismark_bt2_pe_nameSorted.bam \
	-n

deduplicate_bismark \
	--output_dir ${OUT} \
	--paired \
	${IN1}/HMCWKCCXY_s8_1_4981-LF_17_SL334590_val_1_bismark_bt2_pe_nameSorted.bam

### Processing S1
samtools sort \
	${IN2}/S1_cat_R1_val_1_bismark_bt2_pe.bam \
	-o ${IN2}/S1_cat_R1_val_1_bismark_bt2_pe_nameSorted.bam \
	-n

deduplicate_bismark \
	--output_dir ${OUT} \
	--paired \
	${IN2}/S1_cat_R1_val_1_bismark_bt2_pe_nameSorted.bam

### Processing S3
samtools sort \
	${IN3}/S3_cat_R1_val_1_bismark_bt2_pe.bam \
	-o ${IN3}/S3_cat_R1_val_1_bismark_bt2_pe_nameSorted.bam \
	-n

deduplicate_bismark \
	--output_dir ${OUT} \
	--paired \
	${IN3}/S3_cat_R1_val_1_bismark_bt2_pe_nameSorted.bam

date
