#!/bin/bash

#SBATCH --job-name=bismark_deduplicate_S4_V2
#SBATCH --mail-user=shan158538@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --output=bismark_deduplicate_S4_V2_%j.out
#SBATCH --error=bismark_deduplicate_S4_V2_%j.error
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=40gb
#SBATCH --time=1-00:00:00
#SBATCH --qos=soltis-b

date;hostname;pwd

IN=/blue/soltis/shan158538/Methylation/OutPut/SNPsplit_alignment/S4
OUT=/blue/soltis/shan158538/Methylation/OutPut/SNPsplit_alignment/S4_dedup

module purge
module load samtools/1.9
module load bismark/0.22.3

samtools sort \
	${IN}/S4_cat_R1_val_1_bismark_bt2_pe.bam \
	-o ${IN}/S4_cat_R1_val_1_bismark_bt2_pe_nameSorted.bam \
	-n

deduplicate_bismark \
	--output_dir ${OUT} \
	--paired \
	${IN}/S4_cat_R1_val_1_bismark_bt2_pe_nameSorted.bam

date
