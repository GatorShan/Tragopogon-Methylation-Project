#!/bin/bash

#SBATCH --job-name=Fastqc_quality_check.V2
#SBATCH --mail-user=shan158538@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --output=Fastqc_quality_check.V2_%j.out
#SBATCH --error=Fastqc_quality_check.V2_%j.error
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=1gb
#SBATCH --time=3-00:00:00
#SBATCH --qos=soltis

date;hostname;pwd

IN=/orange/soltis/shan158538/Methylation_output/trim_galore
OUT=/orange/soltis/shan158538/Methylation_output/FastQC_result

module purge
module load fastqc/0.11.7

cd ${IN}

for sample in S1 S2 S3; do
	fastqc \
		-o ${OUT} \
		-t 4 \
		${sample}_cat_R1_val_1.fq.gz \
		${sample}_cat_R2_val_2.fq.gz
done

date
