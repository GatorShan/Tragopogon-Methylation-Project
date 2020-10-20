#!/bin/bash

#SBATCH --job-name=HaplotypeCaller_V1
#SBATCH --mail-user=shan158538@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --output=HaplotypeCaller_V1_%j.out
#SBATCH --error=HaplotypeCaller_V1_%j.error
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=7gb
#SBATCH --time=04-00:00:00
#SBATCH --qos=soltis

date;hostname;pwd

### This script is used to call SNPs and indels

INPUT=/blue/soltis/shan158538/Methylation/OutPut/SNP_calling/BWA_alignment/Alignment_output/Bam_sorted_marked_files
REF=/blue/soltis/shan158538/Methylation/OutPut/SNP_calling/BWA_alignment
OUT=/blue/soltis/shan158538/Methylation/OutPut/SNP_calling/Variant_discovery

module purge
module load samtools/1.10
module load gatk/4.1.8.1

cd ${REF}

gatk CreateSequenceDictionary \
	-R Tdub.V1.fasta \
	-O Tdub.V1.dict

samtools \
	faidx \
	Tdub.V1.fasta

gatk --java-options "-Xmx5g" HaplotypeCaller \
	-R ${REF}/Tdub.V1.fasta \
	-I ${INPUT}/Tpr_combined_sorted_marked_AddReadGroup.bam \
	-O ${OUT}/Tpr.g.vcf.gz \
	-ERC GVCF

date
