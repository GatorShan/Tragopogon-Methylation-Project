#!/bin/bash

#SBATCH --job-name=Mask.Tdu.Genome_V1
#SBATCH --mail-user=shan158538@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --output=Mask.Tdu.Genome_V1_%j.out
#SBATCH --error=Mask.Tdu.Genome_V1_%j.error
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=5gb
#SBATCH --time=1-00:00:00
#SBATCH --qos=soltis

date;hostname;pwd

module purge
module load bedtools/2.30.0

IN=/orange/soltis/Tdub_Genome_V1/Tdub.V1.fasta
BED=/blue/soltis/shan158538/Methylation/OutPut/SNP_calling/Variant_filtration/Tpr_combined_filtered.PASS.homo.snps.vcf.gz
OUT=/blue/soltis/shan158538/Methylation/OutPut/SNPsplit/Tdub.V1.masked.fasta

bedtools \
	maskfasta \
	-fi ${IN} \
	-bed ${BED} \
	-fo ${OUT}

date
