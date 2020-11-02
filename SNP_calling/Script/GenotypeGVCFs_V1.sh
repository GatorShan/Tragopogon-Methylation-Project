#!/bin/bash

#SBATCH --job-name=GenotypeGVCFs_V1
#SBATCH --mail-user=shan158538@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --output=GenotypeGVCFs_V1_%j.out
#SBATCH --error=GenotypeGVCFs_V1_%j.error
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=7gb
#SBATCH --time=4-00:00:00
#SBATCH --qos=soltis

date;hostname;pwd

### Perform joint genotyping on one sample

REF=/blue/soltis/shan158538/Methylation/OutPut/SNP_calling/BWA_alignment
INPUT=/blue/soltis/shan158538/Methylation/OutPut/SNP_calling/Variant_discovery

module purge
module load gatk/4.1.8.1

gatk --java-options "-Xmx5g" \
	GenotypeGVCFs \
	-R ${REF}/Tdub.V1.fasta \
	-V ${INPUT}/Tpr_combined.g.vcf.gz \
	-O ${INPUT}/Tpr_combined.vcf.gz

date
