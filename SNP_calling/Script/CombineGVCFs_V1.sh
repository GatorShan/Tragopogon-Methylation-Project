#!/bin/bash

#SBATCH --job-name=CombineGVCFs_V1
#SBATCH --mail-user=shan158538@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --output=CombineGVCFs_V1_%j.out
#SBATCH --error=CombineGVCFs_V1_%j.error
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=7gb
#SBATCH --time=4-00:00:00
#SBATCH --qos=soltis-b

date;hostname;pwd

### Merge many HaplotypeCaller GVCF files into a single GVCF

REF=/blue/soltis/shan158538/Methylation/OutPut/SNP_calling/BWA_alignment
INPUT=/blue/soltis/shan158538/Methylation/OutPut/SNP_calling/Variant_discovery

module purge
module load gatk/4.1.8.1

gatk --java-options "-Xmx5g" \
	CombineGVCFs \
	-R ${REF}/Tdub.V1.fasta \
	--variant ${INPUT}/Tpr_set_1.g.vcf.gz \
	--variant ${INPUT}/Tpr_set_2.g.vcf.gz \
	--variant ${INPUT}/Tpr_set_3.g.vcf.gz \
	--variant ${INPUT}/Tpr_set_4.g.vcf.gz \
	--variant ${INPUT}/Tpr_set_5.g.vcf.gz \
	--variant ${INPUT}/Tpr_set_6.g.vcf.gz \
	--variant ${INPUT}/Tpr_set_7.g.vcf.gz \
	--variant ${INPUT}/Tpr_set_8.g.vcf.gz \
	--variant ${INPUT}/Tpr_set_9.g.vcf.gz \
	--variant ${INPUT}/Tpr_set_10.g.vcf.gz \
	--variant ${INPUT}/Tpr_set_11.g.vcf.gz \
	--variant ${INPUT}/Tpr_set_12.g.vcf.gz \
	--variant ${INPUT}/Tpr_set_13.g.vcf.gz \
	--variant ${INPUT}/Tpr_set_14.g.vcf.gz \
	--variant ${INPUT}/Tpr_set_15.g.vcf.gz \
	--variant ${INPUT}/Tpr_set_16.g.vcf.gz \
	--variant ${INPUT}/Tpr_set_17.g.vcf.gz \
	--variant ${INPUT}/Tpr_set_18.g.vcf.gz \
	--variant ${INPUT}/Tpr_set_19.g.vcf.gz \
	--variant ${INPUT}/Tpr_set_20.g.vcf.gz \
	-O ${INPUT}/Tpr_combined.g.vcf.gz

date
