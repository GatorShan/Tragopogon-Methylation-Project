#!/bin/bash

#SBATCH --job-name=VariantFiltration_V1
#SBATCH --mail-user=shan158538@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --output=VariantFiltration_V1_%j.out
#SBATCH --error=VariantFiltration_V1_%j.error
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=7gb
#SBATCH --time=4-00:00:00
#SBATCH --qos=soltis

date;hostname;pwd

### Filter variant

INPUT=/blue/soltis/shan158538/Methylation/OutPut/SNP_calling/Variant_discovery
OUTPUT=/blue/soltis/shan158538/Methylation/OutPut/SNP_calling/Variant_filtration

module purge
module load gatk/4.1.8.1

gatk --java-options "-Xmx5g" \
	VariantFiltration \
	-V ${INPUT}/Tpr_combined.vcf.gz \
	-filter "QD < 2.0" --filter-name "QD2" \
	-filter "QUAL < 30.0" --filter-name "QUAL30" \
	-filter "SOR > 3.0" --filter-name "SOR3" \
	-filter "FS > 60.0" --filter-name "FS60" \
	-filter "MQ < 40.0" --filter-name "MQ40" \
	-filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
	-filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
	-O ${OUTPUT}/Tpr_combined_filtered.vcf.gz

date
