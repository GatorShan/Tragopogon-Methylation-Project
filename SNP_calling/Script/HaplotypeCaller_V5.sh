#!/bin/bash

#SBATCH --job-name=HaplotypeCaller_V5
#SBATCH --mail-user=shan158538@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --output=HaplotypeCaller_V5_%j.out
#SBATCH --error=HaplotypeCaller_V5_%j.error
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=7gb
#SBATCH --time=4-00:00:00
#SBATCH --qos=soltis-b
#SBATCH --array=1-20%8

date;hostname;pwd

### Array job

INPUT=/blue/soltis/shan158538/Methylation/OutPut/SNP_calling/BWA_alignment/Alignment_output/Bam_sorted_marked_files
REF=/blue/soltis/shan158538/Methylation/OutPut/SNP_calling/BWA_alignment
OUT=/blue/soltis/shan158538/Methylation/OutPut/SNP_calling/Variant_discovery

module purge
module load gatk/4.1.8.1

gatk --java-options "-Xmx5g" \
	HaplotypeCaller \
	-R ${REF}/Tdub.V1.fasta \
	-I ${INPUT}/Tpr_combined_sorted_marked_AddReadGroup.bam \
	-O ${OUT}/Tpr_set_${SLURM_ARRAY_TASK_ID}.g.vcf.gz \
	-L ${OUT}/Scaffold_set_${SLURM_ARRAY_TASK_ID}.list \
	--native-pair-hmm-threads 4 \
	-ERC GVCF

date
