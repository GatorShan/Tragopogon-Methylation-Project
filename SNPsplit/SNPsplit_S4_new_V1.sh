#!/bin/bash

#SBATCH --job-name=SNPsplit_S4_new_V1
#SBATCH --mail-user=shan158538@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --output=SNPsplit_S4_new_V1_%j.out
#SBATCH --error=SNPsplit_S4_new_V1_%j.error
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2gb
#SBATCH --time=7-00:00:00
#SBATCH --qos=soltis

date;hostname;pwd

module purge
module load snpsplit/0.5.0

IN=/blue/soltis/shan158538/Methylation/OutPut/SNPsplit_alignment/S4_dedup
SNP=/blue/soltis/shan158538/Methylation/OutPut/SNP_calling/Variant_filtration/SNP_file_homo_snps_header.txt
OUT=/blue/soltis/shan158538/Methylation/OutPut/SNPsplit_output

SNPsplit \
	--snp_file ${SNP} \
	--paired \
	--bisulfite \
	-o ${OUT} \
	${IN}/S4_cat_R1_val_1_bismark_bt2_pe_nameSorted.deduplicated.bam

date