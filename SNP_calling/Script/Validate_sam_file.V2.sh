#!/bin/bash

#SBATCH --job-name=Validate_sam_file.V2
#SBATCH --mail-user=shan158538@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --output=Validate_sam_file.V2_%j.out
#SBATCH --error=Validate_sam_file.V2_%j.error
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=7gb
#SBATCH --time=01-00:00:00
#SBATCH --qos=soltis

date;hostname;pwd

### This script is used to generate detailed list of ERROR records

INPUT=/blue/soltis/shan158538/Methylation/OutPut/SNP_calling/BWA_alignment/Alignment_output/Bam_sorted_marked_files

module purge
module load picard/2.21.2

java -Xmx5g -jar $HPC_PICARD_DIR/picard.jar \
	ValidateSamFile \
	I=${INPUT}/Tpr_combined_sorted_marked.bam \
	IGNORE_WARNINGS=true \
	O=${INPUT}/Tpr_combined_sorted_marked_bam_validation_ErrorDetail.txt \
	MODE=VERBOSE

date
