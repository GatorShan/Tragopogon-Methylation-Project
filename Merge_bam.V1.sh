#!/bin/bash

#SBATCH --job-name=Merge_bam.V1
#SBATCH --mail-user=shan158538@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --output=Merge_bam.V1_%j.out
#SBATCH --error=Merge_bam.V1_%j.error
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=7gb
#SBATCH --time=01-00:00:00
#SBATCH --qos=soltis

date;hostname;pwd

### This script is used to combine the sorted and markdup BAM files into a single file

INPUT=/blue/soltis/shan158538/Methylation/OutPut/SNP_calling/BWA_alignment/Alignment_output/Bam_sorted_marked_files

module purge
module load picard/2.21.2

cd ${INPUT}

### COMBINE BAM FILES
java -Xmx5g -jar $HPC_PICARD_DIR/picard.jar \
	MergeSamFiles \
	I=Tpr_533_sorted_marked.bam \
	I=Tpr_534_sorted_marked.bam \
	I=Tpr_535_sorted_marked.bam \
	I=Tpr_536_sorted_marked.bam \
	I=Tpr_537_sorted_marked.bam \
	I=Tpr_538_sorted_marked.bam \
	I=Tpr_539_sorted_marked.bam \
	I=Tpr_540_sorted_marked.bam \
	I=Tpr_541_sorted_marked.bam \
	I=Tpr_542_sorted_marked.bam \
	I=Tpr_543_sorted_marked.bam \
	I=Tpr_544_sorted_marked.bam \
	I=Tpr_545_sorted_marked.bam \
	I=Tpr_546_sorted_marked.bam \
	I=Tpr_547_sorted_marked.bam \
	I=Tpr_548_sorted_marked.bam \
	O=Tpr_combined_sorted_marked.bam

### INDEX THE BAM FILE
java -Xmx5g -jar $HPC_PICARD_DIR/picard.jar \
	BuildBamIndex \
	I=Tpr_combined_sorted_marked.bam

date
