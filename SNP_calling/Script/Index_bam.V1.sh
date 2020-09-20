#!/bin/bash

#SBATCH --job-name=Index_bam.V1
#SBATCH --mail-user=shan158538@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --output=Index_bam.V1_%j.out
#SBATCH --error=Index_bam.V1_%j.error
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=3
#SBATCH --mem=10gb
#SBATCH --time=01-00:00:00
#SBATCH --qos=soltis

date;hostname;pwd

### This script is used to mark the duplicates of the BAM files

INPUT=/blue/soltis/shan158538/Methylation/OutPut/SNP_calling/BWA_alignment/Alignment_output/Bam_sorted_marked_files

module purge
module load picard/2.21.2

cd ${INPUT}

for sample in 533 534 535 536 537 538 539 540 541 542 543 544 545 546 547 548
do

### INDEX THE BAM FILE, THE OUTPUT FILES COULD BE VISUALIZED BY USING IGV
echo "Indexing sample ${sample}"
java -Xmx8g -jar $HPC_PICARD_DIR/picard.jar \
	BuildBamIndex \
	I=Tpr_${sample}_sorted_marked.bam
echo "Sample ${sample} has been indexed"

done

date
