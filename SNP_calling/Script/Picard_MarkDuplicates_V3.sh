#!/bin/bash

#SBATCH --job-name=Picard_MarkDuplicates_V3
#SBATCH --mail-user=shan158538@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --output=Picard_MarkDuplicates_V3_%j.out
#SBATCH --error=Picard_MarkDuplicates_V3_%j.error
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=25gb
#SBATCH --time=04-00:00:00
#SBATCH --qos=soltis-b

date;hostname;pwd

### This script is used to mark the duplicates of the BAM files

INPUT=/blue/soltis/shan158538/Methylation/OutPut/SNP_calling/BWA_alignment/Alignment_output
OUTPUT=/blue/soltis/shan158538/Methylation/OutPut/SNP_calling/BWA_alignment/Alignment_output/Samtools_stats/Bam_sorted_marked_files
TMP=/orange/soltis/shan158538/methylpy_temp
METRICS=/blue/soltis/shan158538/Methylation/OutPut/SNP_calling/BWA_alignment/Alignment_output/marked_dup_metrics

module purge
module load picard/2.21.2
module load samtools/1.10
module load multiqc/1.7

cd ${INPUT}

for sample in 533 534 535 536 537 538 539 540 541 542 543 544 545 546 547 548
do

### Sort the bam file
echo "Sorting sample ${sample}"
java -Xmx20g -jar $HPC_PICARD_DIR/picard.jar \
	SortSam \
	I=Tpr_${sample}.bam \
	O=Tpr_${sample}_sorted.bam \
	SORT_ORDER=coordinate \
	TMP_DIR=${TMP} \
	VALIDATION_STRINGENCY=LENIENT
echo "Sample ${sample} has been sorted based on coordinate"

### Mark duplicates
java -Xmx20g -jar $HPC_PICARD_DIR/picard.jar \
	MarkDuplicates \
	I=Tpr_${sample}_sorted.bam \
	O=Tpr_${sample}_sorted_marked.bam \
	M=${METRICS}/Tpr_${sample}_marked_dup_metrics.txt \
	TMP_DIR=${TMP} \
	VALIDATION_STRINGENCY=LENIENT
echo "Duplicates in sample ${sample} have been marked!"

### COLLECT STATISTICS FROM SAM FILES
samtools \
	stats \
	-@ 10 \
	Tpr_${sample}_sorted_marked.bam > ${OUTPUT}/Tpr_${sample}_sorted_marked.bam_stats.txt
echo "The statistics of sample ${sample} has been analyzed"

done

### PRODUCE COMPREHENSIVE STATISTICS FROM ALIGNMENT FILES
cd ${OUTPUT}
multiqc \
	.

date
