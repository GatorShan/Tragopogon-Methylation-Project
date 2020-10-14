#!/bin/bash

#SBATCH --job-name=AddReadGroup.V1
#SBATCH --mail-user=shan158538@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --output=AddReadGroup.V1_%j.out
#SBATCH --error=AddReadGroup.V1_%j.error
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=5gb
#SBATCH --time=01-00:00:00
#SBATCH --qos=soltis

date;hostname;pwd

### This script is used to add read group names to each bam file

INPUT=/blue/soltis/shan158538/Methylation/OutPut/SNP_calling/BWA_alignment/Alignment_output/Bam_sorted_marked_files

module purge
module load picard/2.21.2

for sample in 533 534 535 536 537 538 539 540 541 542 543 544 545 546 547 548
do

java -Xmx3g -jar $HPC_PICARD_DIR/picard.jar \
	AddOrReplaceReadGroups \
	I=${INPUT}/Tpr_${sample}_sorted_marked.bam \
	O=${INPUT}/Tpr_${sample}_sorted_marked_AddReadGroup.bam \
	RGID= CL100078369.L01.${sample}\
	RGLB= LIB-${sample}\
	RGPL=BGISEQ \
	RGPU= CL100078369.L01.${sample}\
	RGSM=Tpr

echo "Read group name has been added to sample ${sample}"

done

date
