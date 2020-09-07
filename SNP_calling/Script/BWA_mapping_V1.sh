#!/bin/bash

#SBATCH --job-name=BWA_mapping_V1
#SBATCH --mail-user=shan158538@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --output=BWA_mapping_V1_%j.out
#SBATCH --error=BWA_mapping_V1_%j.error
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5
#SBATCH --mem=8gb
#SBATCH --time=04-00:00:00
#SBATCH --qos=soltis-b

date;hostname;pwd

### This script is used to map trimmed paired T. pratensis reads to the Tdu reference genome

REF=/blue/soltis/shan158538/Methylation/OutPut/SNP_calling/BWA_alignment/Tdub.V1.fasta
INPUT=/orange/soltis/T.pratensis_genome/TrimmedReads
OUTPUT=/blue/soltis/shan158538/Methylation/OutPut/SNP_calling/BWA_alignment/Alignment_output

module purge
module load bwa/0.7.17

for sample in 533 534 535 536 537 538 539 540 541 542 543 544 545 546 547 548
do

echo "Processing sample ${sample}"

bwa \
	mem \
	-t 5 \
	${REF} \
	${INPUT}/CL100078369_L01_${sample}_1_trimmomatic_paired.fastq \
	${INPUT}/CL100078369_L01_${sample}_2_trimmomatic_paired.fastq \
	> ${OUTPUT}/Tpr_${sample}.sam

echo "sample ${sample} has been mapped to the Tdu ref genome!"

done

date
