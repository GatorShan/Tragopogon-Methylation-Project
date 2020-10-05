#!/bin/bash

#SBATCH --job-name=WgsMetrics.V1
#SBATCH --mail-user=shan158538@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --output=WgsMetrics.V1_%j.out
#SBATCH --error=WgsMetrics.V1_%j.error
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=5gb
#SBATCH --time=01-00:00:00
#SBATCH --qos=soltis

date;hostname;pwd

### This script is used to collect WGS metrics

INPUT=/blue/soltis/shan158538/Methylation/OutPut/SNP_calling/BWA_alignment/Alignment_output/Bam_sorted_marked_files
REF=/blue/soltis/shan158538/Methylation/OutPut/SNP_calling/BWA_alignment/Tdub.V1.fasta

module purge
module load picard/2.21.2

java -Xmx3g -jar $HPC_PICARD_DIR/picard.jar \
	CollectWgsMetrics \
	I=${INPUT}/Tpr_combined_sorted_marked.bam \
	O=${INPUT}/collect_wgs_metrics.txt \
	R=${REF} \
	READ_LENGTH=100

date
