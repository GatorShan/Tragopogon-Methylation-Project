#!/bin/bash

#SBATCH --job-name=methylpy_dataprocessing_Tpr_3058-4-10_v2
#SBATCH --mail-user=shan158538@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --output=methylpy_dataprocessing_Tpr_3058-4-10_v2_%j.out
#SBATCH --error=methylpy_dataprocessing_Tpr_3058-4-10_v2_%j.error
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=30gb
#SBATCH --time=04-00:00:00
#SBATCH --qos=soltis-b

date;hostname;pwd

DB=/ufrc/soltis/shan158538/Methylation/OutPut/Build_reference
IN=/ufrc/soltis/shan158538/Methylation/Raw_data_Dec2019
OUT=/orange/soltis/shan158538/methylpy_output/methylpy_dataprocessing_Tpr_3058-4-10_V2

echo "TMPDIR is: ${TMPDIR}"

module purge
module load gcc/5.2.0
module load methylpy/1.2.9
module load bowtie2/2.3.5.1

methylpy paired-end-pipeline \
	--read1-files ${IN}/S3_cat_R1.fastq.gz \
	--read2-files ${IN}/S3_cat_R2.fastq.gz \
	--sample T.pratensis_3058-4-10 \
	--forward-ref ${DB}/Tdubius_Ref_f \
	--reverse-ref ${DB}/Tdubius_Ref_r \
	--ref-fasta ${DB}/Tdub_lambda.fasta \
	--path-to-output ${OUT} \
	--num-procs 16 \
	--trim-reads True \
	--bowtie2 True \
	--remove-clonal True \
	--path-to-picard="$HPC_PICARD_DIR/" \
	--adapter-seq-read1 AGATCGGAAGAGCACACGTCTGAAC \
	--adapter-seq-read2 AGATCGGAAGAGCGTCGTGTAGGGA \
	--binom-test True \
	--unmethylated-control "NC_001416.1_lambda:" \
	--keep-temp-files True

date
