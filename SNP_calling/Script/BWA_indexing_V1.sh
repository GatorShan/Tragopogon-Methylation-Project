#!/bin/bash

#SBATCH --job-name=BWA_indexing_V1
#SBATCH --mail-user=shan158538@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --output=BWA_indexing_V1_%j.out
#SBATCH --error=BWA_indexing_V1_%j.error
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2gb
#SBATCH --time=05:00:00
#SBATCH --qos=soltis

date;hostname;pwd

### This script is used to index the reference genome of T. dubius, which will be used in the following alignment analysis

REF=/blue/soltis/shan158538/Methylation/OutPut/SNP_calling/BWA_alignment/Tdub.V1.fasta

module purge
module load bwa/0.7.17

bwa \
	index \
	${REF}

date
