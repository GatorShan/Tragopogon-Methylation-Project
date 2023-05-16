#!/bin/bash

#SBATCH --job-name=bismark_genome_prep_V2
#SBATCH --mail-user=shan158538@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --output=bismark_genome_prep_V2_%j.out
#SBATCH --error=bismark_genome_prep_V2_%j.error
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=1gb
#SBATCH --time=01:00:00
#SBATCH --qos=soltis

date;hostname;pwd

## convert the lambda sequence

IN=/orange/soltis/shan158538/Methylation_output/bismark_lambda_genome_prep

module purge
module load bismark/0.22.3

bismark_genome_preparation \
	--parallel 2 \
	--verbose \
	${IN}

## default: create bisulfite indexes for use with Bowtie 2

date
