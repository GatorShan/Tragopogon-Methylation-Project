#!/bin/bash

#SBATCH --job-name=bismark_genome_prep_V1
#SBATCH --mail-user=shan158538@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --output=bismark_genome_prep_V1_%j.out
#SBATCH --error=bismark_genome_prep_V1_%j.error
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=10gb
#SBATCH --time=4-00:00:00
#SBATCH --qos=soltis-b

date;hostname;pwd

IN=/orange/soltis/shan158538/Methylation_output/bismark_genome_prep

module purge
module load bismark/0.22.3

bismark_genome_preparation \
	--parallel 2 \
	--verbose \
	${IN}

## default: create bisulfite indexes for use with Bowtie 2

date
