#!/bin/bash

#SBATCH --job-name=Trim_galore.V1
#SBATCH --mail-user=shan158538@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --output=Trim_galore.V1_%j.out
#SBATCH --error=Trim_galore.V1_%j.error
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=300mb
#SBATCH --time=4-00:00:00
#SBATCH --qos=soltis

date;hostname;pwd

IN=/ufrc/soltis/shan158538/Methylation/Raw_data_Dec2019
OUT=/orange/soltis/shan158538/Methylation_output/trim_galore

module purge
#module load gcc/5.2.0
#module load pigz/2.4
#module load python3/3.6
module load trim_galore/0.5.0

cd ${IN}

trim_galore \
	--paired \
	--trim1 \
	S1_cat_R1.fastq.gz \
	S1_cat_R2.fastq.gz \
	--fastqc \
	--output_dir ${OUT} \
	
## 'trim1' trims one additional base pair from the 3' end of both reads, a step
## that is needed for subsequent alignments of completely overlapping long reads with Bowtie

date
