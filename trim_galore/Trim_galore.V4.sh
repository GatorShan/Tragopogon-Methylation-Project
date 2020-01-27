#!/bin/bash

#SBATCH --job-name=Trim_galore.V4
#SBATCH --mail-user=shan158538@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --output=Trim_galore.V4_%j.out
#SBATCH --error=Trim_galore.V4_%j.error
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=300mb
#SBATCH --time=4-00:00:00
#SBATCH --qos=soltis-b

date;hostname;pwd

IN=/orange/soltis/SequenceBackup/ShengchenShan_MethylSeq_20180731
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
	HMCWKCCXY_s8_1_4981-LF_17_SL334590.fastq.gz \
	HMCWKCCXY_s8_2_4981-LF_17_SL334590.fastq.gz \
	--output_dir ${OUT}
	
## 'trim1' trims one additional base pair from the 3' end of both reads, a step
## that is needed for subsequent alignments of completely overlapping long reads with Bowtie

date
