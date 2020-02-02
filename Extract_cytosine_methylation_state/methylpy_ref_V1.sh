#!/bin/bash

#SBATCH --job-name=methylpy_ref_V1
#SBATCH --mail-user=shan158538@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --output=methylpy_ref_V1_%j.out
#SBATCH --error=methylpy_ref_V1_%j.error
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=12gb
#SBATCH --time=03:00:00

date;hostname;pwd

IN=/orange/soltis/shan158538/Methylation_output/methylpy_ref

cd ${IN}

module purge
module load gcc/5.2.0
module load methylpy/1.2.9

methylpy build-reference \
  --input-files Tdub.V1.fasta \
  --output-prefix Tdu_ref \
  --bowtie2 True \
  --parallel True

date
