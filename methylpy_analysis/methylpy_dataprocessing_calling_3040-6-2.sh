#!/bin/bash

#SBATCH --job-name=methylpy_dataprocessing_calling_3040-6-2
#SBATCH --mail-user=shan158538@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --output=methylpy_dataprocessing_calling_3040-6-2_%j.out
#SBATCH --error=methylpy_dataprocessing_calling_3040-6-2_%j.error
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=25gb
#SBATCH --time=2-00:00:00
#SBATCH --qos=soltis-b

date;hostname;pwd

#DB=/ufrc/soltis/shan158538/Methylation/OutPut/Build_reference
#IN=/orange/soltis/SequenceBackup/ShengchenShan_MethylSeq_20180731
#OUT=/ufrc/soltis/shan158538/Methylation/OutPut/methylpy_processing_data_3.1
CODE=/ufrc/soltis/shan158538/Methylation/SLURM

echo "TMPDIR is: ${TMPDIR}"

module purge
module load gcc/5.2.0
module load methylpy/1.2.9
#module load bowtie2/2.3.5.1
#module load picard/2.18.3

python ${CODE}/run_methylpy_3040-6-2.py

date
