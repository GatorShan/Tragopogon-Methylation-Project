#!/bin/bash

#SBATCH --job-name=Trimmomatic_1.0
#SBATCH --mail-user=shan158538@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --output=Trimmomatic_1.0_%A-%a.out
#SBATCH --error=Trimmomatic_1.0_%A-%a.error
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=5gb
#SBATCH --time=08:00:00
#SBATCH --qos=soltis-b

date;hostname;pwd

IN=/orange/soltis/SequenceBackup/10KP/08.Tragopogon_pratensis/CL100078369_L01
OUT=/ufrc/soltis/shan158538/T.pratensis_Genome/OutPut/Trimmomatic

module load trimmomatic/0.36

trimmomatic \
	PE \
	-threads 4 \
	-phred33 \
	${IN}/CL100078369_L01_${SLURM_ARRAY_TASK_ID}_1.fq.gz \
	${IN}/CL100078369_L01_${SLURM_ARRAY_TASK_ID}_2.fq.gz \
	${OUT}/CL100078369_L01_${SLURM_ARRAY_TASK_ID}_1_trimmomatic_paired.fastq \
	${OUT}/CL100078369_L01_${SLURM_ARRAY_TASK_ID}_1_trimmomatic_unpaired.fastq \
	${OUT}/CL100078369_L01_${SLURM_ARRAY_TASK_ID}_2_trimmomatic_paired.fastq \
	${OUT}/CL100078369_L01_${SLURM_ARRAY_TASK_ID}_2_trimmomatic_unpaired.fastq \
	ILLUMINACLIP:${OUT}/BGIadapter.fa:2:35:4:12:true \
	LEADING:3 TRAILING:3 SLIDINGWINDOW:5:15 MINLEN:50

date
