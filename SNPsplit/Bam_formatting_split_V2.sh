#!/bin/bash

#SBATCH --job-name=Bam_formatting_split_V2
#SBATCH --mail-user=shan158538@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --output=Bam_formatting_split_V2_%j.out
#SBATCH --error=Bam_formatting_split_V2_%j.error
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=5gb
#SBATCH --time=1-00:00:00
#SBATCH --qos=soltis

date;hostname;pwd

module purge
module load samtools/1.9

IN=/blue/soltis/shan158538/Methylation/OutPut/SNPsplit_alignment/S5_dedup

cd ${IN}

# Sort by position
samtools sort \
	-o S5_cat_R1_val_1_bismark_bt2_pe_nameSorted.deduplicated_PosSorted.bam \
	S5_cat_R1_val_1_bismark_bt2_pe_nameSorted.deduplicated.bam

# Index
samtools index \
	S5_cat_R1_val_1_bismark_bt2_pe_nameSorted.deduplicated_PosSorted.bam

date