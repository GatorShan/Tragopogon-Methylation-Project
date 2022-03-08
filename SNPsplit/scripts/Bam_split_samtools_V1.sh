#!/bin/bash

#SBATCH --job-name=Bam_split_samtools_V1
#SBATCH --mail-user=shan158538@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --output=Bam_split_samtools_V1_%j.out
#SBATCH --error=Bam_split_samtools_V1_%j.error
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=5gb
#SBATCH --time=1-00:00:00
#SBATCH --qos=soltis

date;hostname;pwd

module purge
module load samtools/1.9

IN=/blue/soltis/shan158538/Methylation/OutPut/SNPsplit_alignment/S4_dedup
OUT=/blue/soltis/shan158538/Methylation/OutPut/SNPsplit_alignment/S4_dedup/split

for NUM in {1..12310}; do
	
	samtools view \
		-b \
		${IN}/S4_cat_R1_val_1_bismark_bt2_pe_nameSorted.deduplicated_PosSorted.bam \
		Tdub_V1_scaffold_${NUM} > ${OUT}/S4_deduplicated_scaffold_${NUM}.bam
	echo "Scaffold $NUM has been processed"

done

date
