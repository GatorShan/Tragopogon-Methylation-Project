#!/bin/bash

#SBATCH --job-name=Bismark_bam_formatting_V2
#SBATCH --mail-user=shan158538@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --output=Bismark_bam_formatting_V2_%j.out
#SBATCH --error=Bismark_bam_formatting_V2_%j.error
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=30gb
#SBATCH --time=2-00:00:00
#SBATCH --qos=soltis

date;hostname;pwd

module purge
module load samtools/1.9

IN=/blue/soltis/shan158538/Methylation/OutPut/SNPsplit_output

cd ${IN}

for FILE in `ls *genome*bam`; do
	
	echo "$FILE is being processed"
	
	# Sort by position
	samtools sort \
		-o ${FILE%.*}_PosSorted.bam \
		$FILE
	echo -e "\t$FILE has been sorted by position"

	# Index
	samtools index \
		${FILE%.*}_PosSorted.bam
	echo -e "\t${FILE%.*}_PosSorted.bam file has been indexed"

done

date
