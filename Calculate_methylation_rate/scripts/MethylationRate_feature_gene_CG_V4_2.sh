#!/bin/bash

#SBATCH --job-name=MethylationRate_feature_gene_CG_V4_2
#SBATCH --mail-user=shan158538@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --output=MethylationRate_feature_gene_CG_V4_2_%j.out
#SBATCH --error=MethylationRate_feature_gene_CG_V4_2_%j.error
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=3
#SBATCH --mem=500mb
#SBATCH --time=1-00:00:00
#SBATCH --qos=soltis

date;hostname;pwd

module purge
module load python

MTYPE=CG
GFF=/orange/soltis/shan158538/Methylation_output/feature_methylation/Tdub.V1.rm.RENAME.gff

for SAMPLE in S1 S2 S3 S4 S5; do
	echo "Processing sample ${SAMPLE}"
	IN=/orange/soltis/shan158538/Methylation_output/feature_methylation/gene_${MTYPE}/${SAMPLE}
	
	## Generate a list of input allC files; remove the "," at the beginning of the list
	LIST=""
	cd ${IN}
	## Make sure ${IN} is not in the for loop!
	for FILE in *.tsv; do
		SCAFFOLD=$(echo ${FILE} | cut -f 3 -d "_" | cut -f 1 -d ".")
		LIST="${LIST},${SCAFFOLD}"
	done
	LIST=$(echo ${LIST} | cut -c 2-)

	## Calculate the methylation rate of gene regions
	feature_methylation.ss.V1.py \
		-o=${IN}/${SAMPLE}_gene_${MTYPE} \
		-c=${LIST} \
		-p=${SLURM_CPUS_PER_TASK} \
		-m=${MTYPE} \
		-g \
		${GFF} \
		${IN} \
		${SAMPLE}
done

date
