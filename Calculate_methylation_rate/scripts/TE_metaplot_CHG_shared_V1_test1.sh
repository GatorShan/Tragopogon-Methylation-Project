#!/bin/bash

#SBATCH --job-name=TE_metaplot_CHG_shared_V1_test1
#SBATCH --mail-user=shan158538@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --output=TE_metaplot_CHG_shared_V1_test1_%j.out
#SBATCH --error=TE_metaplot_CHG_shared_V1_test1_%j.error
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=500mb
#SBATCH --time=00:10:00
#SBATCH --qos=soltis

date;hostname;pwd

module purge
module load python

MTYPE=CHG
TE_TYPE=TE
GFF=/blue/soltis/shan158538/Methylation/OutPut/TE_annotation/TE_test1.gff

for SAMPLE in DES1; do
	echo "Processing sample ${SAMPLE}"
	IN=/orange/soltis/shan158538/Methylation_output/feature_methylation/gene_${MTYPE}_shared/${SAMPLE}_test1
	
	## Using the default values: 1) number of bins: 20; 2) streamsize; 1000
	TE_metaplot_pe_ss.V2.py \
		-m=${MTYPE} \
		-o=${IN}/${SAMPLE}_${MTYPE}_${TE_TYPE}_test1 \
		${GFF} \
		${IN} \
		${SAMPLE}
done

date
