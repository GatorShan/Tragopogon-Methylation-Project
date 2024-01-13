#!/bin/bash

#SBATCH --job-name=Copia_metaplot_CHH_shared_V3
#SBATCH --mail-user=shan158538@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --output=Copia_metaplot_CHH_shared_V3_%j.out
#SBATCH --error=Copia_metaplot_CHH_shared_V3_%j.error
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=500mb
#SBATCH --time=04:00:00
#SBATCH --qos=soltis

date;hostname;pwd

module purge
module load python

MTYPE=CHH
TE_TYPE=Copia
GFF=/blue/soltis/shan158538/Methylation/OutPut/TE_annotation_new/repeat_annotation_combined_Copia.final.gff

for SAMPLE in DES1 S1 S2 S3 S4 S5; do
	echo "Processing sample ${SAMPLE}"
	IN=/orange/soltis/shan158538/Methylation_output/feature_methylation/gene_${MTYPE}_shared/${SAMPLE}
	
	## Using the new values: 1) number of bins: 20; 2) streamsize; 2000
	TE_metaplot_pe_ss.V4.py \
		-m=${MTYPE} \
		-o=${IN}/${SAMPLE}_${MTYPE}_${TE_TYPE}_new2 \
		${GFF} \
		${IN} \
		${SAMPLE}
done

date
