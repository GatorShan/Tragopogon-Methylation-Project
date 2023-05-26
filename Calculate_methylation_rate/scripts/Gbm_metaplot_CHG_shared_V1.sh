#!/bin/bash

#SBATCH --job-name=Gbm_metaplot_CHG_shared_V1
#SBATCH --mail-user=shan158538@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --output=Gbm_metaplot_CHG_shared_V1_%j.out
#SBATCH --error=Gbm_metaplot_CHG_shared_V1_%j.error
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=500mb
#SBATCH --time=1-00:00:00
#SBATCH --qos=soltis

date;hostname;pwd

module purge
module load python

MTYPE=CHG
GFF=/orange/soltis/shan158538/Methylation_output/feature_methylation/Tdub.V1.rm.RENAME.gff

for SAMPLE in DES1 S1 S2 S3 S4 S5; do
	echo "Processing sample ${SAMPLE}"
	IN=/orange/soltis/shan158538/Methylation_output/feature_methylation/gene_${MTYPE}_shared/${SAMPLE}
	
	## Using the default values: 1) number of bins: 20; 2) streamsize; 1000
	gbm_metaplot_pe_ss.V1.py \
		-m=${MTYPE} \
		-o=${IN}/${SAMPLE}_${MTYPE} \
		${GFF} \
		${IN} \
		${SAMPLE}
done

date