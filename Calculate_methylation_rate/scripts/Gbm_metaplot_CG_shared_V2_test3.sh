#!/bin/bash

#SBATCH --job-name=Gbm_metaplot_CG_shared_V2_test3
#SBATCH --mail-user=shan158538@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --output=Gbm_metaplot_CG_shared_V2_test3_%j.out
#SBATCH --error=Gbm_metaplot_CG_shared_V2_test3_%j.error
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=500mb
#SBATCH --time=1-00:00:00
#SBATCH --qos=soltis

date;hostname;pwd

module purge
module load python

MTYPE=CG
GFF=/orange/soltis/shan158538/Methylation_output/feature_methylation/test_3.gff

for SAMPLE in DES1; do
	echo "Processing sample ${SAMPLE}"
	IN=/orange/soltis/shan158538/Methylation_output/feature_methylation/gene_${MTYPE}_shared/${SAMPLE}_test3
	
	## Using the default values: 1) number of bins: 20; 2) streamsize; 1000
	gbm_metaplot_pe_ss.V3.py \
		-m=${MTYPE} \
		-o=${IN}/${SAMPLE}_${MTYPE}_new \
		${GFF} \
		${IN} \
		${SAMPLE}
done

date
