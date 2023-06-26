#!/bin/bash

#SBATCH --job-name=Gbm_metaplot_CHH_shared_V3
#SBATCH --mail-user=shan158538@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --output=Gbm_metaplot_CHH_shared_V3_%j.out
#SBATCH --error=Gbm_metaplot_CHH_shared_V3_%j.error
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=5gb
#SBATCH --time=03:00:00
#SBATCH --qos=soltis

date;hostname;pwd

module purge
module load python

MTYPE=CHH
GFF=/orange/soltis/shan158538/Methylation_output/feature_methylation/Tdub.V1.rm.RENAME.gff

for SAMPLE in DES1 S1 S2 S3 S4 S5; do
	echo "Processing sample ${SAMPLE}"
	IN=/orange/soltis/shan158538/Methylation_output/feature_methylation/gene_${MTYPE}_shared/${SAMPLE}
	
	## Using the default values: 1) number of bins: 20; 2) streamsize; 1000
	gbm_metaplot_pe_ss.V3.py \
		-m=${MTYPE} \
		-o=${IN}/${SAMPLE}_${MTYPE}_new2 \
		${GFF} \
		${IN} \
		${SAMPLE}
done

date
