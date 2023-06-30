#!/bin/bash

#SBATCH --job-name=TE_metaplot_CHG_shared_V2
#SBATCH --mail-user=shan158538@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --output=TE_metaplot_CHG_shared_V2_%j.out
#SBATCH --error=TE_metaplot_CHG_shared_V2_%j.error
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=500mb
#SBATCH --time=01:00:00
#SBATCH --qos=soltis

date;hostname;pwd

module purge
module load python

MTYPE=CHG
TE_TYPE=TE
GFF=/blue/soltis/shan158538/Methylation/OutPut/TE_annotation/Tdub2_rnd3_50aa.TE.maker.rm.RENAME.gff

for SAMPLE in DES1 S1 S2 S3 S4 S5; do
	echo "Processing sample ${SAMPLE}"
	IN=/orange/soltis/shan158538/Methylation_output/feature_methylation/gene_${MTYPE}_shared/${SAMPLE}
	
	## Using the default values: 1) number of bins: 20; 2) streamsize; 1000
	TE_metaplot_pe_ss.V3.py \
		-m=${MTYPE} \
		-o=${IN}/${SAMPLE}_${MTYPE}_${TE_TYPE}_new \
		${GFF} \
		${IN} \
		${SAMPLE}
done

date
