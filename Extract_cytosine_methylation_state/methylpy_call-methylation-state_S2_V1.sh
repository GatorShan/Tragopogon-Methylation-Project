#!/bin/bash

#SBATCH --job-name=methylpy_call-methylation-state_S2_V1
#SBATCH --mail-user=shan158538@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --output=methylpy_call-methylation-state_S2_V1_%j.out
#SBATCH --error=methylpy_call-methylation-state_S2_V1_%j.error
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=12gb
#SBATCH --time=1-00:00:00
#SBATCH --qos=soltis

date;hostname;pwd

IN=/orange/soltis/shan158538/Methylation_output/bismark_deduplicate
REF=/orange/soltis/shan158538/Methylation_output/methylpy_ref
OUT=/orange/soltis/shan158538/Methylation_output/call_methylation_state

module purge
module load gcc/5.2.0
module load methylpy/1.2.9

methylpy call-methylation-state \
	--input-file ${IN}/S2_cat_R1_val_1_bismark_bt2_pe_nameSorted.deduplicated.bam \
	--paired-end True \
	--sample Tpr_3058-1-2 \
	--ref-fasta ${REF}/Tdub.V1.fasta \
	--path-to-output ${OUT} \
	--num-procs 2 \
	--min-cov 2 \
	--binom-test True

date
