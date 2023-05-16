#!/bin/bash

#SBATCH --job-name=bismark2bedGraph_S2_CHG_V1
#SBATCH --mail-user=shan158538@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --output=bismark2bedGraph_S2_CHG_V1_%j.out
#SBATCH --error=bismark2bedGraph_S2_CHG_V1_%j.error
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=5gb
#SBATCH --time=2-00:00:00
#SBATCH --qos=soltis

date;hostname;pwd

IN=/orange/soltis/shan158538/Methylation_output/bismark_methylation_extractor/S2

module purge
module load bismark/0.22.3

bismark2bedGraph \
	--output S2_CHG \
	--dir ${IN} \
	--CX_context \
	--buffer_size 4G \
	${IN}/CHG_OB_S2_cat_R1_val_1_bismark_bt2_pe_nameSorted.deduplicated.txt.gz \
	${IN}/CHG_OT_S2_cat_R1_val_1_bismark_bt2_pe_nameSorted.deduplicated.txt.gz

date
