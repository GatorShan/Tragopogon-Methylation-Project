#!/bin/bash

#SBATCH --job-name=methylpy_dataprocessing_Tdu_3040-6-2_v1_Picard
#SBATCH --mail-user=shan158538@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --output=methylpy_dataprocessing_Tdu_3040-6-2_v1_Picard_%j.out
#SBATCH --error=methylpy_dataprocessing_Tdu_3040-6-2_v1_Picard_%j.error
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=50gb
#SBATCH --time=4-00:00:00
#SBATCH --qos=soltis-b

# use the output from job 3.1; run picard separately

date;hostname;pwd

module purge
module load picard/2.18.3

FILE=/orange/soltis/shan158538/methylpy_output/methylpy_dataprocessing_Tdu_3040-6-2_V1
TMP=/orange/soltis/shan158538/methylpy_temp

java -Xmx20g -jar \
	/apps/picard/2.18.3/picard.jar \
	MarkDuplicates \
	INPUT=${FILE}/T.dubius_3040-6-2_libA_processed_reads.bam \
	OUTPUT=${FILE}/T.dubius_3040-6-2_libA_processed_reads_no_clonal.bam \
	ASSUME_SORTED=true \
	REMOVE_DUPLICATES=true \
	METRICS_FILE=${FILE}/T.dubius_3040-6-2_libA.metric \
	VALIDATION_STRINGENCY=LENIENT \
	QUIET=true \
	TMP_DIR=${TMP}

date
