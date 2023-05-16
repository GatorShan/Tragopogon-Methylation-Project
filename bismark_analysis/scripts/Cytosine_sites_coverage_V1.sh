#!/bin/bash
#SBATCH --job-name=Cytosine_sites_coverage_V1
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=shan158538@ufl.edu	
#SBATCH --ntasks=1
#SBATCH --mem=60gb
#SBATCH --time=01:00:00
#SBATCH --output=Cytosine_sites_coverage_V1_%j.out 
#SBATCH --error=Cytosine_sites_coverage_V1_%j.error

date; hostname; pwd

module purge
module load python3

SCRIPT=/blue/soltis/shan158538/scripts/Cytosine_coverage_V1.py
LIST=("S1" "S2" "S3" "S4" "S5")

for SAMPLE in "${LIST[@]}"
do
    DIRECTORY=/orange/soltis/shan158538/Methylation_output/bismark_methylation_extractor/${SAMPLE}
    INFILE_1=${DIRECTORY}/${SAMPLE}_cat_R1_val_1_bismark_bt2_pe_nameSorted.deduplicated.bismark.cov
    INFILE_2=${DIRECTORY}/${SAMPLE}_CHG.gz.bismark.cov
    INFILE_3=${DIRECTORY}/${SAMPLE}_CHH.gz.bismark.cov

    echo ${SAMPLE}
    echo "All cytosine sites mean coverage"
    python ${SCRIPT} ${INFILE_1} ${INFILE_2} ${INFILE_3}
done


date
