#!/bin/bash
#SBATCH --job-name=Cytosine_sites_coverage_V2
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=shan158538@ufl.edu	
#SBATCH --ntasks=1
#SBATCH --mem=60gb
#SBATCH --time=01:00:00
#SBATCH --output=Cytosine_sites_coverage_V2_%j.out 
#SBATCH --error=Cytosine_sites_coverage_V2_%j.error

date; hostname; pwd

module purge
module load python3

SCRIPT=/blue/soltis/shan158538/scripts/Cytosine_coverage_V1.py

DIRECTORY=/orange/soltis/shan158538/Methylation_output/bismark_methylation_extractor/DES1
INFILE_1=${DIRECTORY}/HMCWKCCXY_s8_1_4981-LF_17_SL334590_val_1_bismark_bt2_pe_nameSorted.deduplicated.bismark.cov
INFILE_2=${DIRECTORY}/DES1_CHG.gz.bismark.cov
INFILE_3=${DIRECTORY}/DES1_CHH.gz.bismark.cov

echo "All cytosine sites mean coverage for sample DES1"
python ${SCRIPT} ${INFILE_1} ${INFILE_2} ${INFILE_3}

date
