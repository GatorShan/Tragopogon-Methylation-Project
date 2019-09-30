# *Tragopogon* Methylation Project
## Discription
This analysis includs one library from diploid *Tragopogon dibius* (3060-1-4).
It was sequenced on platform Hiseq X, paired-end 150 bp.
The reference genome that I used is from Xiaoxian: `Tdub.V1.fasta`.
[methylpy](https://github.com/yupenghe/methylpy) was used to perform the analysis.

## Method
**Build reference**

1. Combine Tdub.V1.fasta and lambda.fasta (NC_001416.1_lambda) files

  ```bash
  cat Tdub.V1.fasta lambda.fasta >> Tdub_lambda.fasta
  ```

2. Build reference
  ```bash
  #!/bin/bash
  
  #SBATCH --job-name=Build_reference_1.2
  #SBATCH --mail-user=shan158538@ufl.edu
  #SBATCH --mail-type=FAIL,END
  #SBATCH --output=Build_reference_1.2_%j.out
  #SBATCH --error=Build_reference_1.2_%j.error
  #SBATCH --nodes=1
  #SBATCH --ntasks=1
  #SBATCH --cpus-per-task=2
  #SBATCH --mem=12gb
  #SBATCH --time=03:00:00
  
  date;hostname;pwd
  
  module purge
  module load gcc/5.2.0
  module load methylpy/1.2.9
  
  methylpy build-reference \
    --input-files Tdub_lambda.fasta \
    --output-prefix Tdubius_Ref \
    --bowtie2 True
  date
  ```


**Process data**

1. Change the directory for temporary files in ~/.bashrc
  ```bash
  export TMPDIR=/orange/soltis/shan158538/methylpy_temp
  ```

2.  
  ```bash
  source .bashrc
  ```

3. Use the following script (methylpy_processing_data_2.2.sh) to run methylpy

  ```bash
  #!/bin/bash

  #SBATCH --job-name=methylpy_processing_data_2.2
  #SBATCH --mail-user=shan158538@ufl.edu
  #SBATCH --mail-type=FAIL,END
  #SBATCH --output=methylpy_processing_data_2.2_%j.out
  #SBATCH --error=methylpy_processing_data_2.2_%j.error
  #SBATCH --nodes=1
  #SBATCH --ntasks=1
  #SBATCH --cpus-per-task=16
  #SBATCH --mem=64gb
  #SBATCH --time=10-00:00:00
  
  date;hostname;pwd
  
  DB=/ufrc/soltis/shan158538/Methylation/OutPut/Build_reference
  IN=/orange/soltis/SequenceBackup/ShengchenShan_MethylSeq_20180731
  OUT=/ufrc/soltis/shan158538/Methylation/OutPut/methylpy_processing_data_2.2
  
  echo "TMPDIR is: ${TMPDIR}"
  echo "TMP_DIR is: ${TMP_DIR}"
  
  module purge
  module load gcc/5.2.0
  module load methylpy/1.2.9
  module load bowtie2/2.3.5.1
  module load picard/2.18.3
  
  methylpy paired-end-pipeline \
	--read1-files ${IN}/HMCWKCCXY_s8_1_4981-LF_17_SL334590.fastq.gz \
	--read2-files ${IN}/HMCWKCCXY_s8_2_4981-LF_17_SL334590.fastq.gz \
	--sample T.dubius_3060-1-4 \
	--forward-ref ${DB}/Tdubius_Ref_f \
	--reverse-ref ${DB}/Tdubius_Ref_r \
	--ref-fasta ${DB}/Tdub_lambda.fasta \
	--path-to-output ${OUT} \
	--num-procs 16 \
	--trim-reads True \
	--bowtie2 True \
	--remove-clonal True \
	--path-to-picard="$HPC_PICARD_DIR/" \
	--adapter-seq-read1 AGATCGGAAGAGCACACGTCTGAAC \
	--adapter-seq-read2 AGATCGGAAGAGCGTCGTGTAGGGA \
	--binom-test True \
	--unmethylated-control "NC_001416.1_lambda:"

  date
  ```

4. Error messages will show up

  *Exception in thread "main" htsjdk.samtools.util.RuntimeIOException: java.io.IOException: No space left on device*

  *Caused by: java.io.IOException: No space left on device*

5. The output file is:

  Begin splitting reads for libA
  Wed Aug 14 22:40:17 2019

  Begin trimming reads for libA
  Thu Aug 15 01:38:03 2019

  Begin converting reads for libA
  Thu Aug 15 01:55:28 2019

  Begin Running Bowtie2 for libA
  Thu Aug 15 02:06:06 2019

  Processing forward strand hits
  Thu Aug 15 16:39:00 2019

  Processing reverse strand hits
  Fri Aug 16 07:31:57 2019

  Finding multimappers
  Fri Aug 16 08:36:34 2019
  
6. Picard will run seperately to remove duplicates.
  - Input: T.dubius_3060-1-4_libA_processed_reads.bam
  - Output: T.dubius_3060-1-4_libA_processed_reads_no_clonal.bam
  - Location: /ufrc/soltis/shan158538/Methylation/OutPut/methylpy_processing_data_3.1

  ```bash
  #!/bin/bash

  #SBATCH --job-name=methylpy_processing_data_3.1_picard
  #SBATCH --mail-user=shan158538@ufl.edu
  #SBATCH --mail-type=FAIL,END
  #SBATCH --output=methylpy_processing_data_3.1_picard_%j.out
  #SBATCH --error=methylpy_processing_data_3.1_picard_%j.error
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
  
  java -Xmx20g -jar \
	/apps/picard/2.18.3//picard.jar \
	MarkDuplicates \
	INPUT=/ufrc/soltis/shan158538/Methylation/OutPut/methylpy_processing_data_3.1/T.dubius_3060-1-4_libA_processed_reads.bam \
	OUTPUT=/ufrc/soltis/shan158538/Methylation/OutPut/methylpy_processing_data_3.1/T.dubius_3060-1-4_libA_processed_reads_no_clonal.bam \
	ASSUME_SORTED=true \
	REMOVE_DUPLICATES=true \
	METRICS_FILE=/ufrc/soltis/shan158538/Methylation/OutPut/methylpy_processing_data_3.1/T.dubius_3060-1-4_libA.metric \
	VALIDATION_STRINGENCY=LENIENT \
	QUIET=true \
	TMP_DIR=/orange/soltis/shan158538/methylpy_temp

  date
  ```

7. Finish up the following analysis

  ```bash
  #!/bin/bash

  #SBATCH --job-name=methylpy_processing_data_4.0_calling
  #SBATCH --mail-user=shan158538@ufl.edu
  #SBATCH --mail-type=FAIL,END
  #SBATCH --output=methylpy_processing_data_4.0_calling_%j.out
  #SBATCH --error=methylpy_processing_data_4.0_calling_%j.error
  #SBATCH --nodes=1
  #SBATCH --ntasks=1
  #SBATCH --cpus-per-task=16
  #SBATCH --mem=25gb
  #SBATCH --time=2-00:00:00

  date;hostname;pwd

  CODE=/ufrc/soltis/shan158538/Methylation/SLURM

  echo "TMPDIR is: ${TMPDIR}"

  module purge
  module load gcc/5.2.0
  module load methylpy/1.2.9

  python ${CODE}/run_methylpy_1.0.py

  date
  ```





