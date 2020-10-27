# SNP Calling
The T. pratensis (2x) genome sequencing reads are mapped to the reference genome of T. dubius (2x) -- SNPs between the two species are then identified.

## 1. Mapping
### 1.1 Mapping T. pratensis reads to T. dubius reference genome using BWA
Script `BWA_indexing_V1.sh` was used to index the reference genome. Script `BWA_mapping_V1.sh` was used for mapping.

When mapping the reads, I used only the paired-end reads and discarded the unpaired reads (~1.7%). As only a small proportion of reads is unpaired, this is good based on the discussion [here](https://www.biostars.org/p/140318/)

~90% reads mapped to the reference. The statistics of the mapping results could be found [here](https://github.com/GatorShan/Tragopogon-Methylation-Project/blob/master/SNP_calling/MultiQC_results/Bam_file_multiqc_report.html).

Resources:
  - [BWA manual](http://bio-bwa.sourceforge.net/bwa.shtml)
  - [Online tutorial](https://hbctraining.github.io/In-depth-NGS-Data-Analysis-Course/sessionVI/lessons/01_alignment.html)

### 1.2 Convert SAM file to BAM file

Script below was used.
```bash
module purge
module load samtools/1.10

### CONVERT SAM FILE TO BAM FILE
### -b: output BAM
### -S: input format is auto-detected; must be specified as BAM file is expected
### -@: number of threads
samtools \
        view \
        -bS \
        ${INPUT}/Tpr_${sample}.sam > ${INPUT}/Tpr_${sample}.bam \
        -@ 4
```

## 2. Mark duplicates
### 2.1 Sort the BAM file and mark the duplicates
Script `Picard_MarkDuplicates_V3.sh` was used. More information about Picard Markduplicates function can be found [here](https://broadinstitute.github.io/picard/command-line-overview.html#MarkDuplicates)

Tips:
  - Requrest a few more memory (25 Gb) than what is allocated to the java script (-Xmx 20 Gb)
  - Set up a new TMP directory at a local directory (error with the default TMP directory)

![Image_1](https://github.com/GatorShan/Tragopogon-Methylation-Project/blob/master/SNP_calling/images/Bam_sorted_marked_file_multiqc_report_samtools_alignment_plot.png)

![Image_2](https://github.com/GatorShan/Tragopogon-Methylation-Project/blob/master/SNP_calling/images/Bam_sorted_marked_file_multiqc_report_gerneral_statistics.png)

![Image_3](https://github.com/GatorShan/Tragopogon-Methylation-Project/blob/master/SNP_calling/images/Bam_sorted_marked_file_multiqc_report_alignment_matrics.png)

The percentage of duplicate reads is **5-6%**.
![Image_4](https://github.com/GatorShan/Tragopogon-Methylation-Project/blob/master/SNP_calling/images/Duplicate_percentage.png)

### 2.2 Index the BAM files
Then, the BAM files were indexed using the script `Index_bam.V1.sh` (only sorted and indexed bam files can be visualized by IGV).
  - Run graphical user interface (GUI) applications on [HiPerGator](https://help.rc.ufl.edu/doc/GUI_Programs)
  ```bash
  module load gui/2
  
  ### Start a GUI session: --module: environment module(s) to load; -e: application executable; -m: memory (GB); -t: time (hour)
  gui start --module igv -e igv.sh -m 4 -t 2
  
  ### show running sessions; copy a connection URL
  gui show
  
  ### stop a GUI session
  gui stop
  ```
  - IGV
    - Navigate to the dirctory contaning the bam and bai files
    - Launch the GUI program as shown above
    - Load the reference genome
    - Load the bam file
![Image_5](https://github.com/GatorShan/Tragopogon-Methylation-Project/blob/master/SNP_calling/images/IGV_bam_example.png)

### 2.3 Combine individual BAM files
Script `Merge_bam.V1.sh` was used to merge all sorted and dup-marked SAM files into a single SAM file, and then index the output SAM file.

### 2.4 Validates the BAM file
Script `Validate_sam_file.V1.sh` was used to validate the combined BAM file. More information about BAM file validation can be found [here](https://gatk.broadinstitute.org/hc/en-us/articles/360036854731-ValidateSamFile-Picard-) and [here](https://gatk.broadinstitute.org/hc/en-us/articles/360035891231-Errors-in-SAM-or-BAM-files-can-be-diagnosed-with-ValidateSamFile). To generate more detailed information of the error message script `Validate_sam_file.V2.sh` was used.

The summary message is shown below:
```
## HISTOGRAM	java.lang.String
Error Type	Count
ERROR:MISSING_READ_GROUP	1
WARNING:RECORD_MISSING_READ_GROUP	1490169017
```

### 2.5 Coverage
Script `WgsMetrics.V1.sh` was used to collect metrics about coverage and performance of whole genome sequencing (WGS) experiments. To interpret the results use the description shown [here](https://broadinstitute.github.io/picard/picard-metric-definitions.html#CollectWgsMetrics.WgsMetrics).

| Description | Result |
| -- | -- |
| The number of non-N bases in the genome reference over which coverage will be evaluated | 820,481,643 |
| The mean coverage | 43.2 |
| The standard deviation of coverage | 58.7 |
| The median coverage | 21 |
| The total fraction of aligned bases excluded due to all filters | 20.1% |
| The fraction of bases that attained at least 20X sequence coverage in post-filtering bases | 54.0% |
| The fraction of bases that attained at least 30X sequence coverage in post-filtering bases | 48.5% |

### 2.6 Add read group name to individual Bam file
Script `AddReadGroup.V2.sh` was used to add read group names. More info of read group can be found [here](https://gatk.broadinstitute.org/hc/en-us/articles/360035890671?id=6472) and [here](https://gatk.broadinstitute.org/hc/en-us/articles/360037226472-AddOrReplaceReadGroups-Picard-). More info about BGI seq can be found [here](https://en.wikipedia.org/wiki/DNA_nanoball_sequencing).

```bash
java -Xmx3g -jar $HPC_PICARD_DIR/picard.jar \
	AddOrReplaceReadGroups \
	I=${INPUT}/Tpr_${sample}_sorted_marked.bam \
	O=${INPUT}/Tpr_${sample}_sorted_marked_AddReadGroup.bam \
	RGID= CL100078369.L01.${sample}\
	RGLB= LIB-${sample}\
	RGPL=ILLUMINA \
	RGPU= CL100078369.L01.${sample}\
	RGSM=Tpr
```
**Platform is ILLUMINA, although it is actually BGISEQ**. Otherwise, `ERROR:INVALID_PLATFORM_VALUE`.

Based on the information from BGI, each paired fastq files has a unique barcode, and therefore, represents a uniqe library.
```
CL100078369_L01_533_1.fq
CL100078369_L01_533_2.fq
CL100078369_L01_534_1.fq
CL100078369_L01_534_2.fq
...
CL100078369_L01_548_1.fq
CL100078369_L01_548_2.fq
```
`CL100078369_L01_533_1.fq` and `CL100078369_L01_533_2.fq` are from a single library with barcode 533. `CL100078369_L01_534_1.fq` and `CL100078369_L01_534_2.fq` are from another library with barcode 534. All libraries were sequenced on Lane 1 of flow cell CL100078369.

### 2.7 Merge all Bam files (with read group name) and validate the combined Bam file
Script `Merge_bam.V2.sh` was used. The mergered bam file is named as `Tpr_combined_sorted_marked_AddReadGroup.bam`. The validation output is shown below:
```
No errors found
```
## 3. Call variants
Script `HaplotypeCaller_V1.sh` was used. For more information of the method: [here](https://gatk.broadinstitute.org/hc/en-us/articles/360035535932-Germline-short-variant-discovery-SNPs-Indels-) and [here](https://gatk.broadinstitute.org/hc/en-us/articles/360037225632-HaplotypeCaller). The reference fasta file need to be indexed and a dictionary file need to be created, for more [info](https://gatk.broadinstitute.org/hc/en-us/articles/360035531652-FASTA-Reference-genome-format).

**However**, because of the long running time when analyzing all scaffolds at once (e.g. after running for 5 days with 4 CPUs and 60 gb memory, scaffold 861 is being processed), we split the scaffolds and running different jobs at the same time. For more info of the method, please find it [here](https://gatk.broadinstitute.org/hc/en-us/articles/360035531852-Intervals-and-interval-lists).

Based on the size distribution (inside the image folder) of the scaffolds, 20 sublist (`split_scaffolds.xlsx`) was gernerated. Script `HaplotypeCaller_V5.sh` was used.
