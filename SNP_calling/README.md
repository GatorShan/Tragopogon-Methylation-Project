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
Script `Validate_sam_file.V1.sh` was used to validate the combined BAM file. More information about BAM file validation can be found [here](https://gatk.broadinstitute.org/hc/en-us/articles/360036854731-ValidateSamFile-Picard-) and [here](https://gatk.broadinstitute.org/hc/en-us/articles/360035891231-Errors-in-SAM-or-BAM-files-can-be-diagnosed-with-ValidateSamFile)

The summary message is shown below:
```
## HISTOGRAM	java.lang.String
Error Type	Count
ERROR:MISSING_READ_GROUP	1
WARNING:RECORD_MISSING_READ_GROUP	1490169017
```
