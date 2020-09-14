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
Script `Picard_MarkDuplicates_V3.sh` was used to sort the bam file and mark the duplicates.

Tips:
  - Requrest a few more memory (25 Gb) than what is allocated to the java script (-Xmx 20 Gb)
  - Set up a new TMP directory at a local directory (error with the default TMP directory)

![Image_1](https://github.com/GatorShan/Tragopogon-Methylation-Project/blob/master/SNP_calling/images/Bam_sorted_marked_file_multiqc_report_samtools_alignment_plot.png)

![Image_2](https://github.com/GatorShan/Tragopogon-Methylation-Project/blob/master/SNP_calling/images/Bam_sorted_marked_file_multiqc_report_gerneral_statistics.png)

![Image_3](https://github.com/GatorShan/Tragopogon-Methylation-Project/blob/master/SNP_calling/images/Bam_sorted_marked_file_multiqc_report_alignment_matrics.png)
