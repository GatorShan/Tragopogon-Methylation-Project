# methyKit Analysis
## 1. Introduction
**methylKit** is an R package for DNA methylation analysis and annotation from high-throughput bisulfite sequencing. Following instructions from its  [github page](https://github.com/al2na/methylKit) and [online manual](https://bioconductor.org/packages/release/bioc/vignettes/methylKit/inst/doc/methylKit.html#23_Reading_the_methylation_calls_from_sorted_Bismark_alignments).

Scripts:
```r
myDiff25p.hyper_1 = getMethylDiff(myDiff_1,difference=25,qvalue=0.01,type="hyper", save.db = TRUE)
myDiff25p.hypo_1 = getMethylDiff(myDiff_1,difference=25,qvalue=0.01,type="hypo", save.db = TRUE)
myDiff25p_1 = getMethylDiff(myDiff_1,difference=25,qvalue=0.01, save.db = TRUE)
```


"difference: cutoff for ABSOLUTE VALUE of methylation percentage change between test and control (default:25)"; <span style="background-color: #FFFF00">THIS NUMBER SHOULD BE CHANGED FOR DIFFERENT CONTEXT!!!</span>

"type: one of the "hyper","hypo" or "all" strings. Specifies what type of differentially menthylated bases/regions should be returned. For retrieving Hyper-methylated regions/bases type="hyper", for hypo-methylated type="hypo" (default:"all")"


## 2. Reading the methylation calls from sorted Bismark alignments
### 2.1 Sort and index bam files
Script `Bismark_bam_formatting_V1.sh` was used. Input files are deduplicated.
### 2.2 Read Bismark alignment
Script `processBismarkAln_CG_V1.r` and `processBismarkAln_CG_V1.sh` were used.

|Sample.id|Original file|
|--|--|
|Tdu_1|HMCWKCCXY_s8_1_4981-LF_17_SL334590_val_1_bismark_bt2_pe_nameSorted.deduplicated_PosSorted.bam|
|Tdu_2|S1_cat_R1_val_1_bismark_bt2_pe_nameSorted.deduplicated_PosSorted.bam|
|Tpr_1|S2_cat_R1_val_1_bismark_bt2_pe_nameSorted.deduplicated_PosSorted.bam|
|Tpr_2|S3_cat_R1_val_1_bismark_bt2_pe_nameSorted.deduplicated_PosSorted.bam|
|Tms_1|S4_cat_R1_val_1_bismark_bt2_pe_nameSorted.deduplicated_PosSorted.bam|
|Tms_2|S5_cat_R1_val_1_bismark_bt2_pe_nameSorted.deduplicated_PosSorted.bam|

Example saved/output files:
```
Tdu_1_CpG.txt
Tdu_1_CpG_conversionStats.txt
Tdu_2_CpG.txt
Tdu_2_CpG_conversionStats.txt
...
```
### 2.3 Descriptive statistics on samples
Percent methylation and read coverage information were calculated using scripts `methylKit_DescriptiveStatistics.r` and `methylKit_DescriptiveStatistics.sh`. Results could be found in files `methylKit_DescriptiveStatistics.pdf` and `methylKit_DescriptiveStatistics_17957860.out`.

## 3. Comparative analysis
### 3.1 DMR between Tdu and Tpr (CpG context)
To identify DMRs between Tdu and Tpr, scripts `MethylDiff_Tdu_Tpr.r` and `MethylDiff_Tdu_Tpr.sh` were used.
