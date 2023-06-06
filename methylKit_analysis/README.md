# methylKit Analysis
## 1. Introduction
**methylKit** is an R package for DNA methylation analysis and annotation from high-throughput bisulfite sequencing. Following instructions from its  [github page](https://github.com/al2na/methylKit) and [online manual](https://bioconductor.org/packages/release/bioc/vignettes/methylKit/inst/doc/methylKit.html#23_Reading_the_methylation_calls_from_sorted_Bismark_alignments).

**Using methylKit, this section is focusing on reading methylation calls from sorted Bismark alignment**. DMR analysis results could be found **[here](https://github.com/GatorShan/Tragopogon-Methylation-Project/tree/master/DMR_analysis_methylKit#dmr-analysis-using-methylkit)**.

## 2. Reading the methylation calls from sorted Bismark alignments

### IMPORTANT RESULTS: Read Bismark alignment for T. dubius, T. pratensis, T. miscellus, and two subgenomes; mincov = 3
|Sample.id|Original file|
|--|--|
|Tdu_1|HMCWKCCXY_s8_1_4981-LF_17_SL334590_val_1_bismark_bt2_pe_nameSorted.deduplicated_PosSorted.bam|
|Tdu_2|S1_cat_R1_val_1_bismark_bt2_pe_nameSorted.deduplicated_PosSorted.bam|
|Tpr_1|S2_cat_R1_val_1_bismark_bt2_pe_nameSorted.deduplicated_PosSorted.bam|
|Tpr_2|S3_cat_R1_val_1_bismark_bt2_pe_nameSorted.deduplicated_PosSorted.bam|
|Tms_1|S4_cat_R1_val_1_bismark_bt2_pe_nameSorted.deduplicated_PosSorted.bam|
|Tms_2|S5_cat_R1_val_1_bismark_bt2_pe_nameSorted.deduplicated_PosSorted.bam|
| Tms_1_du | S4_cat_R1_val_1_bismark_bt2_pe_nameSorted.deduplicated.genome1_PosSorted.bam |
| Tms_2_du | S5_cat_R1_val_1_bismark_bt2_pe_nameSorted.deduplicated.genome1_PosSorted.bam |
| Tms_1_pr | S4_cat_R1_val_1_bismark_bt2_pe_nameSorted.deduplicated.genome2_PosSorted.bam |
| Tms_2_pr | S5_cat_R1_val_1_bismark_bt2_pe_nameSorted.deduplicated.genome2_PosSorted.bam |
#### 2.1 CpG methylation
For T. dubius, T. pratensis, and T. miscellus, scripts `processBismarkAln_CG_minCOV-3_V1.r` and `processBismarkAln_CG_minCOV-3_V1.sh` were used. For the subgneomes of T. miscellus, scripts `processBismarkAln_subgenome_compare_CG_minCOV-3_V1.r` and `processBismarkAln_subgenome_compare_CG_minCOV-3_V1.sh` were used.

Example saved/output files:
```
Tdu_1_cov3_CpG.txt
Tdu_2_cov3_CpG.txt
Tms_1_du_cov3_CpG.txt
Tms_1_pr_cov3_CpG.txt
```

**We then identify overlapping bases cross diploid parents and subgenomes. The output from this step is used for downstream DMR comparing analysis**. Scirpt `Overlap_bases.V2.py` and `Overlap_bases.V2-1.sh` were used. The output is located at folder `/blue/soltis/shan158538/Methylation/OutPut/Overlap_bases/CpG_mincov_3`. **In total, there are 2,048,198 shared bases.**

Example outputs:
```
Tdu_1_cov3_CpG_overlap.txt
Tdu_2_cov3_CpG_overlap.txt
Tms_1_du_cov3_CpG_overlap.txt
```

#### 2.2 CHG methylation
Scripts `processBismarkAln_CHG_minCOV-3_V1.r` and `processBismarkAln_CHG_minCOV-3_V1.sh` were used. For the subgneomes of T. miscellus, scripts `processBismarkAln_subgenome_compare_CHG_minCOV-3_V1.r` and `processBismarkAln_subgenome_compare_CHG_minCOV-3_V1.sh` were used. To identify overlapped bases, scirpt `Overlap_bases.V2-1_CHG.sh` was used. The output is located at folder `/blue/soltis/shan158538/Methylation/OutPut/Overlap_bases/CHG_mincov_3`. **In total, there are 1,899,609 shared bases.**


#### 2.3 CHH methylation
Scripts `processBismarkAln_CHH_minCOV-3_V1.r` and `processBismarkAln_CHH_minCOV-3_V1.sh` were used. For the subgneomes of T. miscellus, scripts `processBismarkAln_subgenome_compare_CHH_minCOV-3_V1.r` and `processBismarkAln_subgenome_compare_CHH_minCOV-3_V1.sh` were used. To identify overlapped bases, scirpt `Overlap_bases.V2-1_CHH.sh` was used. The output is located at folder `/blue/soltis/shan158538/Methylation/OutPut/Overlap_bases/CHH_mincov_3`. **In total, there are 11,441,306 shared bases.**

#### 2.4 Additional notes
The coverage and methylation rate showed some difference between the **methylKit** and the **bismark methylation profile extration** pipelines. Both pipelines used the same deduplicated bam files from the bismark alignment (nameSorted.deduplicated and nameSorted.deduplicated_PosSorted for bismark extractor and methylKit, respectively). I believe the difference is resulted from the extra filters from the [bismark extractor](https://www.bioinformatics.babraham.ac.uk/projects/bismark/Bismark_User_Guide.pdf): `--no-overlap` and `--ignore_r2 2`. Results from [methylKit](https://www.rdocumentation.org/packages/methylKit/versions/0.99.2/topics/processBismarkAln) were used to identify DFR, and results from bismark extractor were used for constructing metaplot etc. Therefore, I believe the differences won't introduce biases to the results.


## 3. [Unused previous analyses](https://github.com/GatorShan/Tragopogon-Methylation-Project/tree/master/methylKit_analysis/Unused_previous_analyses)
