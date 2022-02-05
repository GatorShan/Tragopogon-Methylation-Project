# Identify Differentially Methylated Regions (DMRs)
## 1. Introduction
1. The R package [dmrseq](https://github.com/kdkorthauer/dmrseq/blob/master/vignettes/dmrseq.Rmd) was used to identify DMRs.
2. Input files are from Bismark.
* example: `S1_cat_R1_val_1_bismark_bt2_pe_nameSorted.deduplicated.bismark.cov` from folder `bismark_methylation_extractor`
* only **CpG** methylation results are included in the above file

## 2. Data analysis
### 2.1 Load input files
```r
> library(dmrseq)
> setwd("/orange/soltis/shan158538/Methylation_output/bismark_methylation_extractor")
> bismarkBSseq <- read.bismark(files = c("DES1/HMCWKCCXY_s8_1_4981-LF_17_SL334590_val_1_bismark_bt2_pe_nameSorted.deduplicated.bismark.cov",
+                                       "S1/S1_cat_R1_val_1_bismark_bt2_pe_nameSorted.deduplicated.bismark.cov",
+                                       "S2/S2_cat_R1_val_1_bismark_bt2_pe_nameSorted.deduplicated.bismark.cov",
+                                       "S3/S3_cat_R1_val_1_bismark_bt2_pe_nameSorted.deduplicated.bismark.cov",
+                                       "S4/S4_cat_R1_val_1_bismark_bt2_pe_nameSorted.deduplicated.bismark.cov",
+                                       "S5/S5_cat_R1_val_1_bismark_bt2_pe_nameSorted.deduplicated.bismark.cov"
+                                       ),
+                              rmZeroCov = TRUE,
+                              strandCollapse = FALSE,
+                              verbose = TRUE)
> pData(bismarkBSseq)$species <- c("Tdu", "Tdu", "Tpr", "Tpr","Tms", "Tms")
> pData(bismarkBSseq)$individual <- c("3060-1-4", "3040-6-2", "3058-1-2", "3058-4-10", "3059-7-7", "3059-21-5")
> loci.idx <- which(DelayedMatrixStats::rowSums2(getCoverage(bismarkBSseq, type="Cov")==0) == 0)
> bismarkBSseq.filtered <- bismarkBSseq[loci.idx]
> save(bismarkBSseq.filtered, file = "bismarkBSseq.filtered")
```
