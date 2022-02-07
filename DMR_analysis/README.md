# Identify Differentially Methylated Regions (DMRs)
## 1. Introduction
1. The R package [dmrseq](https://github.com/kdkorthauer/dmrseq/blob/master/vignettes/dmrseq.Rmd) was used to identify DMRs.
2. Input files are from Bismark.
* example: `S1_cat_R1_val_1_bismark_bt2_pe_nameSorted.deduplicated.bismark.cov` from folder `bismark_methylation_extractor`
* only **CpG** methylation results are included in the above file

## 2. Data analysis
### 2.1 Load input files
The following codes were executed at [RStudio Server](https://help.rc.ufl.edu/doc/RStudio_Server) from UFRC. **The memory requsted is 32 Gb** -- if using small amount of memory, an error will pop up.

```r
# Load the dmrseq package
> library(dmrseq)

# set the dirctory of the input files
> setwd("/orange/soltis/shan158538/Methylation_output/bismark_methylation_extractor")

# Bismark input loading
# Tdu: DES1 and S1; Tpr: S2 and S3; Tms: S4 and S5
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

# Add sample metadata; pData contains the covariate of interest
> pData(bismarkBSseq)$species <- c("Tdu", "Tdu", "Tpr", "Tpr","Tms", "Tms")
> pData(bismarkBSseq)$individual <- c("3060-1-4", "3040-6-2", "3058-1-2", "3058-4-10", "3059-7-7", "3059-21-5")

# Filtering bismarkBSseq object; remove all CpGs that have no coverage in at least one sample
> loci.idx <- which(DelayedMatrixStats::rowSums2(getCoverage(bismarkBSseq, type="Cov")==0) == 0)
> bismarkBSseq.filtered <- bismarkBSseq[loci.idx]

# Save the filtered object
> save(bismarkBSseq.filtered, file = "bismarkBSseq.filtered")
```
### 2.2 DMR analysis between Tdu and Tpr
Scripts `DMR_Tdu_Tpr.r` and `DMR_Tdu_Tpr.sh` were used. It took 1 d and 10 hours to finish running the program. 60 gb memeory was requested. Downstream analyses were performed at RStudio server.
```r
# Load the DMR results; DMR between Tdu and Tpr
load("regions_Tdu_Tpr")
regions_Tdu_Tpr

# How many retions were significantly at the DFR cutoff at 0.05?
sum(regions_Tdu_Tpr$qval < 0.05)

# Select just the regions below FDR 0.05 and place in a new data.frame
sigRegions_Tdu_Tpr <- regions_Tdu_Tpr[regions_Tdu_Tpr$qval < 0.05,]
sigRegions_Tdu_Tpr

# Hypo- or Hyper- methylation; determin the proportion of regions with hyper- or hypo-metylation
# In this case, condition: Tpr vs Tdu; >0 means hyper in Tpr; <0 means hyper is Tdu
# Hyper in Tpr
sum(sigRegions_Tdu_Tpr$stat > 0)
sum(sigRegions_Tdu_Tpr$stat > 0) / length(sigRegions_Tdu_Tpr)
# Hyper in Tdu
sum(sigRegions_Tdu_Tpr$stat < 0)
sum(sigRegions_Tdu_Tpr$stat < 0) / length(sigRegions_Tdu_Tpr)

# Plot DMRs; here plots the top 4 regions
plotDMRs(bismarkBSseq.filtered.Tdu_Tpr, regions=sigRegions_Tdu_Tpr[1:4,], testCovariate="species")

# Exporting results to CSV files
write.csv(as.data.frame(sigRegions_Tdu_Tpr), 
          file="sigDMR_Tdu_Tpr_results.csv")
```
Some results
| Category | Number |
| - | - |
| Sig. regions (FDR < 0.05) | 18,298 |
| Hyper methylated regions in Tdu | 13,419 (73.3%) |
| Hyper methylated regions in Tpr | 4,879 (26.7%) |





