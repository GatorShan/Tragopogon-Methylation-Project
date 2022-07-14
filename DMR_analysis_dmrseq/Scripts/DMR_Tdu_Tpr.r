library(dmrseq)
setwd("/orange/soltis/shan158538/Methylation_output/bismark_methylation_extractor")

# Load the object containing the methylation data from Tdu, Tpr, Tms
load("bismarkBSseq.filtered")
bismarkBSseq.filtered
pData(bismarkBSseq.filtered)

# Filter the object, and only contain samples from Tdu and Tpr
loci.idx <- which(DelayedMatrixStats::rowSums2(getCoverage(bismarkBSseq.filtered, type="Cov")==0) == 0)
sample.idx <- which(pData(bismarkBSseq.filtered)$species %in% c("Tdu", "Tpr"))
# When filerting the object, seems need both loci.idx and sample.idx; although the loci.idx has been applied
bismarkBSseq.filtered.Tdu_Tpr <- bismarkBSseq.filtered[loci.idx, sample.idx]
bismarkBSseq.filtered.Tdu_Tpr
pData(bismarkBSseq.filtered.Tdu_Tpr)

# Differentially mehtylated regions
testCovariate <- "species"
regions_Tdu_Tpr <- dmrseq(bs=bismarkBSseq.filtered.Tdu_Tpr,
                  cutoff = 0.05,
                  testCovariate=testCovariate)

# Save the DMR results
save(regions_Tdu_Tpr, file = "regions_Tdu_Tpr")
