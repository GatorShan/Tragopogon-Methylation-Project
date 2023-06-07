library(methylKit)
setwd("/blue/soltis/shan158538/Methylation/OutPut/Overlap_bases/CHH_mincov_3")

file.list = list("Tms_1_du_cov3_CHH_overlap.txt",
                 "Tms_2_du_cov3_CHH_overlap.txt",
                 "Tms_1_pr_cov3_CHH_overlap.txt",
                 "Tms_2_pr_cov3_CHH_overlap.txt")
                 
sample.list = list("Tms_1_du_cov3_CHH_overlap_v4", "Tms_2_du_cov3_CHH_overlap_v4", "Tms_1_pr_cov3_CHH_overlap_v4", "Tms_2_pr_cov3_CHH_overlap_v4")

myobj_lowCov = methRead(location = file.list,
                 sample.id = sample.list,
                 assembly="Tdu_ref_genome",
                 treatment=c(1,1,0,0),
                 mincov = 3,
                 context="CHH")

tiles = tileMethylCounts(myobj_lowCov,win.size=300,step.size=300,cov.bases=10)

# set destrand to FALSE for CHH
meth=unite(tiles, destrand=FALSE)

# cutoff for CHH is 10
myDiff=calculateDiffMeth(meth)
myDiff25p.hyper=getMethylDiff(myDiff,difference=10,qvalue=0.01,type="hyper", save.db = TRUE)
myDiff25p.hypo=getMethylDiff(myDiff,difference=10,qvalue=0.01,type="hypo", save.db = TRUE)
myDiff25p=getMethylDiff(myDiff,difference=10,qvalue=0.01, save.db = TRUE)
