library(methylKit)
setwd("/blue/soltis/shan158538/Methylation/OutPut/Overlap_bases/CpG_mincov_5")

file.list = list("Tms_1_du_cov5_CpG_overlap.txt",
                 "Tms_2_du_cov5_CpG_overlap.txt",
                 "Tms_1_pr_cov5_CpG_overlap.txt",
                 "Tms_2_pr_cov5_CpG_overlap.txt")
                 
sample.list = list("Tms_1_du_cov5_overlap", "Tms_2_du_cov5_overlap", "Tms_1_pr_cov5_overlap", "Tms_2_pr_cov5_overlap")

myobj_lowCov = methRead(location = file.list,
                 sample.id = sample.list,
                 assembly="Tdu_ref_genome",
                 treatment=c(1,1,0,0),
                 context="CpG")

tiles = tileMethylCounts(myobj_lowCov,win.size=300,step.size=300,cov.bases = 1)

meth=unite(tiles, destrand=TRUE)

myDiff=calculateDiffMeth(meth)
myDiff25p.hyper=getMethylDiff(myDiff,difference=25,qvalue=0.01,type="hyper", save.db = TRUE)
myDiff25p.hypo=getMethylDiff(myDiff,difference=25,qvalue=0.01,type="hypo", save.db = TRUE)
myDiff25p=getMethylDiff(myDiff,difference=25,qvalue=0.01, save.db = TRUE)