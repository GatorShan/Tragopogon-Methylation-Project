library(methylKit)
setwd("/blue/soltis/shan158538/Methylation/OutPut/Overlap_bases/CpG_mincov_3")

file.list = list("Tdu_1_cov3_CpG_overlap.txt",
                 "Tdu_2_cov3_CpG_overlap.txt",
                 "Tpr_1_cov3_CpG_overlap.txt",
                 "Tpr_2_cov3_CpG_overlap.txt")
                 
sample.list = list("Tdu_1_cov3_overlap_v4", "Tdu_2_cov3_overlap_v4", "Tpr_1_cov3_overlap_v4", "Tpr_2_cov3_overlap_v4")

myobj_lowCov = methRead(location = file.list,
                 sample.id = sample.list,
                 assembly="Tdu_ref_genome",
                 treatment=c(1,1,0,0),
                 mincov = 3,
                 context="CpG")

tiles = tileMethylCounts(myobj_lowCov,win.size=300,step.size=300,cov.bases=10)

meth=unite(tiles, destrand=TRUE)

myDiff=calculateDiffMeth(meth)
myDiff25p.hyper=getMethylDiff(myDiff,difference=35,qvalue=0.01,type="hyper", save.db = TRUE)
myDiff25p.hypo=getMethylDiff(myDiff,difference=35,qvalue=0.01,type="hypo", save.db = TRUE)
myDiff25p=getMethylDiff(myDiff,difference=35,qvalue=0.01, save.db = TRUE)
