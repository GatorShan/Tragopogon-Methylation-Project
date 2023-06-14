library(methylKit)
setwd("/blue/soltis/shan158538/Methylation/OutPut/Overlap_bases/CHG_mincov_3")

file.list = list("Tdu_1_cov3_CHG_overlap.txt",
                 "Tdu_2_cov3_CHG_overlap.txt",
                 "Tpr_1_cov3_CHG_overlap.txt",
                 "Tpr_2_cov3_CHG_overlap.txt")
                 
sample.list = list("Tdu_1_cov3_CHG_overlap_v4", "Tdu_2_cov3_CHG_overlap_v4", "Tpr_1_cov3_CHG_overlap_v4", "Tpr_2_cov3_CHG_overlap_v4")

myobj_lowCov = methRead(location = file.list,
                 sample.id = sample.list,
                 assembly="Tdu_ref_genome",
                 treatment=c(1,1,0,0),
                 mincov = 3,
                 context="CHG")

tiles = tileMethylCounts(myobj_lowCov,win.size=300,step.size=300,cov.bases=10)

# destrand set as FALSE for CHG context
meth=unite(tiles, destrand=FALSE)

# difference cutoff set as 25 for CHG context (Bob suggested this value on 06/14/2023)
myDiff=calculateDiffMeth(meth)
myDiff25p.hyper=getMethylDiff(myDiff,difference=25,qvalue=0.01,type="hyper", save.db = TRUE)
myDiff25p.hypo=getMethylDiff(myDiff,difference=25,qvalue=0.01,type="hypo", save.db = TRUE)
myDiff25p=getMethylDiff(myDiff,difference=25,qvalue=0.01, save.db = TRUE)
