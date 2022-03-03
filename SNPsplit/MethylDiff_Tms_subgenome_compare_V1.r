library(methylKit)
setwd("/blue/soltis/shan158538/Methylation/OutPut/SNPsplit_output")

file.list = list("Tms_1_du_CpG.txt",
                 "Tms_2_du_CpG.txt",
                 "Tms_1_pr_CpG.txt",
                 "Tms_2_pr_CpG.txt")
                 
sample.list = list("Tms_1_du", "Tms_2_du", "Tms_1_pr", "Tms_2_pr")

myobj_lowCov = methRead(location = file.list,
                 sample.id = sample.list,
                 assembly="Tdu_ref_genome",
                 treatment=c(1,1,0,0),
                 context="CpG",
                 mincov= 3)

tiles = tileMethylCounts(myobj_lowCov,win.size=1000,step.size=1000,cov.bases = 10)

meth=unite(tiles, destrand=TRUE)

myDiff=calculateDiffMeth(meth)
myDiff25p.hyper=getMethylDiff(myDiff,difference=25,qvalue=0.01,type="hyper", save.db = TRUE)
myDiff25p.hypo=getMethylDiff(myDiff,difference=25,qvalue=0.01,type="hypo", save.db = TRUE)
myDiff25p=getMethylDiff(myDiff,difference=25,qvalue=0.01, save.db = TRUE)
diffMethPerChr(myDiff,plot=TRUE,qvalue.cutoff=0.01, meth.cutoff=25)