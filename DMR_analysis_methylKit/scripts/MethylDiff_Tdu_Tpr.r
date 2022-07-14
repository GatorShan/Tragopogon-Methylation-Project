library(methylKit)
setwd("/orange/soltis/shan158538/Methylation_output/bismark_deduplicate")

file.list = list("Tdu_1_CpG.txt",
                 "Tdu_2_CpG.txt",
                 "Tpr_1_CpG.txt",
                 "Tpr_2_CpG.txt")
                 
sample.list = list("Tdu_1", "Tdu_2", "Tpr_1", "Tpr_2")

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
