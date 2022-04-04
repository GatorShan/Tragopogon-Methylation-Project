library(methylKit)


file.list = list("/orange/soltis/shan158538/Methylation_output/bismark_deduplicate/Tdu_1_CpG.txt",
				 "/orange/soltis/shan158538/Methylation_output/bismark_deduplicate/Tdu_2_CpG.txt",
				 "/orange/soltis/shan158538/Methylation_output/bismark_deduplicate/Tpr_1_CpG.txt",
                 "/orange/soltis/shan158538/Methylation_output/bismark_deduplicate/Tpr_2_CpG.txt")
                 
sample.list = list("Tdu_1", "Tdu_2", "Tpr_1", "Tpr_2")


# Annalyzing CpG context between Tdu and Tpr; permutation 0 0 1 1
setwd("/blue/soltis/shan158538/Methylation/OutPut/DMR_methylkit/permutation_CpG_Tdu_Tpr/0011")

myobj_lowCov_1 = methRead(location = file.list,
                 sample.id = sample.list,
                 assembly="Tdu_ref_genome",
                 treatment=c(0,0,1,1),
                 context="CpG",
                 mincov= 3)

tiles_1 = tileMethylCounts(myobj_lowCov_1,win.size=1000,step.size=1000,cov.bases = 10)
meth_1 = unite(tiles_1, destrand=FALSE)
myDiff_1 = calculateDiffMeth(meth_1)
myDiff25p.hyper_1 = getMethylDiff(myDiff_1,difference=25,qvalue=0.01,type="hyper", save.db = TRUE)
myDiff25p.hypo_1 = getMethylDiff(myDiff_1,difference=25,qvalue=0.01,type="hypo", save.db = TRUE)
myDiff25p_1 = getMethylDiff(myDiff_1,difference=25,qvalue=0.01, save.db = TRUE)


# Annalyzing CpG context between Tdu and Tpr; permutation 0 1 0 1
setwd("/blue/soltis/shan158538/Methylation/OutPut/DMR_methylkit/permutation_CpG_Tdu_Tpr/0101")

myobj_lowCov_2 = methRead(location = file.list,
                 sample.id = sample.list,
                 assembly="Tdu_ref_genome",
                 treatment=c(0,1,0,1),
                 context="CpG",
                 mincov= 3)

tiles_2 = tileMethylCounts(myobj_lowCov_2,win.size=1000,step.size=1000,cov.bases = 10)
meth_2 = unite(tiles_2, destrand=FALSE)
myDiff_2 = calculateDiffMeth(meth_2)
myDiff25p.hyper_2 = getMethylDiff(myDiff_2,difference=25,qvalue=0.01,type="hyper", save.db = TRUE)
myDiff25p.hypo_2 = getMethylDiff(myDiff_2,difference=25,qvalue=0.01,type="hypo", save.db = TRUE)
myDiff25p_2 = getMethylDiff(myDiff_2,difference=25,qvalue=0.01, save.db = TRUE)                 
                 

# Annalyzing CpG context between Tdu and Tpr; permutation 0 1 1 0
setwd("/blue/soltis/shan158538/Methylation/OutPut/DMR_methylkit/permutation_CpG_Tdu_Tpr/0110")

myobj_lowCov_3 = methRead(location = file.list,
                 sample.id = sample.list,
                 assembly="Tdu_ref_genome",
                 treatment=c(0,1,1,0),
                 context="CpG",
                 mincov= 3)
                 
tiles_3 = tileMethylCounts(myobj_lowCov_3,win.size=1000,step.size=1000,cov.bases = 10)
meth_3 = unite(tiles_3, destrand=FALSE)
myDiff_3 = calculateDiffMeth(meth_3)
myDiff25p.hyper_3 = getMethylDiff(myDiff_3,difference=25,qvalue=0.01,type="hyper", save.db = TRUE)
myDiff25p.hypo_3 = getMethylDiff(myDiff_3,difference=25,qvalue=0.01,type="hypo", save.db = TRUE)
myDiff25p_3 = getMethylDiff(myDiff_3,difference=25,qvalue=0.01, save.db = TRUE)        


# Annalyzing CpG context between Tdu and Tpr; permutation 1 0 0 1
setwd("/blue/soltis/shan158538/Methylation/OutPut/DMR_methylkit/permutation_CpG_Tdu_Tpr/1001")

myobj_lowCov_4 = methRead(location = file.list,
                 sample.id = sample.list,
                 assembly="Tdu_ref_genome",
                 treatment=c(1,0,0,1),
                 context="CpG",
                 mincov= 3)
                 
tiles_4 = tileMethylCounts(myobj_lowCov_4,win.size=1000,step.size=1000,cov.bases = 10)
meth_4 = unite(tiles_4, destrand=FALSE)
myDiff_4 = calculateDiffMeth(meth_4)
myDiff25p.hyper_4 = getMethylDiff(myDiff_4,difference=25,qvalue=0.01,type="hyper", save.db = TRUE)
myDiff25p.hypo_4 = getMethylDiff(myDiff_4,difference=25,qvalue=0.01,type="hypo", save.db = TRUE)
myDiff25p_4 = getMethylDiff(myDiff_4,difference=25,qvalue=0.01, save.db = TRUE)


# Annalyzing CpG context between Tdu and Tpr; permutation 1 0 1 0
setwd("/blue/soltis/shan158538/Methylation/OutPut/DMR_methylkit/permutation_CpG_Tdu_Tpr/1010")

myobj_lowCov_5 = methRead(location = file.list,
                 sample.id = sample.list,
                 assembly="Tdu_ref_genome",
                 treatment=c(1,0,1,0),
                 context="CpG",
                 mincov= 3)

tiles_5 = tileMethylCounts(myobj_lowCov_5,win.size=1000,step.size=1000,cov.bases = 10)
meth_5 = unite(tiles_5, destrand=FALSE)
myDiff_5 = calculateDiffMeth(meth_5)
myDiff25p.hyper_5 = getMethylDiff(myDiff_5,difference=25,qvalue=0.01,type="hyper", save.db = TRUE)
myDiff25p.hypo_5 = getMethylDiff(myDiff_5,difference=25,qvalue=0.01,type="hypo", save.db = TRUE)
myDiff25p_5 = getMethylDiff(myDiff_5,difference=25,qvalue=0.01, save.db = TRUE)