library(methylKit)

# Annalyzing CHH context between Tdu and Tpr
setwd("/blue/soltis/shan158538/Methylation/OutPut/DMR_methylkit/CHH_Tdu_vs_Tpr")

file.list_1 = list("/orange/soltis/shan158538/Methylation_output/bismark_deduplicate/Tdu_1_CHH.txt",
				 "/orange/soltis/shan158538/Methylation_output/bismark_deduplicate/Tdu_2_CHH.txt",
				 "/orange/soltis/shan158538/Methylation_output/bismark_deduplicate/Tpr_1_CHH.txt",
                 "/orange/soltis/shan158538/Methylation_output/bismark_deduplicate/Tpr_2_CHH.txt")
                 
sample.list_1 = list("Tdu_1", "Tdu_2", "Tpr_1", "Tpr_2")

myobj_lowCov_1 = methRead(location = file.list_1,
                 sample.id = sample.list_1,
                 assembly="Tdu_ref_genome",
                 treatment=c(1,1,0,0),
                 context="CHH",
                 mincov= 3)

tiles_1 = tileMethylCounts(myobj_lowCov_1,win.size=1000,step.size=1000,cov.bases = 10)

meth_1 = unite(tiles_1, destrand=TRUE)

myDiff_1 = calculateDiffMeth(meth_1)
myDiff25p.hyper_1 = getMethylDiff(myDiff_1,difference=25,qvalue=0.01,type="hyper", save.db = TRUE)
myDiff25p.hypo_1 = getMethylDiff(myDiff_1,difference=25,qvalue=0.01,type="hypo", save.db = TRUE)
myDiff25p_1 = getMethylDiff(myDiff_1,difference=25,qvalue=0.01, save.db = TRUE)


# Annalyzing CHH context between Tms-du and Tms-pr
setwd("/blue/soltis/shan158538/Methylation/OutPut/DMR_methylkit/CHH_Tms-du_vs_Tms-pr")

file.list_2 = list("/blue/soltis/shan158538/Methylation/OutPut/SNPsplit_output/Tms_1_du_CHH.txt",
				 "/blue/soltis/shan158538/Methylation/OutPut/SNPsplit_output/Tms_2_du_CHH.txt",
				 "/blue/soltis/shan158538/Methylation/OutPut/SNPsplit_output/Tms_1_pr_CHH.txt",
                 "/blue/soltis/shan158538/Methylation/OutPut/SNPsplit_output/Tms_2_pr_CHH.txt")
                 
sample.list_2 = list("Tms_1_du", "Tms_2_du", "Tms_1_pr", "Tms_2_pr")

myobj_lowCov_2 = methRead(location = file.list_2,
                 sample.id = sample.list_2,
                 assembly="Tdu_ref_genome",
                 treatment=c(1,1,0,0),
                 context="CHH",
                 mincov= 3)
                 
tiles_2 = tileMethylCounts(myobj_lowCov_2,win.size=1000,step.size=1000,cov.bases = 10)

meth_2 = unite(tiles_2, destrand=TRUE)

myDiff_2 = calculateDiffMeth(meth_2)
myDiff25p.hyper_2 = getMethylDiff(myDiff_2,difference=25,qvalue=0.01,type="hyper", save.db = TRUE)
myDiff25p.hypo_2 = getMethylDiff(myDiff_2,difference=25,qvalue=0.01,type="hypo", save.db = TRUE)
myDiff25p_2 = getMethylDiff(myDiff_2,difference=25,qvalue=0.01, save.db = TRUE)                 
                 
                 
# Annalyzing CHH context between Tdu and Tms-du
setwd("/blue/soltis/shan158538/Methylation/OutPut/DMR_methylkit/CHH_Tdu_vs_Tms-du")

file.list_3 = list("/orange/soltis/shan158538/Methylation_output/bismark_deduplicate/Tdu_1_CHH.txt",
				 "/orange/soltis/shan158538/Methylation_output/bismark_deduplicate/Tdu_2_CHH.txt",
				 "/blue/soltis/shan158538/Methylation/OutPut/SNPsplit_output/Tms_1_du_CHH.txt",
				 "/blue/soltis/shan158538/Methylation/OutPut/SNPsplit_output/Tms_2_du_CHH.txt") 

sample.list_3 = list("Tdu_1", "Tdu_2", "Tms_1_du", "Tms_2_du")

myobj_lowCov_3 = methRead(location = file.list_3,
                 sample.id = sample.list_3,
                 assembly="Tdu_ref_genome",
                 treatment=c(1,1,0,0),
                 context="CHH",
                 mincov= 3)
                 
tiles_3 = tileMethylCounts(myobj_lowCov_3,win.size=1000,step.size=1000,cov.bases = 10)

meth_3 = unite(tiles_3, destrand=TRUE)

myDiff_3 = calculateDiffMeth(meth_3)
myDiff25p.hyper_3 = getMethylDiff(myDiff_3,difference=25,qvalue=0.01,type="hyper", save.db = TRUE)
myDiff25p.hypo_3 = getMethylDiff(myDiff_3,difference=25,qvalue=0.01,type="hypo", save.db = TRUE)
myDiff25p_3 = getMethylDiff(myDiff_3,difference=25,qvalue=0.01, save.db = TRUE)        


# Annalyzing CHH context between Tpr and Tms-pr
setwd("/blue/soltis/shan158538/Methylation/OutPut/DMR_methylkit/CHH_Tpr_vs_Tms-pr")

file.list_4 = list("/orange/soltis/shan158538/Methylation_output/bismark_deduplicate/Tpr_1_CHH.txt",
				 "/orange/soltis/shan158538/Methylation_output/bismark_deduplicate/Tpr_2_CHH.txt",
				 "/blue/soltis/shan158538/Methylation/OutPut/SNPsplit_output/Tms_1_pr_CHH.txt",
				 "/blue/soltis/shan158538/Methylation/OutPut/SNPsplit_output/Tms_2_pr_CHH.txt") 

sample.list_4 = list("Tpr_1", "Tpr_2", "Tms_1_pr", "Tms_2_pr")

myobj_lowCov_4 = methRead(location = file.list_4,
                 sample.id = sample.list_4,
                 assembly="Tdu_ref_genome",
                 treatment=c(1,1,0,0),
                 context="CHH",
                 mincov= 3)
                 
tiles_4 = tileMethylCounts(myobj_lowCov_4,win.size=1000,step.size=1000,cov.bases = 10)

meth_4 = unite(tiles_4, destrand=TRUE)

myDiff_4 = calculateDiffMeth(meth_4)
myDiff25p.hyper_4 = getMethylDiff(myDiff_4,difference=25,qvalue=0.01,type="hyper", save.db = TRUE)
myDiff25p.hypo_4 = getMethylDiff(myDiff_4,difference=25,qvalue=0.01,type="hypo", save.db = TRUE)
myDiff25p_4 = getMethylDiff(myDiff_4,difference=25,qvalue=0.01, save.db = TRUE)

