library(methylKit)
setwd("/orange/soltis/shan158538/Methylation_output/bismark_deduplicate")

file.list = list("Tdu_1_CpG.txt",
                 "Tdu_2_CpG.txt",
                 "Tpr_1_CpG.txt",
                 "Tpr_2_CpG.txt",
                 "Tms_1_CpG.txt",
                 "Tms_2_CpG.txt")
                 
sample.list = list("Tdu_1", "Tdu_2", "Tpr_1", "Tpr_2", "Tms_1", "Tms_2")

myobj = methRead(location = file.list,
                 sample.id = sample.list,
                 assembly="Tdu_ref_genome",
                 treatment=c(2,2,1,1,0,0))
                 
x <- c(1,2,3,4,5,6)
for (val in x) {
    print(paste("Processing sample:", val))
    getMethylationStats(myobj[[val]],plot=FALSE,both.strands=FALSE)
    getMethylationStats(myobj[[val]],plot=TRUE,both.strands=FALSE)
    getCoverageStats(myobj[[val]],plot=TRUE,both.strands=FALSE)
}