library(methylKit)
setwd("/blue/soltis/shan158538/Methylation/OutPut/SNPsplit_output")

# Reading the methylation calls from sorted Bismark alignments
my.methRaw=processBismarkAln(location=list("S4_cat_R1_val_1_bismark_bt2_pe_nameSorted.deduplicated.genome1_PosSorted.bam",
                                        	"S5_cat_R1_val_1_bismark_bt2_pe_nameSorted.deduplicated.genome1_PosSorted.bam",
                                            "S4_cat_R1_val_1_bismark_bt2_pe_nameSorted.deduplicated.genome2_PosSorted.bam",
                                            "S5_cat_R1_val_1_bismark_bt2_pe_nameSorted.deduplicated.genome2_PosSorted.bam"),
                              sample.id=list("Tms_1_du_cov3", "Tms_2_du_cov3", "Tms_1_pr_cov3", "Tms_2_pr_cov3"),
                              assembly="Tdu_ref_genome", 
                              read.context="CpG",
                              mincov = 3,
                              save.folder=getwd(),
                              treatment=c(1,1,0,0))
