library(methylKit)
setwd("/orange/soltis/shan158538/Methylation_output/bismark_deduplicate")

# Reading the methylation calls from sorted Bismark alignments
my.methRaw=processBismarkAln(location=list("HMCWKCCXY_s8_1_4981-LF_17_SL334590_val_1_bismark_bt2_pe_nameSorted.deduplicated_PosSorted.bam",
                                            "S1_cat_R1_val_1_bismark_bt2_pe_nameSorted.deduplicated_PosSorted.bam",
                                            "S2_cat_R1_val_1_bismark_bt2_pe_nameSorted.deduplicated_PosSorted.bam",
                                            "S3_cat_R1_val_1_bismark_bt2_pe_nameSorted.deduplicated_PosSorted.bam",
                                            "S4_cat_R1_val_1_bismark_bt2_pe_nameSorted.deduplicated_PosSorted.bam",
                                            "S5_cat_R1_val_1_bismark_bt2_pe_nameSorted.deduplicated_PosSorted.bam"),
                              sample.id=list("Tdu_1_cov3", "Tdu_2_cov3", "Tpr_1_cov3", "Tpr_2_cov3", "Tms_1_cov3", "Tms_2_cov3"),
                              assembly="Tdu_ref_genome",
                              read.context="CpG",
                              mincov = 3,
                              save.folder=getwd(),
                              treatment=c(2,2,1,1,0,0))
