# DMR analysis using methylKit
## 1. Introduction
**methylKit** is an R package for DNA methylation analysis and annotation from high-throughput bisulfite sequencing. Following instructions from its  [github page](https://github.com/al2na/methylKit) and [online manual](https://bioconductor.org/packages/release/bioc/vignettes/methylKit/inst/doc/methylKit.html#23_Reading_the_methylation_calls_from_sorted_Bismark_alignments).

This section is focusing on DMR analysis and comparison.

Scripts:
```r
myDiff25p.hyper_1 = getMethylDiff(myDiff_1,difference=25,qvalue=0.01,type="hyper", save.db = TRUE)
myDiff25p.hypo_1 = getMethylDiff(myDiff_1,difference=25,qvalue=0.01,type="hypo", save.db = TRUE)
myDiff25p_1 = getMethylDiff(myDiff_1,difference=25,qvalue=0.01, save.db = TRUE)
```


"difference: cutoff for ABSOLUTE VALUE of methylation percentage change between test and control (default:25)"; **THIS NUMBER SHOULD BE CHANGED FOR DIFFERENT CONTEXT!!!**

"type: one of the "hyper","hypo" or "all" strings. Specifies what type of differentially menthylated bases/regions should be returned. For retrieving Hyper-methylated regions/bases type="hyper", for hypo-methylated type="hypo" (default:"all")"
