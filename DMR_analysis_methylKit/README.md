# DMR analysis using methylKit
## 1. Introduction
**methylKit** is an R package for DNA methylation analysis and annotation from high-throughput bisulfite sequencing. Following instructions from its  [github page](https://github.com/al2na/methylKit) and [online manual](https://bioconductor.org/packages/release/bioc/vignettes/methylKit/inst/doc/methylKit.html#23_Reading_the_methylation_calls_from_sorted_Bismark_alignments).

**This section is focusing on DMR analysis and comparison. The results of methylation calling using methylKit could be found [here](https://github.com/GatorShan/Tragopogon-Methylation-Project/blob/master/methylKit_analysis/README.md#2-reading-the-methylation-calls-from-sorted-bismark-alignments).**

Scripts:
```r
meth=unite(tiles, destrand=TRUE)
### By default, unite function produces bases/regions covered in all samples

myDiff=calculateDiffMeth(meth)
myDiff25p.hyper = getMethylDiff(myDiff,difference=25,qvalue=0.01,type="hyper", save.db = TRUE)
myDiff25p.hypo = getMethylDiff(myDiff,difference=25,qvalue=0.01,type="hypo", save.db = TRUE)
myDiff25p = getMethylDiff(myDiff,difference=25,qvalue=0.01, save.db = TRUE)
```

-- "difference: cutoff for ABSOLUTE VALUE of methylation percentage change between test and control (default:25)"; **THIS NUMBER SHOULD BE CHANGED FOR DIFFERENT CONTEXT!!!**

-- "type: one of the "hyper","hypo" or "all" strings. Specifies what type of differentially menthylated bases/regions should be returned. For retrieving Hyper-methylated regions/bases type="hyper", for hypo-methylated type="hypo" (default:"all")"

-- "Setting destrand=TRUE (the default is FALSE) will merge reads on both strands of a CpG dinucleotide. This provides better coverage, but only advised when looking at CpG methylation (for CpH methylation this will cause wrong results in subsequent analyses)."

## 2. Results
| Context | processBismarkAln mincov | using shared bases? | methRead mincov | Window size | Window step | Window cov.bases | Difference | DMR between parents | DMR between subgenomes | Overlapping DMR (both direaction) | Script |
|--|--|--|--|--|--|--|--|--|--|--|--|
| CpG | 3 | yes | 3 | 300 | 300 | 10 | 35% | 3,549 | 3,262 | 2,360 | `MethylDiff_Tdu-Tpr_CpG_mincov3_V4.r` and `MethylDiff_Tms_subgenome_compare_mincov3_V4.r` |
| CHG | 3 | yes | 3 | 300 | 300 | 10 | 10% |  |  |  | `MethylDiff_Tdu-Tpr_CHG_mincov3_V4.r` and `MethylDiff_Tms_subgenome_compare_CHG_mincov3_V4.r` |
| CHH | 3 | yes | 3 | 300 | 300 | 10 | 10% |  |  |  | `MethylDiff_Tdu-Tpr_CHH_mincov3_V4.r` and `MethylDiff_Tms_subgenome_compare_CHH_mincov3_V4.r` |

## 3. [Previous analyses](https://github.com/GatorShan/Tragopogon-Methylation-Project/tree/master/DMR_analysis_methylKit/Previous_analyses)
