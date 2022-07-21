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


"difference: cutoff for ABSOLUTE VALUE of methylation percentage change between test and control (default:25)"; **THIS NUMBER SHOULD BE CHANGED FOR DIFFERENT CONTEXT!!!**

"type: one of the "hyper","hypo" or "all" strings. Specifies what type of differentially menthylated bases/regions should be returned. For retrieving Hyper-methylated regions/bases type="hyper", for hypo-methylated type="hypo" (default:"all")"

## 2. CpG methylation
### 2.1 DMR between Tdu and Tpr, and DMR between two subgenomes
**SUMMARY OF DIFFERENT SETTINGS**
| processBismarkAln mincov | using shared bases? | methRead mincov | Window size | Window step | Window cov.bases | Difference | DMR between parents | DMR between subgenomes | Overlapping DMR (both direaction) | Script |
|--|--|--|--|--|--|--|--|--|--|--|
| 10 | no | 3 | 1000 | 1000 | 10 | 25% | 10,558 | 1,890 | 855 | `MethylDiff_Tdu_Tpr.r` and `MethylDiff_Tms_subgenome_compare_V1.r` |
| 10 | no | 3 | 1000 | 1000 | setting removed | 25% | 52,319 | 23,088 |  | `MethylDiff_Tdu_Tpr_V2.r` and `MethylDiff_Tms_subgenome_compare_V2.r` |
| 5 | yes | 10 | 300 | 300 | 1 | 25% | | | | `MethylDiff_Tdu-Tpr_CpG_mincov5_V1.r` and `MethylDiff_Tms_subgenome_compare_mincov5_V1.r` |
| 10 | yes | 10 | 1000 | 1000 | 10 | 25% | | | | `MethylDiff_Tdu-Tpr_CpG_mincov10_V1.r` and `MethylDiff_Tms_subgenome_compare_mincov10_V1.r` |
| 3 | yes | 3 | 1000 | 1000 | 10 | 25% | 7,319 | 5,996 | 4,179 | `MethylDiff_Tdu-Tpr_CpG_mincov3_V1.r` and `MethylDiff_Tms_subgenome_compare_mincov3_V1.r` |
| 3 | yes | 3 | 1000 | 1000 | 10 | 35% | 4,763 | 3,919 | 2,771 | `MethylDiff_Tdu-Tpr_CpG_mincov3_V2.r` and `MethylDiff_Tms_subgenome_compare_mincov3_V2.r` |
| 3 | yes | 3 | 1000 | 1000 | 20 | 35% | 1,487 | 1,285 | 914 | `MethylDiff_Tdu-Tpr_CpG_mincov3_V3.r` and `MethylDiff_Tms_subgenome_compare_mincov3_V3.r` |
| 3 | yes | 3 | 300 | 300 | 10 | 35% | 3,549 | 3,262 | 2,360 | `MethylDiff_Tdu-Tpr_CpG_mincov3_V4.r` and `MethylDiff_Tms_subgenome_compare_mincov3_V4.r` |

**The results from scripts `MethylDiff_Tdu-Tpr_CpG_mincov3_V4.r` and `MethylDiff_Tms_subgenome_compare_mincov3_V4.r` look reliable, so use these results for downstream analysis.**


### 2.2 DMR between diploid Tdu (2x) and the Tdu subgenome in T. miscellus (4x)
#### 2.2.1 Setting 1
| processBismarkAln mincov | using shared bases? | methRead mincov | Window size | Window step | Window cov.bases |
|--|--|--|--|--|--|
| 10 | no | 3 | 1000 | 1000 | 10 |

Scripts `MethylDiff_CpG_Tdu_vs_Tms-du_V1.r` and `MethylDiff_CpG_Tdu_vs_Tms-du_V1.sh` were used.

### 2.3 DMR between diploid Tpr (2x) and the Tpr subgenome in T. miscellus (4x)
#### 2.3.1 Setting 1
| processBismarkAln mincov | using shared bases? | methRead mincov | Window size | Window step | Window cov.bases |
|--|--|--|--|--|--|
| 10 | no | 3 | 1000 | 1000 | 10 |

Scripts `MethylDiff_CpG_Tpr_vs_Tms-pr_V1.r` and `MethylDiff_CpG_Tpr_vs_Tms-pr_V1.sh` were used.


## 3. CHG methylation
### 3.1 DMR (CHG) between two parents, two subgenomes, diploid Tdu (2x) and the Tdu subgenome, and diploid Tpr (2x) and the Tpr subgenome
Setting destrand=TRUE (the default is FALSE) will merge reads on both strands of a CpG dinucleotide. This provides better coverage, but only advised when looking at CpG methylation (for CpH methylation this will cause wrong results in subsequent analyses).

Scripts `MethylDiff_CHG_V1.r` and `MethylDiff_CHG_V1.sh` were used.

## 4. CHH methylation
### 4.1 DMR (CHH) between two parents, two subgenomes, diploid Tdu (2x) and the Tdu subgenome, and diploid Tpr (2x) and the Tpr subgenome
Setting destrand=TRUE (the default is FALSE) will merge reads on both strands of a CpG dinucleotide. This provides better coverage, but only advised when looking at CpG methylation (for CpH methylation this will cause wrong results in subsequent analyses). **This result is wrong, since we should use a lower cutoff for DMR analysis in CHH context. See [Section 1. Introduction](https://github.com/GatorShan/Tragopogon-Methylation-Project/blob/master/DMR_analysis_methylKit/README.md#1-introduction) for more information.**

Scripts `MethylDiff_CHH_V1.r` and `MethylDiff_CHH_V1.sh` were used.

