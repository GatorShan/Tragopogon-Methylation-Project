# DMR analysis using methylKit
## 1. Introduction
**methylKit** is an R package for DNA methylation analysis and annotation from high-throughput bisulfite sequencing. Following instructions from its  [github page](https://github.com/al2na/methylKit) and [online manual](https://bioconductor.org/packages/release/bioc/vignettes/methylKit/inst/doc/methylKit.html#23_Reading_the_methylation_calls_from_sorted_Bismark_alignments).

This section is focusing on DMR analysis and comparison. The results of methylation calling using methylKit could be found [here](https://github.com/GatorShan/Tragopogon-Methylation-Project/blob/master/methylKit_analysis/README.md#2-reading-the-methylation-calls-from-sorted-bismark-alignments).

Scripts:
```r
myDiff25p.hyper_1 = getMethylDiff(myDiff_1,difference=25,qvalue=0.01,type="hyper", save.db = TRUE)
myDiff25p.hypo_1 = getMethylDiff(myDiff_1,difference=25,qvalue=0.01,type="hypo", save.db = TRUE)
myDiff25p_1 = getMethylDiff(myDiff_1,difference=25,qvalue=0.01, save.db = TRUE)
```


"difference: cutoff for ABSOLUTE VALUE of methylation percentage change between test and control (default:25)"; **THIS NUMBER SHOULD BE CHANGED FOR DIFFERENT CONTEXT!!!**

"type: one of the "hyper","hypo" or "all" strings. Specifies what type of differentially menthylated bases/regions should be returned. For retrieving Hyper-methylated regions/bases type="hyper", for hypo-methylated type="hypo" (default:"all")"

## 2. CpG methylation
### 2.1 DMR between Tdu and Tpr
#### 2.1.1 Setting 1
| processBismarkAln mincov | using shared bases? | methRead mincov | Window size | Window step | Window cov.bases |
|--|--|--|--|--|--|
| 10 | no | 3 | 1000 | 1000 | 10 |

Scripts `MethylDiff_Tdu_Tpr.r` and `MethylDiff_Tdu_Tpr.sh` were used.

#### 2.1.2 Setting 2
| processBismarkAln mincov | using shared bases? | methRead mincov | Window size | Window step | Window cov.bases |
|--|--|--|--|--|--|
| 5 | yes | 10 (default) | 300 | 300 | 1 |

Scripts `MethylDiff_Tdu-Tpr_CpG_mincov5_V1.r` and ``MethylDiff_Tdu-Tpr_CpG_mincov5_V1.sh` were used.

### 2.2 DMR between two subgenomes
#### 2.2.1 Setting 1
| processBismarkAln mincov | using shared bases? | methRead mincov | Window size | Window step | Window cov.bases |
|--|--|--|--|--|--|
| 10 | no | 3 | 1000 | 1000 | 10 |

Scripts `MethylDiff_Tms_subgenome_compare_V1.r` and `MethylDiff_Tms_subgenome_compare_V1.sh` were used.

#### 2.2.2 Setting 2
| processBismarkAln mincov | using shared bases? | methRead mincov | Window size | Window step | Window cov.bases |
|--|--|--|--|--|--|
| 5 | yes | 10 (default) | 300 | 300 | 1 |

Scripts `MethylDiff_Tms_subgenome_compare_mincov5_V1.r` and `MethylDiff_Tms_subgenome_compare_mincov5_V1.sh` were used.

#### 2.2.3 Setting 3
| processBismarkAln mincov | using shared bases? | methRead mincov | Window size | Window step | Window cov.bases |
|--|--|--|--|--|--|
| 3 | yes | 3 | 1000 | 1000 | 10 |

Scripts `MethylDiff_Tms_subgenome_compare_mincov3_V1.r` and `MethylDiff_Tms_subgenome_compare_mincov3_V1.sh` were used.

### 2.3 DMR between diploid Tdu (2x) and the Tdu subgenome in T. miscellus (4x)
#### 2.3.1 Setting 1
| processBismarkAln mincov | using shared bases? | methRead mincov | Window size | Window step | Window cov.bases |
|--|--|--|--|--|--|
| 10 | no | 3 | 1000 | 1000 | 10 |

Scripts `MethylDiff_CpG_Tdu_vs_Tms-du_V1.r` and `MethylDiff_CpG_Tdu_vs_Tms-du_V1.sh` were used.

### 2.4 DMR between diploid Tpr (2x) and the Tpr subgenome in T. miscellus (4x)
#### 2.4.1 Setting 1
| processBismarkAln mincov | using shared bases? | methRead mincov | Window size | Window step | Window cov.bases |
|--|--|--|--|--|--|
| 10 | no | 3 | 1000 | 1000 | 10 |

Scripts `MethylDiff_CpG_Tpr_vs_Tms-pr_V1.r` and `MethylDiff_CpG_Tpr_vs_Tms-pr_V1.sh` were used.



