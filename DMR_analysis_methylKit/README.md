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
| Context | processBismarkAln mincov | using shared bases? | methRead mincov | Window size | Window step | Window cov.bases | Difference | DMR between parents (A) | DMR between subgenomes | Overlapping DMR (both direaction) | Overlapping DMR (same direaction) (B) | Polyploid DMRs that are inherited from diploids (B/A%) | Scripts |
|--|--|--|--|--|--|--|--|--|--|--|--|--|--|
| CpG | 3 | yes | 3 | 300 | 300 | 10 | 35% | 3,549 | 3,262 | 2,360 | 2,300 | 64.8% | `MethylDiff_Tdu-Tpr_CpG_mincov3_V4.r` and `MethylDiff_Tms_subgenome_compare_mincov3_V4.r` |
| CHG | 3 | yes | 3 | 300 | 300 | 10 | 25% | 4,886 | 4,041 | 2,994 | 2,838 | 58.1% | `MethylDiff_Tdu-Tpr_CHG_mincov3_V5.r` and `MethylDiff_Tms_subgenome_compare_CHG_mincov3_V5.r` |
| CHH | 3 | yes | 3 | 300 | 300 | 10 | 10% | 18,143 | 14,535 | 7,676 | 7,440 | 41.0% | `MethylDiff_Tdu-Tpr_CHH_mincov3_V4.r` and `MethylDiff_Tms_subgenome_compare_CHH_mincov3_V4.r` |

Notes:

-- To decompress methylDiff output: `gunzip -c methylDiff_hyper.txt.bgz > methylDiff_hyper.txt`

-- To identify shared and unique DMRs, script `Overlap_DMR_V1.py` was used. Outputs are located at `/blue/soltis/shan158538/Methylation/OutPut/DMR_methylkit`.

## 3. DMR heatmap
<img src="https://github.com/GatorShan/Tragopogon-Methylation-Project/blob/master/DMR_analysis_methylKit/images/Heatmap_DMRs_updated.png" width=900 height=600>

Scripts `DMR_heatmap_v1.ipynb`, `DMR_heatmap_v1_CHG-Copy1.ipynb`, `DMR_heatmap_v1_CHH.ipynb`, `DMR_heatmap_v2.ipynb`, `DMR_heatmap_v2-CHG-Copy1.ipynb`, and `DMR_heatmap_v2-CHH.ipynb` were used to format the input files; scripts `DMR_heatmap.R`, `DMR_CHG_heatmap.R`, and `DMR_CHH_heatmap.R` were used to generate the heatmap. Settings when saving the figure from R can be found [here](https://github.com/GatorShan/Tragopogon-Methylation-Project/blob/master/DMR_analysis_methylKit/images/Settings_saving_figure_from_R.png).

## 4. Quantitative analysis of DMRs showing parental legacy
<img src="https://github.com/GatorShan/Tragopogon-Methylation-Project/blob/master/DMR_analysis_methylKit/images/Quantitative_analysis_DMR_parental-legacy.png" width=900 height=300>

We quantitatively examined these DMRs in each cytosine context to answer the question: how does the difference in methylation level between the diploids compare to that between the two subgenomes in the polyploid? For each DMR showing parental legacy, the difference in methylation level between the two subgenomes (represented by A) was subtracted from the difference between the diploid parents (B); the resulting absolute value (|B-A|) was used to construct the density distribution plot (Fig. S1). In CG, CHG, and CHH contexts, 8.2%, 16.4%, and 14.2% of the DMRs showing parental legacy exhibited substantial alteration in methylation level differences following polyploidy, respectively (Fig. S1). That is, |B-A| was greater than the cutoff defining DMR (i.e., 35%, 25%, and 10% in CG, CHG, and CHH contexts, respectively). Script `DMR_parental_legacy_abs_diff_density.V2.R` was used to construct the figure.

## 5. [Previous analyses (not used in the manuscript)](https://github.com/GatorShan/Tragopogon-Methylation-Project/tree/master/DMR_analysis_methylKit/Previous_analyses)
