# *Tragopogon* Methylation Project
<p align="center">
<img src="https://github.com/GatorShan/Tragopogon-Methylation-Project/blob/master/Calculate_methylation_rate/images/Sampling.png" width=320 height=300>
</p>

<p align="center">
<img src="https://github.com/GatorShan/Tragopogon-Methylation-Project/blob/master/Calculate_methylation_rate/images/WGBS.png" width=750 height=330>
</p>

## 1. Materials
| Species | Sequence ID |
| -- | -- |
| *T. dubius* (2*x*); (3060-1-4; Garfield) | DES1 |
| *T. dubius* (2*x*); (3040-6-2; Pullman) | S1 |
| *T. pratensis* (2*x*); (3058-1-2; Garfield) | S2 |
| *T. pratensis* (2*x*); (3058-4-10; Garfield) | S3 |
| *T. miscellus* (4*x*); (3059-7-7; Garfield) | S4 |
| *T. miscellus* (4*x*); (3059-21-5; Garfield) | S5 |

## 2. Data trimming
The trimming process and statistics can be found in the [trim_galore section](https://github.com/GatorShan/Tragopogon-Methylation-Project/tree/master/trim_galore).

## 3. Alignment
[Bismark analysis section](https://github.com/GatorShan/Tragopogon-Methylation-Project/tree/master/bismark_analysis) showed the alignment process.

## 4. Methylation level
Section [Calculate_methylation_rate](https://github.com/GatorShan/Tragopogon-Methylation-Project/tree/master/Calculate_methylation_rate) includes analyses of genome-wide methylation rate, proportions of methylated cytosines in different contexts, methylation metaplots at gene bodies and TEs (including various types) and their flanking regions.

## 5. Differentiation of the subgenome origin of *T. miscellus* reads
Section [SNP_calling](https://github.com/GatorShan/Tragopogon-Methylation-Project/tree/master/SNP_calling) identified SNPs between *T. dubius* and *T. pratensis*. Section [SNPsplit](https://github.com/GatorShan/Tragopogon-Methylation-Project/tree/master/SNPsplit) differentiate *T. dubius*- and *T.pratensis*-derived reads in *T. miscellus*.

## 6. DMR and DMG analysis
Section [methylKit_analysis](https://github.com/GatorShan/Tragopogon-Methylation-Project/tree/master/methylKit_analysis) described the process of reading methylation calls from sorted Bismark alignment files using methylKit. Section [DMR_analysis_methylKit](https://github.com/GatorShan/Tragopogon-Methylation-Project/tree/master/DMR_analysis_methylKit) described the process of DMR analysis and comparison. Section [DMG](https://github.com/GatorShan/Tragopogon-Methylation-Project/tree/master/DMG) includes the identification of DMGs and GO enrichment analysis of these DMGs.
