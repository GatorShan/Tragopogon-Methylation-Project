# Calculate methylation rate
### After discussing with co-authors, to eliminate the mapping efficiency biases between different species, we think the calculation should only include loci that are shared by all samples. Previous analyses using all loci can be found [here](https://github.com/GatorShan/Tragopogon-Methylation-Project/tree/master/Calculate_methylation_rate/previous_analyses).

The starting bismark files are from [bismark analysis step 4](https://github.com/GatorShan/Tragopogon-Methylation-Project/tree/master/bismark_analysis#4-extract-bismark-methylation-profiles). For each individual, there are three files, e.g.:
  - `S2_CpG.gz.bismark.cov`
  - `S2_CHG.gz.bismark.cov`
  - `S2_CHH.gz.bismark.cov`

Format of each cov file:

`<chromosome> <start position> <end position> <methylation percentage> <count methylated> <count unmethylated>`

```
Tdub_V1_scaffold_1	436	436	100	12	0
Tdub_V1_scaffold_1	444	444	100	1	0
Tdub_V1_scaffold_1	445	445	94.1176470588235	16	1
Tdub_V1_scaffold_1	496	496	100	1	0
Tdub_V1_scaffold_1	497	497	96.2962962962963	26	1
```

## 1. Identify shared loci
I identified loci with **minimum coverage of 3 across all samples. Therefore, the different efficiency of alignment between species won't affect the results from downstream analyses**. Script `overlap_min-cov_4.py` and `overlap_min-cov_CG.sh`, for example, were used. Outputs can be found at `/orange/soltis/shan158538/Methylation_output/bismark_coverage_files`. Example:
```
DES1_CpG.gz.bismark_shared_filtered.cov
S1_CpG.gz.bismark_shared_filtered.cov
S2_CpG.gz.bismark_shared_filtered.cov
S3_CpG.gz.bismark_shared_filtered.cov
S4_CpG.gz.bismark_shared_filtered.cov
S5_CpG.gz.bismark_shared_filtered.cov
```

**Summary**
| Sample ID | CG_loci | CG_loci_shared | CHG_loci | CHG_loci_shared | CHH_loci | CHH_loci_shared |
| -- | -- | -- | -- | -- | -- | -- |
| DES1 | 40,707,284 | 8,378,772 | 34,991,145 | 7,308,535 | 200,938,529 | 39,359,361 |
| S1 | 41,715,737 | 8,378,772 | 35,899,048 | 7,308,535 | 209,104,719 | 39,359,361 |
| S2 | 17,066,881 | 8,378,772 | 14,922,781 | 7,308,535 | 83,856,112 | 39,359,361 |
| S3 | 17,373,091 | 8,378,772 | 15,145,705 | 7,308,535 | 85,113,050 | 39,359,361 |
| S4 | 41,890,596 | 8,378,772 | 36,066,828 | 7,308,535 | 210,162,748 | 39,359,361 |
| S5 | 41,680,189 | 8,378,772 | 35,885,330 | 7,308,535 | 209,237,140 | 39,359,361 |
* Note: shared loci were found in all six samples and with a min coverage of 3 in each sample; the overall loci number is from outputs in bismark analysis step 4.

## 2. Genome-wide weighted methylation level
**Using only shared loci with minimum coverage of 3 as input**, script `Genome_wide_methylation_rate_V2.ipynb` was used to calculate the genome-wide weighted methylation level. The following results were included in the manuscript. Additionally, following arcsine square root data transformation, ANOVA and post hoc Tukey analyses were employed to determine if genome-wide weighted DNA methylation levels were significantly different between species. The methylation level for T. miscellus was compared with the mid-parent value (MPV; the average of T. dubius and T. pratensis) at different cytosine contexts using a one-sample t-test. Script `Statistical_test_genome-wide_methylation_levels_transformed_data.ipynb` was used.
| Sample ID | DES1 | S1 | S2 | S3 | S4 | S5 |
| -- | -- | -- | -- | -- | -- | -- |
| CpG | 90.2% | 89.1% | 84.5% | 85.9% | 87.1% | 86.6% |
| CHG | 74.0% | 72.4% | 67.4% | 69.3% | 68.9% | 68.2% |
| CHH | 10.0% | 11.6% | 9.4% | 10.1% | 10.6% | 10.3% |

## 3. Proportions of methylated cytosine in different contexts
Extract the results from `/orange/soltis/shan158538/Methylation_output/bismark_coverage_files/Genome_wide_methylation.report_shared_loci.txt`. For the statistical test, script `Statistical_test_proportion_methylated_cytosines_transformed_data.ipynb` was used.

| Sample ID | Total number of methylated C's | Total methylated C's in CpG context | Total methylated C's in CHG context | Total methylated C's in CHH context | %CpG | %CHG | %CHH |
| -- | -- | -- | -- | -- | -- | -- | -- |
| DES1 | 495871480| 231239678| 159681869| 104949933| 46.6%| 32.2%| 21.2%|
| S1 | 812668675 | 357807030| 246129920| 208731725| 44.0%| 30.3%| 25.7%|
| S2 | 402793108 | 183181182| 127033908| 92578018| 45.5%| 31.5%| 23.0%|
| S3 | 477609409| 214420539 | 149453906| 113734964| 44.9%| 31.3%| 23.8%|
| S4 | 942064832 | 419643657 | 285325595| 237095580| 44.5%| 30.3%| 25.2%|
| S5 | 659780963 | 295351703 | 200702819| 163726441| 44.8%| 30.4%| 24.8%|

## 4. File reformatting for the following metaplot analysis
### 4.1 Bismark format to allC format
For example, script `Reformat_CG_cov_files_shared_V1.sh` was used.
### 4.2 Rename the output files
Since I got the `argument too long` error, the following script is used to shortern filenames.

```
for file in *.tsv; do mv $file "allc_DES1_$(echo $file | cut -f 6 -d "_")"; done
```

Output: e.g. `allc_DES1_10000.tsv`, in which `10000` is the scaffold ID.

## 5. Gene body methylation metaplot
![CDS_metaplot_demo](https://github.com/GatorShan/Tragopogon-Methylation-Project/blob/master/Calculate_methylation_rate/images/CDS_metaplot_demo.png)

Script **`gbm_metaplot_pe_ss.V3.py`** is used to generate the metaplot for CDS methylation level (for a gene, only count sites within CDS regions), which is modified from [here](https://github.com/bhofmei/analysis-scripts/blob/master/methyl/gbm_metaplot_pe.py) (**the changes include: 1. use all cytosine sites for upstream and downstream calculations; 2. correct the method of combining pos and neg strands methylation status in upstream and downstream regions**). Default: 1 kb upstream and downstream of gene body (delimited by the start and stop site of a gene, but not the CDS), and 20 bins for each region.

To confirm that `gbm_metaplot_pe_ss.V3.py` is correct, I tested the script using selected genes. For more information, see `Gbm_metaplot_CG_shared_V2_test3.sh` and `metaplot_test_v3.ipynb`. 

Scripts `Gbm_metaplot_CG_shared_V3.sh`, `Gbm_metaplot_CHG_shared_V3.sh`, and `Gbm_metaplot_CHH_shared_V3.sh` were used to generate the metaplot using shared sites. Input files are from step 4.2. Example outputs are `DES1_CG_new2_gbm_metaplot.tsv`, `DES1_CHG_new2_gbm_metaplot.tsv`, and `DES1_CHH_new2_gbm_metaplot.tsv`.

<img src="https://github.com/GatorShan/Tragopogon-Methylation-Project/blob/master/Calculate_methylation_rate/images/Metaplot_analysis_gene_TE_updated_01122024.png" width=750 height=900>

### Statistical test
The methylation levels were compared among the three species within each of the three regions (i.e., upstream, gene body, and downstream) using one-way ANOVA and post hoc Tukey analyses (following arcsine square root data transformation). The methylation level of T. miscellus was compared with the MPV using a one-sample t-test. Script `Statistical_test_gbm_transformed_data.ipynb` was used.

## 6. TE methylation metaplot
### 6.1 Prepare TE gff file
These results are based on the new analysis from Jon Spoelhof. 

gff files are (`/blue/soltis/shan158538/Methylation/OutPut/TE_annotation_new/`):
```
repeat_annotation_combined_TE.final.gff (all types of TEs included)
repeat_annotation_combined_Copia.final.gff
repeat_annotation_combined_Gypsy.final.gff
repeat_annotation_combined_LINE.final.gff
repeat_annotation_combined_DNA.final.gff (these are DNA transposons)
```

### 6.2 Metaplot analysis
To confirm that `TE_metaplot_pe_ss.V3.py` is correct, I tested the script using selected genes. For more information, see `TE_metaplot_CHG_shared_V1_test1.sh` and `metaplot_TE_test_v1.ipynb`. **For the actual analysis below, I used script `TE_metaplot_pe_ss.V4.py` (the upstream and downstream window size is 2,000 bp, which is 1,000 bp in V3)**.

**TE as a whole**

Scripts `TE_metaplot_CG_shared_V3.sh`, `TE_metaplot_CHG_shared_V3.sh`, and `TE_metaplot_CHH_shared_V3.sh` were used to generate the metaplot using shared sites. Example outputs are `DES1_CG_TE_new2_gbm_metaplot.tsv`, `DES1_CHG_TE_new2_gbm_metaplot.tsv`, and `DES1_CHH_TE_new2_gbm_metaplot.tsv`.

**Each type of TE**

TE types include Copia elements, Gypsy elements, LINEs, and DNA transposons. Example scripts: `Copia_metaplot_CHG_shared_V3.sh` and `DNA_metaplot_CG_shared_V3.sh`. Example outputs are `DES1_CG_LINE_new2_gbm_metaplot.tsv` and `DES1_CG_Gypsy_new2_gbm_metaplot.tsv`.

### 6.3 Statistical test
**Comparison between species within each region & between the polyploid and the MPV**

Following arcsine square root data transformation, one-way ANOVA and post hoc Tukey tests were used to compare the methylation levels among species within each region (i.e., upstream, TE, and downstream). A one-sample t-test was conducted to compare the methylation level of T. miscellus to the MPV. The following scripts were used `Statistical_test_TE_metaplot_transformed_data.ipynb`, `Statistical_test_Copia_metaplot_transformed_data.ipynb`, `Statistical_test_Gypsy_metaplot_transformed_data.ipynb`, `Statistical_test_LINE_metaplot_transformed_data.ipynb`, `Statistical_test_DNA_metaplot_transformed_data.ipynb`.

**Comparison between TE body and the flanking regions**

The following scripts were used: `Statistical_test_TE_body_elevation.ipynb`, `Statistical_test_Copia_body_elevation.ipynb`, `Statistical_test_Gypsy_body_elevation.ipynb`, `Statistical_test_LINE_body_elevation.ipynb`, `Statistical_test_DNA_body_elevation.ipynb`.







