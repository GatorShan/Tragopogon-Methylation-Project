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
**Using only shared loci with minimum coverage of 3 as input**, script `Genome_wide_methylation_rate_V2.ipynb` was used to calculate genome-wide weighted methylation level. The following results were included in the manuscript.
| Sample ID | DES1 | S1 | S2 | S3 | S4 | S5 |
| -- | -- | -- | -- | -- | -- | -- |
| CpG | 90.2% | 89.1% | 84.5% | 85.9% | 87.1% | 86.6% |
| CHG | 74.0% | 72.4% | 67.4% | 69.3% | 68.9% | 68.2% |
| CHH | 10.0% | 11.6% | 9.4% | 10.1% | 10.6% | 10.3% |

## 3. File reformatting
### 3.1 Bismark format to allC format
For example, script `Reformat_CG_cov_files_shared_V1.sh` was used.
### 3.2 Rename the output files
Since I got the `argument too long` error, the following script is used to shortern filenames.

```
for file in *.tsv; do mv $file "allc_DES1_$(echo $file | cut -f 6 -d "_")"; done
```

Output: e.g. `allc_DES1_10000.tsv`, in which `10000` is the scaffold ID.

## 4. Gene body methylation metaplot
![CDS_metaplot_demo](https://github.com/GatorShan/Tragopogon-Methylation-Project/blob/master/Calculate_methylation_rate/images/CDS_metaplot_demo.png)

Script `gbm_metaplot_pe_ss.V1.py` is used to generate the metaplot for CDS methylation level (for a gene, only count sites within CDS regions), which is modified from [here](https://github.com/bhofmei/analysis-scripts/blob/master/methyl/gbm_metaplot_pe.py). Default: 1 kb upstream and downstream of gene body (delimited by the start and stop site of a gene, but not the CDS), and 20 bins for each region.

Scripts `Gbm_metaplot_CG_shared_V1.sh`, `Gbm_metaplot_CHG_shared_V1.sh`, and `Gbm_metaplot_CHH_shared_V1.sh` were used. Input files are from step 3.2

### CG gene metaplot

<img src="https://github.com/GatorShan/Tragopogon-Methylation-Project/blob/master/Calculate_methylation_rate/images/CG_gene_metaplot_shared_loci.png" width=600 height=400>

### CHG gene metaplot

<img src="https://github.com/GatorShan/Tragopogon-Methylation-Project/blob/master/Calculate_methylation_rate/images/CHG_gene_metaplot_shared_loci.png" width=600 height=400>

### CHH gene metaplot

<img src="https://github.com/GatorShan/Tragopogon-Methylation-Project/blob/master/Calculate_methylation_rate/images/CHH_gene_metaplot_shared_loci.png" width=600 height=400>

### Statistical test
The methylation level difference between 3 species at each bin for each context was tested using **one-way ANOVA**. The script and results are in `gbm_ANOVA_CG.ipynb`, `gbm_ANOVA_CHG.ipynb`, and `gbm_ANOVA_CHH.ipynb`. Using 0.01 as P-value cutoff, **the three species showed no difference in all bins for all context**.
