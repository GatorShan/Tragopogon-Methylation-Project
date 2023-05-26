# Calculate methylation rate
The output files of bismark were used. For each individual, there are three files (in folder `/orange/soltis/shan158538/Methylation_output/bismark_coverage_files`):
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

Additionally, I identified the loci with **minimum coverage of 3 across all samples. Therefore, the different efficiency of alignment between species won't affect the results from downstream analyses**. Script `overlap_min-cov_4.py` and `overlap_min-cov_CG.sh`, for example, were used. Outputs are (using CG context as example):
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
* Note: shared loci were found in all six samples and with a min coverage of 3 in each sample.

## 1. Genome-wide weighted methylation level
Script `Genome_wide_methylation_rate_V1.ipynb` was used (using all loci), which is modified from [here](https://github.com/niederhuth/Widespread-natural-variation-of-DNA-methylation-within-angiosperms/blob/c9966e4e9df6d37649c3923509874bce0dd3ad80/mC_pyTools.py#L31); the results are consistent with the report from bismark!

| Sample ID | DES1 | S1 | S2 | S3 | S4 | S5 |
| -- | -- | -- | -- | -- | -- | -- |
| CpG | 91.6% | 90.3% | 83.0% | 84.9% | 87.7% | 87.3% |
| CHG | 76.3% | 74.7% | 65.7% | 68.1% | 69.8% | 69.1% |
| CHH | 11.1% | 12.5% | 9.3% | 10.1% | 11.0% | 10.5% |

**Using only shared loci with minimum coverage of 3 as input**, script `Genome_wide_methylation_rate_V2.ipynb` was used to calculate genome-wide weighted methylation level. The following results were included in the manuscript.
| Sample ID | DES1 | S1 | S2 | S3 | S4 | S5 |
| -- | -- | -- | -- | -- | -- | -- |
| CpG | 90.2% | 89.1% | 84.5% | 85.9% | 87.1% | 86.6% |
| CHG | 74.0% | 72.4% | 67.4% | 69.3% | 68.9% | 68.2% |
| CHH | 10.0% | 11.6% | 9.4% | 10.1% | 10.6% | 10.3% |


## 2. Methylation rate at different genetic features

| Sample ID | DES1 | S1 | S2 | S3 | S4 | S5 |
| -- | -- | -- | -- | -- | -- | -- |
| gene_CpG | 71.2% | 65.1% | 52.4% | 59.3% | 60.6% | 57.7% |
| gene_CHG | % | % | % | % | % | % |
| gene_CHH | % | % | % | % | % | % |

### 2.1 Reformatting
Scripts `Reformat_bismark_cov_V3.py` and `Reformat_CG_cov_files_V1.sh` was used. The format is changed from `bismark` to `allC`; output files are written into separate files based on the scaffold number (correct format for downstream analysis). **When only include shared loci with minimum coverage of 3**, script `Reformat_CG_cov_files_shared_V1.sh` was used. Example script:

```
cd /orange/soltis/shan158538/Methylation_output/bismark_coverage_files
Reformat_bismark_cov_V3.py DES1_CpG.gz.bismark.cov -o=/orange/soltis/shan158538/Methylation_output/feature_methylation/gene_CG/DES1/
```
Input:

`<chromosome> <start position> <end position> <methylation percentage> <count methylated> <count unmethylated>`
```
Tdub_V1_scaffold_1	14	14	66.6666666666667	2	1
Tdub_V1_scaffold_1	37	37	100	2	0
Tdub_V1_scaffold_1	38	38	100	15	0
Tdub_V1_scaffold_1	49	49	100	2	0
Tdub_V1_scaffold_1	50	50	100	17	0
```

Output (in file `allc_DES1_Tdub_V1_scaffold_1.tsv`):

`<chromosome> <position> <strand> <sequence context> <mc> <cov> <methylated>`

```
Tdub_V1_scaffold_1      14      NA      CG      2       3       NA
Tdub_V1_scaffold_1      37      NA      CG      2       2       NA
Tdub_V1_scaffold_1      38      NA      CG      15      15      NA
Tdub_V1_scaffold_1      49      NA      CG      2       2       NA
Tdub_V1_scaffold_1      50      NA      CG      17      17      NA
```
#### 2.1.1 Rename the output files
Since I got the `argument too long` error, the following script is used to shortern filenames.

```
for file in *.tsv; do mv $file "allc_DES1_$(echo $file | cut -f 6 -d "_")"; done
```

Output: e.g. `allc_DES1_10000.tsv`, in which `10000` is the scaffold ID.

#### 2.1.2 Reformat the scaffold name in the gff file
Since the `*.tsv` files have shorterned filenames, scaffold name in gff files should be changed accordingly for downstream analysis.

```
sed 's/Tdub_V1_scaffold_//' Tdub.V1.rm.gff > Tdub.V1.rm.RENAME.gff
```

### 2.2 Calculate the methylation rate within gene regions across the genome
Script `feature_methylation.ss.V1.py` (modifiled from [here](https://github.com/bhofmei/analysis-scripts/blob/master/methyl/feature_methylation.py#L264)) was used.

`MethylationRate_feature_gene_CG_V4.sh` were used, for example, to submit the job. Example ouput `DES1_gene_CG_genes.tsv`:

```
DES1	genes	CG	0.712482
```

## 3. Gene body methylation metaplot
Script `gbm_metaplot_pe_ss.V1.py` is used to generate the metaplot for CDS methylation level (for a gene, only count sites within CDS regions), which is modified from [here](https://github.com/bhofmei/analysis-scripts/blob/master/methyl/gbm_metaplot_pe.py). Default: 1 kb upstream and downstream of gene body (delimited by the start and stop site of a gene, but not the CDS), and 20 bins for each region.

Scripts ( e.g. `Gbm_metaplot_CHH_V1.sh` ) were used to submit jobs. Input files are from step 2.1.1

**When only include shared loci with minimum coverage of 3**, scripts `Gbm_metaplot_CG_shared_V1.sh`, `Gbm_metaplot_CHG_shared_V1.sh`, and `Gbm_metaplot_CHH_shared_V1.sh` were used.

**Metaplot of CG methylation**
![CG_methylation](https://github.com/GatorShan/Tragopogon-Methylation-Project/blob/master/Calculate_methylation_rate/images/CG_Feb2020.png)


**Metaplot of CHG methylation**
![CHG_methylation](https://github.com/GatorShan/Tragopogon-Methylation-Project/blob/master/Calculate_methylation_rate/images/CHG_Feb2020.png)


**Metaplot of CHH methylation**
![CHH_methylation](https://github.com/GatorShan/Tragopogon-Methylation-Project/blob/master/Calculate_methylation_rate/images/CHH_Feb2020.png)
