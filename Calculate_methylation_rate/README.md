# Calculate methylation rate

### The following text showed the results from previous analysis. After discussing with co-authors, to eliminate the mapping efficiency biases between different species, we think the calculation should only include loci that are shared by all samples. The updated results can be found in directory [Shared_min-cov_sites](https://github.com/GatorShan/Tragopogon-Methylation-Project/tree/master/Calculate_methylation_rate/Shared_min-cov_sites).

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

## 1. Genome-wide weighted methylation level
Script `Genome_wide_methylation_rate_V1.ipynb` was used (using all loci), which is modified from [here](https://github.com/niederhuth/Widespread-natural-variation-of-DNA-methylation-within-angiosperms/blob/c9966e4e9df6d37649c3923509874bce0dd3ad80/mC_pyTools.py#L31); the results are consistent with the report from bismark!

| Sample ID | DES1 | S1 | S2 | S3 | S4 | S5 |
| -- | -- | -- | -- | -- | -- | -- |
| CpG | 91.6% | 90.3% | 83.0% | 84.9% | 87.7% | 87.3% |
| CHG | 76.3% | 74.7% | 65.7% | 68.1% | 69.8% | 69.1% |
| CHH | 11.1% | 12.5% | 9.3% | 10.1% | 11.0% | 10.5% |

## 2. Methylation rate at different genetic features

| Sample ID | DES1 | S1 | S2 | S3 | S4 | S5 |
| -- | -- | -- | -- | -- | -- | -- |
| gene_CpG | 71.2% | 65.1% | 52.4% | 59.3% | 60.6% | 57.7% |
| gene_CHG | % | % | % | % | % | % |
| gene_CHH | % | % | % | % | % | % |

### 2.1 Reformatting
Scripts `Reformat_bismark_cov_V3.py` and `Reformat_CG_cov_files_V1.sh` was used. The format is changed from `bismark` to `allC`; output files are written into separate files based on the scaffold number (correct format for downstream analysis). Example script:

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
![CDS_metaplot_demo](https://github.com/GatorShan/Tragopogon-Methylation-Project/blob/master/Calculate_methylation_rate/images/CDS_metaplot_demo.png)

Script `gbm_metaplot_pe_ss.V1.py` is used to generate the metaplot for CDS methylation level (for a gene, only count sites within CDS regions), which is modified from [here](https://github.com/bhofmei/analysis-scripts/blob/master/methyl/gbm_metaplot_pe.py). Default: 1 kb upstream and downstream of gene body (delimited by the start and stop site of a gene, but not the CDS), and 20 bins for each region.

Scripts ( e.g. `Gbm_metaplot_CHH_V1.sh` ) were used to submit jobs. Input files are from step 2.1.1

**Metaplot of CG methylation**
![CG_methylation](https://github.com/GatorShan/Tragopogon-Methylation-Project/blob/master/Calculate_methylation_rate/images/CG_Feb2020.png)


**Metaplot of CHG methylation**
![CHG_methylation](https://github.com/GatorShan/Tragopogon-Methylation-Project/blob/master/Calculate_methylation_rate/images/CHG_Feb2020.png)


**Metaplot of CHH methylation**
![CHH_methylation](https://github.com/GatorShan/Tragopogon-Methylation-Project/blob/master/Calculate_methylation_rate/images/CHH_Feb2020.png)
