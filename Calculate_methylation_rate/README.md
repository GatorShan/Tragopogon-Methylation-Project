# Calculate methylation rate
The output files of bismark were used. For each individual, there are three files (in folder `bismark_coverage_files`):
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

## 1. Genome-wide methylation rate
Script `Genome_wide_methylation_rate_V1.ipynb` was used, which is modified from [here](https://github.com/niederhuth/Widespread-natural-variation-of-DNA-methylation-within-angiosperms/blob/c9966e4e9df6d37649c3923509874bce0dd3ad80/mC_pyTools.py#L31); the results are consistent with the report from bismark!

| Sample ID | DES1 | S1 | S2 | S3 | S4 | S5 |
| -- | -- | -- | -- | -- | -- | -- |
| CpG | 91.6% | 90.3% | 83.0% | 84.9% | 87.7% | 87.3% |
| CHG | 76.3% | 74.7% | 65.7% | 68.1% | 69.8% | 69.1% |
| CHH | 11.1% | 12.5% | 9.3% | 10.1% | 11.0% | 10.5% |

## 2. Methylation rate at different genetic features
### 2.1 Reformatting
Script `Reformat_bismark_cov_V3.py` was used. The format is changed from `bismark` to `allC`; output files are written into separate files based on the scaffold number (correct format for downstream analysis). Example script:

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
