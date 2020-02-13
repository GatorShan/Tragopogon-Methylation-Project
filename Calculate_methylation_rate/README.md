# Calculate methylation rate
The output files of bismark were used. For each individual, there are three files (in folder `bismark_coverage_files`):
  - `S2_CpG.gz.bismark.cov`
  - `S2_CHG.gz.bismark.cov`
  - `S2_CHH.gz.bismark.cov`

Format of each cov file:
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

