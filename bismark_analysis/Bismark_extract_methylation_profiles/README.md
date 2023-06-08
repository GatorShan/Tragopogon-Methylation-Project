# Extract methylation profiles
## 1. CpG context
Script `bismark_methylation_extractor_S2_V2.sh` etc. was used.

Example output files
```
-rw-r--r-- 1 shan158538 soltis 718M Feb  3 13:36 CHG_OB_S2_cat_R1_val_1_bismark_bt2_pe_nameSorted.deduplicated.txt.gz
-rw-r--r-- 1 shan158538 soltis 728M Feb  3 13:36 CHG_OT_S2_cat_R1_val_1_bismark_bt2_pe_nameSorted.deduplicated.txt.gz
-rw-r--r-- 1 shan158538 soltis 2.6G Feb  3 13:36 CHH_OB_S2_cat_R1_val_1_bismark_bt2_pe_nameSorted.deduplicated.txt.gz
-rw-r--r-- 1 shan158538 soltis 2.6G Feb  3 13:36 CHH_OT_S2_cat_R1_val_1_bismark_bt2_pe_nameSorted.deduplicated.txt.gz
-rw-r--r-- 1 shan158538 soltis 766M Feb  3 13:36 CpG_OB_S2_cat_R1_val_1_bismark_bt2_pe_nameSorted.deduplicated.txt.gz
-rw-r--r-- 1 shan158538 soltis 777M Feb  3 13:36 CpG_OT_S2_cat_R1_val_1_bismark_bt2_pe_nameSorted.deduplicated.txt.gz
-rw-r--r-- 1 shan158538 soltis 105M Feb  3 18:26 S2_cat_R1_val_1_bismark_bt2_pe_nameSorted.deduplicated.bedGraph.gz
-rw-r--r-- 1 shan158538 soltis 108M Feb  3 18:26 S2_cat_R1_val_1_bismark_bt2_pe_nameSorted.deduplicated.bismark.cov.gz
-rw-r--r-- 1 shan158538 soltis  28K Feb  3 13:36 S2_cat_R1_val_1_bismark_bt2_pe_nameSorted.deduplicated.M-bias.txt
-rw-r--r-- 1 shan158538 soltis  918 Feb  3 13:36 S2_cat_R1_val_1_bismark_bt2_pe_nameSorted.deduplicated_splitting_report.txt
```

File `S2_cat_R1_val_1_bismark_bt2_pe_nameSorted.deduplicated.bismark.cov.gz` could be used for downstream analysis of **CG-DMRs**; By default, this mode will **only consider cytosines in CpG context**, but it can be extended to cytosines in any sequence context by using the option --CX

`<chromosome> <start position> <end position> <methylation percentage> <count methylated> <count unmethylated>`

```
Tdub_V1_scaffold_1	374	374	100	1	0
Tdub_V1_scaffold_1	382	382	100	1	0
Tdub_V1_scaffold_1	425	425	100	1	0
Tdub_V1_scaffold_1	426	426	100	10	0
Tdub_V1_scaffold_1	435	435	100	1	0
Tdub_V1_scaffold_1	436	436	100	12	0
Tdub_V1_scaffold_1	444	444	100	1	0
Tdub_V1_scaffold_1	445	445	94.1176470588235	16	1
Tdub_V1_scaffold_1	496	496	100	1	0
Tdub_V1_scaffold_1	497	497	96.2962962962963	26	1
```
## 2. CHG context
For example, script `bismark2bedGraph_S2_CHG_V1.sh` etc. were used.

Example outputs:
  - `S2_CHG.gz`
  - `S2_CHG.gz.bismark.cov.gz`

## 3. CHH context
For example, script `bismark2bedGraph_S1_S3_CHH_V1.sh` etc. were used.

Example outputs:
  - `DES1_CHH.gz`
  - `DES1_CHH.gz.bismark.cov`
