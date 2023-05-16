# Bismark analysis
Bismark is a set of tools for the time-efficient analysis of Bisulfite-Seq (BS-Seq) data. Bismark performs alignments of bisulfite-treated reads to a reference genome and cytosine methylation calls at the same time. The protocol of the analysis could be found [here](https://github.com/FelixKrueger/Bismark/tree/master/Docs).

## 1. Genome indexing
Bismark will create two individual folders: one for a C->T converted genome and the other one for the G->A converted genome. After creating C->T and G->A versions of the genome they will be indexed in parallel.

Script `bismark_genome_prep_V1.sh` was used.

## 2. Alignment

**Here are the summary report of the alignment. We here used the default settings, for more detail about comparing different settings find section 2.1**

| Sample ID | DES1 | S1 | S2 | S3 | S4 | S5 |
| -- | -- | -- | -- | -- | -- | -- |
| Species | T. dubius (2x); (3060-1-4; Garfield) | T. dubius (2x); (3040-6-2; Pullman) | T. pratensis (2x); (3058-1-2; Garfield) | T. pratensis (2x); (3058-4-10; Garfield) | T. miscellus (4x); (3059-7-7; Garfield) | T. miscellus (4x); (3059-21-5; Garfield) |
| Script used | `bismark_alignment_DES1_V1.sh` | `bismark_alignment_S1_V1.sh` | `bismark_alignment_S2_V1.sh` | `bismark_alignment_S3_V1.sh` | `bismark_alignment_S4_V1.sh` | `bismark_alignment_S5_V1.sh` |
| Mapping efficiency | 44.4% | 47.4% | 17.8% | 17.3% | 32.2% | 31.0% |
| Time | 2d 10h | 3d 8h | 1d 8h | 1d 15h | 4d 1h | 1d 0h |
| Results | C methylated in CpG context:    91.1% C methylated in CHG context:    76.1% C methylated in CHH context:    11.2% C methylated in unknown context (CN or CHN):    15.5% | C methylated in CpG context:    89.7% C methylated in CHG context:    74.3% C methylated in CHH context:    12.4% C methylated in unknown context (CN or CHN):    16.2% | C methylated in CpG context:    82.6% C methylated in CHG context:    65.4% C methylated in CHH context:    9.1% C methylated in unknown context (CN or CHN):    14.2% | C methylated in CpG context:    84.5% C methylated in CHG context:    67.8% C methylated in CHH context:    9.9% C methylated in unknown context (CN or CHN):    15.2% | C methylated in CpG context:    87.2% C methylated in CHG context:    69.5% C methylated in CHH context:    10.9% C methylated in unknown context (CN or CHN):    15.1% | C methylated in CpG context:    86.7% C methylated in CHG context:    68.8% C methylated in CHH context:    10.4% C methylated in unknown context (CN or CHN):    14.6% |


### 2.1 Use less strigent mapping method for S1 and S2
#### 2.2.1 Subsampling
Script `Subsampling_V2.sh` was used to extract 10% of reads randomly from the original fastq.gz files.

Output:
  - `S1_cat_R2_val_2_subset_0.1.fq.gz` and `S1_cat_R1_val_1_subset_0.1.fq.gz`
  - `S2_cat_R1_val_1_subset_0.1.fq.gz` and `S2_cat_R2_val_2_subset_0.1.fq.gz`

#### 2.2.2 bismak alignment

For S1:

| Alignment job ID | S1_V1 | S1_V2 | S1_V3 |
| -- | -- | -- | -- |
| Parameters | default (100% reads) | -N 1 (10% reads) | --score_min L,0,-0.6 (10% reads) |
| Mapping efficiency | 47.4% | 44.2% | 65.8% |
| Time | 3d 8h 19m | 8h 26m | 15h 34m |
| Results | C methylated in CpG context:    89.7% C methylated in CHG context:    74.3% C methylated in CHH context:    12.4% C methylated in unknown context (CN or CHN):    16.2% | C methylated in CpG context:    89.7% C methylated in CHG context:    74.3% C methylated in CHH context:    12.4% C methylated in unknown context (CN or CHN):    16.5% | C methylated in CpG context:    86.6% C methylated in CHG context:    71.4% C methylated in CHH context:    12.3% C methylated in unknown context (CN or CHN):    17.5% |


For S2:

| Alignment job ID | S2_V1 | S2_V2 | S2_V3 |
| -- | -- | -- | -- |
| Parameters | default (100% reads) | -N 1 (10% reads) | --score_min L,0,-0.6 (10% reads) |
| Mapping efficiency | 17.8% | 13.4% | 50.9% |
| Time | 1d 8h 52m | 3h 49m | 6h 20m |
| Results | C methylated in CpG context:    82.6% C methylated in CHG context:    65.4% C methylated in CHH context:    9.1% C methylated in unknown context (CN or CHN):    14.2% | C methylated in CpG context:    81.5% C methylated in CHG context:    63.8% C methylated in CHH context:    8.9% C methylated in unknown context (CN or CHN):    13.8% | C methylated in CpG context:    78.7% C methylated in CHG context:    62.6% C methylated in CHH context:    9.7% C methylated in unknown context (CN or CHN):    16.1% |

**After talking with Bob, we prefer to use a more conservative mapping paramter; so results from V1 mapping scripts will be used for now!**

**"I would suggest you proceed cautiously.  More is not always better, especially if reads are being misplaced. This is up the user and dependent on the questions you are asking. In general, we mostly take a very conservative approach"**

## 3. Deduplication
Remove alignments to the same position in the genome from the Bismark mapping output, which can arise by e.g excessive PCR amplification; by default, the first alignment to a given position will be used irrespective of its methylation call.

Script `bismark_deduplicate_S2_V3.sh` etc. were used.

| Sample ID | DES1 | S1 | S2 | S3 | S4 | S5 |
| -- | -- | -- | -- | -- | -- | -- |
| % of deduplicated leftover sequences | 65.82% | 72.50% | 86.51% | 83.89% | 72.30% | 81.51% |

## 4. Extract Bismark methylation profiles
### 4.1 CpG context
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
### 4.2 CHG context
For example, script `bismark2bedGraph_S2_CHG_V1.sh` etc. were used.

Example outputs:
  - `S2_CHG.gz`
  - `S2_CHG.gz.bismark.cov.gz`

### 4.3 CHH context
For example, script `bismark2bedGraph_S1_S3_CHH_V1.sh` etc. were used.

Example outputs:
  - `DES1_CHH.gz`
  - `DES1_CHH.gz.bismark.cov`

## 5. Average depth at cytosine sites
The python code `Cytosine_coverage_V1.py` and associated bash code `Cytosine_sites_coverage_V1.sh` and `Cytosine_sites_coverage_V2.sh` were used.

| Sample ID | DES1 | S1 | S2 | S3 | S4 | S5 |
| -- | -- | -- | -- | -- | -- | -- |
| Average depth at cytosine sites | 17.0 | 27.2 | 14.0 | 15.5 | 22.8 | 16.4 |
 
## 6. Non-conversion rate
Spiked in unmethylated lambda DNA was used to calculate non-conversion rate.

