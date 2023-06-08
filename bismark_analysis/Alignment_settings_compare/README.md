## Use less strigent mapping method for S1 and S2
### 1. Subsampling
Script `Subsampling_V2.sh` was used to extract 10% of reads randomly from the original fastq.gz files.

Output:
  - `S1_cat_R2_val_2_subset_0.1.fq.gz` and `S1_cat_R1_val_1_subset_0.1.fq.gz`
  - `S2_cat_R1_val_1_subset_0.1.fq.gz` and `S2_cat_R2_val_2_subset_0.1.fq.gz`

### 2. Bismak alignment

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

