# Bismark analysis
Bismark is a set of tools for the time-efficient analysis of Bisulfite-Seq (BS-Seq) data. Bismark performs alignments of bisulfite-treated reads to a reference genome and cytosine methylation calls at the same time. The protocol of the analysis could be found [here](https://github.com/FelixKrueger/Bismark/tree/master/Docs).

## 1. Genome indexing
Bismark will create two individual folders: one for a C->T converted genome and the other one for the G->A converted genome. After creating C->T and G->A versions of the genome they will be indexed in parallel.

Script `bismark_genome_prep_V1.sh` was used.

## 2. Alignment
For S1, script `bismark_alignment_S1_V1.sh` was used.
- ```
    Final Alignment report
    ======================
    Sequence pairs analysed in total:       614845591
    Number of paired-end alignments with a unique best hit: 291340988
    Mapping efficiency:     47.4% 
    Sequence pairs with no alignments under any condition:  277397113
    Sequence pairs did not map uniquely:    46107490
    Sequence pairs which were discarded because genomic sequence could not be extracted:    20738
    ```

For S2, script `bismark_alignment_S2_V1.sh` was used.
  - ```
    Final Alignment report
    ======================
    Sequence pairs analysed in total:       258486470
    Number of paired-end alignments with a unique best hit: 46100560
    Mapping efficiency:     17.8% 
    Sequence pairs with no alignments under any condition:  198923920
    Sequence pairs did not map uniquely:    13461990
    Sequence pairs which were discarded because genomic sequence could not be extracted:    2760
    ```

**As the mapping efficiency is very different between Tdu and Tpr samples, I want to use a less strigent mapping method, and see if the Tpr mapping efficiency would increase a little bit**

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
| Mapping efficiency | 47.4% | 44.2% | |
| Time | 3d 8h 19m | 8h 26m | |
| Results | C methylated in CpG context:    89.7% C methylated in CHG context:    74.3% C methylated in CHH context:    12.4% C methylated in unknown context (CN or CHN):    16.2% | C methylated in CpG context:    89.7% C methylated in CHG context:    74.3% C methylated in CHH context:    12.4% C methylated in unknown context (CN or CHN):    16.5% | |


For S2:

| Alignment job ID | S2_V1 | S2_V2 | S2_V3 |
| -- | -- | -- | -- |
| Parameters | default (100% reads) | -N 1 (10% reads) | --score_min L,0,-0.6 (10% reads) |
| Mapping efficiency | 17.8% | 13.4% | 50.9% |
| Time | 1d 8h 52m | 3h 49m | 6h 20m |
| Results | C methylated in CpG context:    82.6% C methylated in CHG context:    65.4% C methylated in CHH context:    9.1% C methylated in unknown context (CN or CHN):    14.2% | C methylated in CpG context:    81.5% C methylated in CHG context:    63.8% C methylated in CHH context:    8.9% C methylated in unknown context (CN or CHN):    13.8% | C methylated in CpG context:    78.7% C methylated in CHG context:    62.6% C methylated in CHH context:    9.7% C methylated in unknown context (CN or CHN):    16.1% |

