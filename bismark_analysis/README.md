# Bismark analysis
Bismark is a set of tools for the time-efficient analysis of Bisulfite-Seq (BS-Seq) data. Bismark performs alignments of bisulfite-treated reads to a reference genome and cytosine methylation calls at the same time. The protocol of the analysis could be found [here](https://github.com/FelixKrueger/Bismark/tree/master/Docs).

## 1. Genome indexing
Bismark will create two individual folders: one for a C->T converted genome and the other one for the G->A converted genome. After creating C->T and G->A versions of the genome they will be indexed in parallel.

Script `bismark_genome_prep_V1.sh` was used.

## 2. Alignment
For S1, script `bismark_alignment_S1_V1.sh` was used.


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
