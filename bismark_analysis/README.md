# Bismark analysis
Bismark is a set of tools for the time-efficient analysis of Bisulfite-Seq (BS-Seq) data. Bismark performs alignments of bisulfite-treated reads to a reference genome and cytosine methylation calls at the same time. The protocol of the analysis could be found [here](https://github.com/FelixKrueger/Bismark/tree/master/Docs).

## 1. Genome indexing
Bismark will create two individual folders: one for a C->T converted genome and the other one for the G->A converted genome. After creating C->T and G->A versions of the genome they will be indexed in parallel.

Script `bismark_genome_prep_V1.sh` was used.

## 2. Alignment

**Here are the summary report of the alignment. We here used the default settings, for more detail about comparing different settings see [here](https://github.com/GatorShan/Tragopogon-Methylation-Project/tree/master/bismark_analysis/Alignment_settings_compare)**

| Sample ID | DES1 | S1 | S2 | S3 | S4 | S5 |
| -- | -- | -- | -- | -- | -- | -- |
| Species | T. dubius (2x); (3060-1-4; Garfield) | T. dubius (2x); (3040-6-2; Pullman) | T. pratensis (2x); (3058-1-2; Garfield) | T. pratensis (2x); (3058-4-10; Garfield) | T. miscellus (4x); (3059-7-7; Garfield) | T. miscellus (4x); (3059-21-5; Garfield) |
| Script used | `bismark_alignment_DES1_V1.sh` | `bismark_alignment_S1_V1.sh` | `bismark_alignment_S2_V1.sh` | `bismark_alignment_S3_V1.sh` | `bismark_alignment_S4_V1.sh` | `bismark_alignment_S5_V1.sh` |
| Mapping efficiency | 44.4% | 47.4% | 17.8% | 17.3% | 32.2% | 31.0% |
| Time | 2d 10h | 3d 8h | 1d 8h | 1d 15h | 4d 1h | 1d 0h |
| Results | C methylated in CpG context:    91.1% C methylated in CHG context:    76.1% C methylated in CHH context:    11.2% C methylated in unknown context (CN or CHN):    15.5% | C methylated in CpG context:    89.7% C methylated in CHG context:    74.3% C methylated in CHH context:    12.4% C methylated in unknown context (CN or CHN):    16.2% | C methylated in CpG context:    82.6% C methylated in CHG context:    65.4% C methylated in CHH context:    9.1% C methylated in unknown context (CN or CHN):    14.2% | C methylated in CpG context:    84.5% C methylated in CHG context:    67.8% C methylated in CHH context:    9.9% C methylated in unknown context (CN or CHN):    15.2% | C methylated in CpG context:    87.2% C methylated in CHG context:    69.5% C methylated in CHH context:    10.9% C methylated in unknown context (CN or CHN):    15.1% | C methylated in CpG context:    86.7% C methylated in CHG context:    68.8% C methylated in CHH context:    10.4% C methylated in unknown context (CN or CHN):    14.6% |

## 3. Deduplication
Remove alignments to the same position in the genome from the Bismark mapping output, which can arise by e.g excessive PCR amplification; by default, the first alignment to a given position will be used irrespective of its methylation call.

Script `bismark_deduplicate_S2_V3.sh` etc. were used.

| Sample ID | DES1 | S1 | S2 | S3 | S4 | S5 |
| -- | -- | -- | -- | -- | -- | -- |
| % of deduplicated leftover sequences | 65.82% | 72.50% | 86.51% | 83.89% | 72.30% | 81.51% |

## 4. Extract Bismark methylation profiles
For methods, please find it [here](https://github.com/GatorShan/Tragopogon-Methylation-Project/tree/master/bismark_analysis/Bismark_extract_methylation_profiles). The results can be found in `/orange/soltis/shan158538/Methylation_output/bismark_coverage_files`. Important outputs are:
```
DES1_CpG.gz.bismark.cov
DES1_CHG.gz.bismark.cov
DES1_CHH.gz.bismark.cov

S1_CpG.gz.bismark.cov
S1_CHG.gz.bismark.cov
S1_CHH.gz.bismark.cov

S2_CpG.gz.bismark.cov
S2_CHG.gz.bismark.cov
S2_CHH.gz.bismark.cov

S3_CpG.gz.bismark.cov
S3_CHG.gz.bismark.cov
S3_CHH.gz.bismark.cov

S4_CpG.gz.bismark.cov
S4_CHG.gz.bismark.cov
S4_CHH.gz.bismark.cov

S5_CpG.gz.bismark.cov
S5_CHG.gz.bismark.cov
S5_CHH.gz.bismark.cov
```

## 5. Average depth at cytosine sites
The python code `Cytosine_coverage_V1.py` and associated bash code `Cytosine_sites_coverage_V1.sh` and `Cytosine_sites_coverage_V2.sh` were used.

| Sample ID | DES1 | S1 | S2 | S3 | S4 | S5 |
| -- | -- | -- | -- | -- | -- | -- |
| Average depth at cytosine sites | 17.0 | 27.2 | 14.0 | 15.5 | 22.8 | 16.4 |
 
## 6. Bisulfite conversion rate
| Sample ID | DES1 | S1 | S2 | S3 | S4 | S5 |
| -- | -- | -- | -- | -- | -- | -- |
| Bisulfite conversion rate | 99.8% | 99.6% | 99.6% | 99.6% | 99.6% | 99.6% |

Spiked in unmethylated [lambda DNA](https://github.com/GatorShan/Tragopogon-Methylation-Project/blob/master/bismark_analysis/lambda.fa) was used to calculate bisulfite conversion rate. If 100% of cytosines in the lambda DNA were converted to T after bisulfite treatment, the methylation level should be 0%.

- Script `bismark_genome_prep_V2.sh` was used to convert the lambda reference DNA into two different bisulfite converted versions and index them for alignment with Bowtie 2
- Script `bismark_alignment_DES1_lambda.sh` etc. were used to map the reads from each sample to the lambda reference DNA
- The alignment reports, e.g. `HMCWKCCXY_s8_1_4981-LF_17_SL334590_val_1_bismark_bt2_PE_report.txt`, were used to calculte the bisulfite conversion rate for each sample

An example alignment report (from sample DES1):
```
Total number of C's analysed:	14970774

Total methylated C's in CpG context:	7745
Total methylated C's in CHG context:	8952
Total methylated C's in CHH context:	18384
Total methylated C's in Unknown context:	1

Total unmethylated C's in CpG context:	3428445
Total unmethylated C's in CHG context:	3592318
Total unmethylated C's in CHH context:	7914930
Total unmethylated C's in Unknown context:	1589
```
The bisulfite conversion rate = Total unmethylated C's / Total number of C's analysed = (3428445+3592318+7914930)/14970774 = 99.8%

Total number of C's analysed = 14970774 = 7745+8952+18384+3428445+3592318+7914930 (methylated C's in Unknown context were not counted)
