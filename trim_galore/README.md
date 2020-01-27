# Quality control using trim galore
## 1. Description
The method of quality control is following the protocol [here](https://www.epigenesys.eu/en/protocols/bio-informatics/483-quality-control-trimming-and-alignment-of-bisulfite-seq-data-prot-57)

## 2. Procedure
Scripts `Trim_galore.V1.sh`, `Trim_galore.V2.sh`, `Trim_galore.V3.sh`, and `Trim_galore.V4.sh` were used.

```
Trimming mode: paired-end
Trim Galore version: 0.5.0
Cutadapt version: 1.18
Quality Phred score cutoff: 20
Quality encoding type selected: ASCII+33
Adapter sequence: 'AGATCGGAAGAGC' (Illumina TruSeq, Sanger iPCR; auto-detected)
Maximum trimming error rate: 0.1 (default)
Minimum required adapter overlap (stringency): 1 bp
Minimum required sequence length for both reads before a sequence pair gets removed: 20 bp
All sequences will be trimmed by 1 bp on their 3' end to avoid problems with invalid paired-end alignments with Bowtie 1
```

Lots of reads contain adapters!

```
HMCWKCCXY_s8_1_4981-LF_17_SL334590.fastq.gz_trimming_report.txt:Reads with adapters:               172,831,678 (40.9%)
HMCWKCCXY_s8_2_4981-LF_17_SL334590.fastq.gz_trimming_report.txt:Reads with adapters:               177,479,219 (42.0%)
S1_cat_R1.fastq.gz_trimming_report.txt:Reads with adapters:               250,624,136 (40.7%)
S1_cat_R2.fastq.gz_trimming_report.txt:Reads with adapters:               278,618,910 (45.3%)
S2_cat_R1.fastq.gz_trimming_report.txt:Reads with adapters:                98,456,178 (38.1%)
S2_cat_R2.fastq.gz_trimming_report.txt:Reads with adapters:               121,617,822 (47.0%)
S3_cat_R1.fastq.gz_trimming_report.txt:Reads with adapters:               127,517,814 (39.4%)
S3_cat_R2.fastq.gz_trimming_report.txt:Reads with adapters:               149,321,905 (46.2%)
S4_cat_R1.fastq.gz_trimming_report.txt:Reads with adapters:               300,021,289 (40.3%)
S4_cat_R2.fastq.gz_trimming_report.txt:Reads with adapters:               349,357,798 (47.0%)
S5_cat_R1.fastq.gz_trimming_report.txt:Reads with adapters:               186,971,237 (39.3%)
S5_cat_R2.fastq.gz_trimming_report.txt:Reads with adapters:               224,163,517 (47.2%)
```

Output files:
```
-rw-r--r-- 1 shan158538 soltis 30G Jan 24 14:00 HMCWKCCXY_s8_1_4981-LF_17_SL334590_val_1.fq.gz
-rw-r--r-- 1 shan158538 soltis 34G Jan 24 14:00 HMCWKCCXY_s8_2_4981-LF_17_SL334590_val_2.fq.gz
-rw-r--r-- 1 shan158538 soltis 43G Jan 20 03:32 S1_cat_R1_val_1.fq.gz
-rw-r--r-- 1 shan158538 soltis 39G Jan 20 03:32 S1_cat_R2_val_2.fq.gz
-rw-r--r-- 1 shan158538 soltis 18G Jan 19 05:04 S2_cat_R1_val_1.fq.gz
-rw-r--r-- 1 shan158538 soltis 17G Jan 19 05:04 S2_cat_R2_val_2.fq.gz
-rw-r--r-- 1 shan158538 soltis 23G Jan 20 02:51 S3_cat_R1_val_1.fq.gz
-rw-r--r-- 1 shan158538 soltis 21G Jan 20 02:51 S3_cat_R2_val_2.fq.gz
-rw-r--r-- 1 shan158538 soltis 52G Jan 20 12:46 S4_cat_R1_val_1.fq.gz
-rw-r--r-- 1 shan158538 soltis 47G Jan 20 12:46 S4_cat_R2_val_2.fq.gz
-rw-r--r-- 1 shan158538 soltis 33G Jan 21 21:43 S5_cat_R1_val_1.fq.gz
-rw-r--r-- 1 shan158538 soltis 31G Jan 21 21:43 S5_cat_R2_val_2.fq.gz
```

## 3. FastQC analysis of trimmed files
Script `Fastqc_quality_check.V2.sh` etc. was used.

## 4. Statistics of the trimmed files
| Species | Sequence ID | No. of paired end reads before trimming | After trim_glore |
| -- | -- | -- | -- |
| T. dubius (2x); (3060-1-4; Garfield) | DES-0001 | 422567268 | 421664516 |
| T. dubius (2x); (3040-6-2; Pullman) | S1 | 615060200 | 614845591 |
| T. pratensis (2x); (3058-1-2; Garfield) | S2 | 258651069 | 258486470 |
| T. pratensis (2x); (3058-4-10; Garfield) | S3 | 323286190 | 323032184 |
| T. miscellus (4x); (3059-7-7; Garfield) | S4 | 743663007 | 743190004 |
| T. miscellus (4x); (3059-21-5; Garfield) | S5 | 475271766 | 474830346 |

