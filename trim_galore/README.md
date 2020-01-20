# Quality control using trim galore
## 1. Description
The method of quality control is following the protocol [here](https://www.epigenesys.eu/en/protocols/bio-informatics/483-quality-control-trimming-and-alignment-of-bisulfite-seq-data-prot-57)

## 2. Procedure
Scripts `Trim_galore.V1.sh`, `Trim_galore.V2.sh`, and `Trim_galore.V3.sh` were used.

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
S1_cat_R1.fastq.gz_trimming_report.txt:Reads with adapters:               250,624,136 (40.7%)
S1_cat_R2.fastq.gz_trimming_report.txt:Reads with adapters:               278,618,910 (45.3%)
S2_cat_R1.fastq.gz_trimming_report.txt:Reads with adapters:                98,456,178 (38.1%)
S2_cat_R2.fastq.gz_trimming_report.txt:Reads with adapters:               121,617,822 (47.0%)
S3_cat_R1.fastq.gz_trimming_report.txt:Reads with adapters:               127,517,814 (39.4%)
S3_cat_R2.fastq.gz_trimming_report.txt:Reads with adapters:               149,321,905 (46.2%)
S4_cat_R1.fastq.gz_trimming_report.txt:Reads with adapters:               300,021,289 (40.3%)
S4_cat_R2.fastq.gz_trimming_report.txt:Reads with adapters:               349,357,798 (47.0%)
```

Fastqc analysis was performed separately using script `Fastqc_quality_check.V2.sh`.

Output files:
```
43G Jan 20 03:32 S1_cat_R1_val_1.fq.gz
39G Jan 20 03:32 S1_cat_R2_val_2.fq.gz
18G Jan 19 05:04 S2_cat_R1_val_1.fq.gz
17G Jan 19 05:04 S2_cat_R2_val_2.fq.gz
23G Jan 20 02:51 S3_cat_R1_val_1.fq.gz
21G Jan 20 02:51 S3_cat_R2_val_2.fq.gz
```
