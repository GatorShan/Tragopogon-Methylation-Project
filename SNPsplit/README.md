# SNPsplit
[SNPsplit](https://github.com/FelixKrueger/SNPsplit) is an allele-specific alignment sorter which is designed to read alignment files in SAM/ BAM format and determine the allelic origin of reads that cover known SNP3.1 positions.
It was used to examine the methylation profiles in two subgenomes of T. miscellus (4x).
## 1. Prepare the SNP file
A specific [format](https://github.com/FelixKrueger/SNPsplit/blob/master/SNPsplit_User_Guide.md#3-storing-snp-positions) of SNP file is expected. [bcftools](https://samtools.github.io/bcftools/bcftools.html#query) was used to transform
the VCF file  into user-defined formats.

Input VCF file: `Tpr_combined_filtered.PASS.vcf.gz`; from section [SNP_calling/Filter variants](https://github.com/GatorShan/Tragopogon-Methylation-Project/tree/master/SNP_calling#4-filter-variants)

The following commond was used to generate the SNP file
```bash
bcftools query -i 'GT="1/1" | GT="2/2" | GT="3/3"' -f '%ID\t%CHROM\t%POS\t%ID\t%REF/%ALT\n' Tpr_combined_filtered.PASS.vcf.gz > SNP_file.txt
### query -i: select matching samples with either GT="1/1", or GT="2/2", or GT="3/3"; GT represents genotype
### -f: format; ID (value won't be used), Chromoson, SNP position, ID (value won't be used), reference sequence/alternate sequence
```

Output file: `SNP_file.txt`. The first few lines:
```
.	Tdub_V1_scaffold_1	754	.	A/G
.	Tdub_V1_scaffold_1	2610	.	C/G
.	Tdub_V1_scaffold_1	3431	.	A/T
.	Tdub_V1_scaffold_1	3439	.	A/AT
.	Tdub_V1_scaffold_1	3447	.	T/A
```

Add a header to the `SNP_file.txt`, and the output file is `SNP_file_header.txt`, which will be used in downstream analysis.
```bash
echo -e "SNP-ID\tChromosome\tPosition\tStrand\tRef/SNP" | cat - SNP_file.txt > SNP_file_header.txt
### cat will interpret - as standard input, and will insert the output of echo before adding on the contents of SNP_file.txt
```
```
SNP-ID	Chromosome	Position	Strand	Ref/SNP
.	Tdub_V1_scaffold_1	754	.	A/G
.	Tdub_V1_scaffold_1	2610	.	C/G
.	Tdub_V1_scaffold_1	3431	.	A/T
.	Tdub_V1_scaffold_1	3439	.	A/AT
```

## 2. Mask the Tdu reference genome
### 2.1 Prepare the VCF file
Only accept homozygous alternative alleles in the VCF file
```bash
bcftools filter -i 'GT="1/1" | GT="2/2" | GT="3/3"' Tpr_combined_filtered.PASS.vcf.gz > Tpr_combined_filtered.PASS.homo.vcf.gz
```
Only accept snps from the VCF file; remove indels (the example file `mgp.v5.merged.snps_all.dbSNP142.vcf.gz` doesn't include indels)
```bash
bcftools view --types snps Tpr_combined_filtered.PASS.homo.vcf.gz > Tpr_combined_filtered.PASS.homo.snps.vcf.gz
```
Statistics of the PASS, homozygous, snps only VCF file
```
SN      0       number of samples:      1
SN      0       number of records:      2423349
SN      0       number of no-ALTs:      0
SN      0       number of SNPs: 2100512
SN      0       number of MNPs: 0
SN      0       number of indels:       322837
SN      0       number of others:       0
SN      0       number of multiallelic sites:   0
SN      0       number of multiallelic SNP sites:       0
```
### 2.2 Mask the Tdu genome
[bedtools](https://bedtools.readthedocs.io/en/latest/content/tools/maskfasta.html) was used to mask the Tdu genome. The script is `Mask.Tdu.Genome_V1.sh`(very fast!).
Output file is `Tdub.V1.masked.fasta`.

## 3. Bismark mapping
### 3.1 Index the masked Tdu genome
Script `bismark_genome_prep_V2.sh` was used.
### 3.2 Alignment
Scripts `bismark_alignment_S4_V2.sh` and `bismark_alignment_S5_V2.sh` were used.

Alignment report for S4:
```
Final Alignment report
======================
Sequence pairs analysed in total:       743190004
Number of paired-end alignments with a unique best hit: 244864209
Mapping efficiency:     32.9% 
Sequence pairs with no alignments under any condition:  452401417
Sequence pairs did not map uniquely:    45924378
Sequence pairs which were discarded because genomic sequence could not be extracted:    18884

Final Cytosine Methylation Report
=================================
Total number of C's analysed:   11553845192

Total methylated C's in CpG context:    1499546151
Total methylated C's in CHG context:    1020986189
Total methylated C's in CHH context:    906489816
Total methylated C's in Unknown context:        20063120

Total unmethylated C's in CpG context:  220227247
Total unmethylated C's in CHG context:  448522862
Total unmethylated C's in CHH context:  7458072927
Total unmethylated C's in Unknown context:      48189205

C methylated in CpG context:    87.2%
C methylated in CHG context:    69.5%
C methylated in CHH context:    10.8%
C methylated in unknown context (CN or CHN):    29.4%
```

Alignment report for S5:
```
Final Alignment report
======================
Sequence pairs analysed in total:       474830346
Number of paired-end alignments with a unique best hit: 151115840
Mapping efficiency:     31.8% 
Sequence pairs with no alignments under any condition:  295055651
Sequence pairs did not map uniquely:    28658855
Sequence pairs which were discarded because genomic sequence could not be extracted:    11358

Final Cytosine Methylation Report
=================================
Total number of C's analysed:   7301285444

Total methylated C's in CpG context:    931114288
Total methylated C's in CHG context:    635993158
Total methylated C's in CHH context:    547809871
Total methylated C's in Unknown context:        12420757

Total unmethylated C's in CpG context:  143822815
Total unmethylated C's in CHG context:  289407940
Total unmethylated C's in CHH context:  4753137372
Total unmethylated C's in Unknown context:      31546186

C methylated in CpG context:    86.6%
C methylated in CHG context:    68.7%
C methylated in CHH context:    10.3%
C methylated in unknown context (CN or CHN):    28.3%
```

### 3.3 Deduplication
Scripts `bismark_deduplicate_S4_V2.sh` and `bismark_deduplicate_S5_V2.sh` were used.

Output: `S4_cat_R1_val_1_bismark_bt2_pe_nameSorted.deduplicated.bam` and `S5_cat_R1_val_1_bismark_bt2_pe_nameSorted.deduplicated.bam`.


Report for S4:
```
Total number duplicated alignments removed:     68022816 (27.78%)
Duplicated alignments were found at:    49483432 different position(s)
Total count of deduplicated leftover sequences: 176822509 (72.22% of total)
```


Report for S5:
```
Total number duplicated alignments removed:     28086892 (18.59%)
Duplicated alignments were found at:    22611322 different position(s)
Total count of deduplicated leftover sequences: 123017590 (81.41% of total)
```

## 4. Split bam files from Bismark alignment
**Bam files don't have to be splited for downstream analysis. The purpose of remaining this section is showing an example of how to split a bam file.**

**When runing SNPsplit, remove flag --verbose will accelerate the speed a LOT.**

As it took forever to run SNPsplit when using a single huge bam file (**this is not ture; see above**), the output bam file from Bismark alignment was divided into small subfiles based on scaffold ID. There are 12,310 scaffolds in total, so there are 12,310 subfiles.

### 4.1 Sort and index the bam file
Scripts `Bam_formatting_split_V1.sh` and `Bam_formatting_split_V2.sh` were used.
### 4.2 Split the bam file
Scripts `Bam_split_samtools_V1.sh` and `Bam_split_samtools_V2.sh` were used.

## 5. Running SNPsplit
Scripts `SNPsplit_S4_new_V1.sh` and `SNPsplit_S5_new_V1.sh` were used.

**Important outputs** are:
```
S4_cat_R1_val_1_bismark_bt2_pe_nameSorted.deduplicated.genome1.bam
S4_cat_R1_val_1_bismark_bt2_pe_nameSorted.deduplicated.genome2.bam
S5_cat_R1_val_1_bismark_bt2_pe_nameSorted.deduplicated.genome1.bam
S5_cat_R1_val_1_bismark_bt2_pe_nameSorted.deduplicated.genome2.bam
```
Notes:
| ID | Meaning |
| - | - |
| S4 | T. miscellus (4x); (3059-7-7; Garfield) |
| S5 | T. miscellus (4x); (3059-21-5; Garfield) |
| genome 1 | Tdu-homeolog; reference |
| genome 2 | Tpr-homeolog |


Allele-tagging report for S4:
```
Processed 353645018 read alignments in total
Reads were unaligned and hence skipped: 0 (0.00%)
292910221 reads were unassignable (82.83%)
33659280 reads were specific for genome 1 (9.52%)
25647596 reads were specific for genome 2 (7.25%)
13680625 reads that were unassignable contained C>T SNPs preventing the assignment
763067 reads did not contain one of the expected bases at known SNP positions (0.22%)
1427921 contained conflicting allele-specific SNPs (0.40%)
```

Allele-tagging report for S5:
```
Processed 246035180 read alignments in total
Reads were unaligned and hence skipped: 0 (0.00%)
203138440 reads were unassignable (82.56%)
23342868 reads were specific for genome 1 (9.49%)
18563887 reads were specific for genome 2 (7.55%)
9724257 reads that were unassignable contained C>T SNPs preventing the assignment
527280 reads did not contain one of the expected bases at known SNP positions (0.21%)
989985 contained conflicting allele-specific SNPs (0.40%)
```
## 6. DMR analysis
[methylkit](https://github.com/al2na/methylKit) was used to perform the DMR analysis
### 6.1 Sort and index bam files
Script `Bismark_bam_formatting_V2.sh` was used. Outputs are:
```
S4_cat_R1_val_1_bismark_bt2_pe_nameSorted.deduplicated.genome1_PosSorted.bam
S5_cat_R1_val_1_bismark_bt2_pe_nameSorted.deduplicated.genome1_PosSorted.bam
S4_cat_R1_val_1_bismark_bt2_pe_nameSorted.deduplicated.genome2_PosSorted.bam
S5_cat_R1_val_1_bismark_bt2_pe_nameSorted.deduplicated.genome2_PosSorted.bam
```
### 6.2 Read the methylation calls from sorted Bismark alignments
| Sample.id | Original file |
| - | - |
| Tms_1_du | S4_cat_R1_val_1_bismark_bt2_pe_nameSorted.deduplicated.genome1_PosSorted.bam |
| Tms_2_du | S5_cat_R1_val_1_bismark_bt2_pe_nameSorted.deduplicated.genome1_PosSorted.bam |
| Tms_1_pr | S4_cat_R1_val_1_bismark_bt2_pe_nameSorted.deduplicated.genome2_PosSorted.bam |
| Tms_2_pr | S5_cat_R1_val_1_bismark_bt2_pe_nameSorted.deduplicated.genome2_PosSorted.bam |

#### 6.2.1 CpG methylation
Scripts `processBismarkAln_subgenome_compare_CG_V1.r` and `processBismarkAln_subgenome_compare_CG_V1.sh` were used. By default, minimum read coveage per base is 10.

### 6.3 DMR (CpG) between two subgenomes in T. miscellus (4x)
Scripts `MethylDiff_Tms_subgenome_compare_V1.r` and `MethylDiff_Tms_subgenome_compare_V1.sh` were used.

### 6.4 DMR (CpG) between diploid Tdu (2x) and the Tdu subgenome in T. miscellus (4x)
Scripts `MethylDiff_CpG_Tdu_vs_Tms-du_V1.r` and `MethylDiff_CpG_Tdu_vs_Tms-du_V1.sh` were used.

### 6.5 DMR (CpG) between diploid Tpr (2x) and the Tpr subgenome in T. miscellus (4x)
Scripts `MethylDiff_CpG_Tpr_vs_Tms-pr_V1.r` and `MethylDiff_CpG_Tpr_vs_Tms-pr_V1.sh` were used.




