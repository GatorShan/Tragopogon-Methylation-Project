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
### 3.3 Deduplication
Scripts `bismark_deduplicate_S4_V2.sh` and `bismark_deduplicate_S5_V2.sh` were used.

Output: `S4_cat_R1_val_1_bismark_bt2_pe_nameSorted.deduplicated.bam` and `S5_cat_R1_val_1_bismark_bt2_pe_nameSorted.deduplicated.bam`.
