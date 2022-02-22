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
