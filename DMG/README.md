# Differentially methylated genes (DMGs) identification and GO enrichment analysis of DMGs
## 1. DMG identification
**Genes that are overlapped with differentially methylated regions (DMRs; 300 bp windows) are defined as DMGs.**
### 1.1 Inputs
Input files are located at `/blue/soltis/shan158538/Methylation/OutPut/DMR_methylkit`. For instance, files containing CHG-DMRs include:
```
Overlap_same_direction_parents_CHG_methylDiff_all_subgenomes_CHG_methylDiff_all.txt
Unique_parents_CHG_methylDiff_all.txt
Unique_subgenomes_CHG_methylDiff_all.txt
Overlap_diff_direction_parents_CHG_methylDiff_all_subgenomes_CHG_methylDiff_all.txt
```
### 1.2 Identify DMGs
Script `DMG_V1.py` was used. Usage: `DMG_V1.py Tdub.V1.rm.gff Unique_subgenomes_methylDiff_all.txt`. The output is `DMG_Unique_subgenomes_methylDiff_all.txt`:
```
TragDub24505-RA
TragDub19346-RA
TragDub29136-RA
TragDub30060-RA
```
## 2. Prepare files for GO enrichment analysis
GOseq pipeline included in Trinity (version r20180213-2.6.5) was used for GO enrichment analysis by using GO terms derived from the annotation of T. dubius reference genome (FDR < 0.05). The method has been used in Trag inflorescence transcriptome paper and can be found [here](https://github.com/GatorShan/Tragopogon-Inflorescence-RNA-seq-Analysis/tree/master/Annotation/GO_enrichment#gene-ontology-enrichment-analysis).
### 2.1 Extract GO assignment
Xiaoxian's previous work has annotated the Tdu ref genome using Trinotate, and the result can be found in `Tdub.trinotate.annotation.xls`.
```
module load trinotate/3.0.1
$HPC_TRINOTATE_DIR/util/extract_GO_assignments_from_Trinotate_xls.pl \
  --Trinotate_xls Tdub.trinotate.annotation.xls \
  -G \
  > Tdu_go_annotation_no_ancestral.txt &
```
### 2.2 Generate gene length file
```
module load gcc/5.2.0
module load trinity/r20180213-2.6.5
${TRINITY_HOME}/util/misc/fasta_seq_length.pl Tdub.V1.transcripts.fasta > Tdub.V1.transcripts_length.txt
```
### 2.3 Generate background gene list
```
cut -f 1 Tdub.V1.transcripts_length.txt > Background_gene_id.txt
```
The output file contains all 30,325 genes.
## 3. GO enrichment analysis on DMGs
The following script was used (need to load R version 3.6 to avoid error):
```
module load gcc/5.2.0
module load trinity/r20180213-2.6.5
module load R/3.6
${TRINITY_HOME}/Analysis/DifferentialExpression/run_GOseq.pl \
    --genes_single_factor DMG_file.txt \
    --GO_assignments ../Tdu_go_annotation_no_ancestral.txt \
    --lengths ../Tdub.V1.transcripts_length.txt \
    --background ../Background_gene_id.txt
```
### 3.1 CG context
