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
## 2. GO enrichment analysis of DMGs
