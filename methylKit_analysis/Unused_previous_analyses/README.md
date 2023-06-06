## 3. Unused previous analyses

### 3.1 Read Bismark alignment for T. dubius, T. pratensis, and T. miscellus; mincov = 10
Sort and index bam files first; script `Bismark_bam_formatting_V1.sh` was used. Input files are deduplicated. **By default, the minimum read coverage to call a mehtylation status for a base is 10 (mincov = 10)**

#### 3.1.1 CpG methylation
Script `processBismarkAln_CG_V1.r` and `processBismarkAln_CG_V1.sh` were used.

Example saved/output files:
```
Tdu_1_CpG.txt
Tdu_1_CpG_conversionStats.txt
Tdu_2_CpG.txt
Tdu_2_CpG_conversionStats.txt
```
Descriptive statistics on samples: percent methylation and read coverage information were calculated using scripts `methylKit_DescriptiveStatistics.r` and `methylKit_DescriptiveStatistics.sh`. Results could be found in files `methylKit_DescriptiveStatistics.pdf` and `methylKit_DescriptiveStatistics_17957860.out`.

### 3.2 Read Bismark alignment for two subgenomes; mincov = 10
Start with **results from [SNPsplit](https://github.com/GatorShan/Tragopogon-Methylation-Project/blob/master/SNPsplit/README.md#5-running-snpsplit)**, script `Bismark_bam_formatting_V2.sh` was used to sort and index bam files. **By default, minimum read coveage per base is 10.**

#### 3.2.1 CpG methylation
Scripts `processBismarkAln_subgenome_compare_CG_V1.r` and `processBismarkAln_subgenome_compare_CG_V1.sh` were used. 
#### 3.2.2 CHG methylation
Scripts `processBismarkAln_subgenome_compare_CHG_V1.r` and `processBismarkAln_subgenome_compare_CHG_V1.sh` were used.
#### 3.2.3 CHH methylation
Scripts `processBismarkAln_subgenome_compare_CHH_V1.r` and `processBismarkAln_subgenome_compare_CHH_V1.sh` were used.

### 3.3 Read Bismark alignment for T. dubius, T. pratensis, T. miscellus, and two subgenomes; mincov = 5
The minimum read coverage to call a mehtylation status for a base is 5 (mincov = 5)
#### 3.3.1 CpG methylation
For T. dubius, T. pratensis, and T. miscellus, scripts `processBismarkAln_CG_minCOV-5_V1.r` and `processBismarkAln_CG_minCOV-5_V1.sh` were used. For the subgneomes of T. miscellus, scripts `processBismarkAln_subgenome_compare_CG_minCOV-5_V1.r` and `processBismarkAln_subgenome_compare_CG_minCOV-5_V1.sh` were used.

Example output files:
```
Tdu_1_cov5_CpG.txt
Tdu_2_cov5_CpG.txt
Tms_1_du_cov5_CpG.txt
Tms_2_du_cov5_CpG.txt
```
#### 3.3.2 CpG methylation; identify overlapping bases cross diploid parents and subgenomes
The output from this step is used for downstream DMR comparing analysis. The output is located at folder `/blue/soltis/shan158538/Methylation/OutPut/Overlap_bases/CpG_mincov_5`. The script that I used is: `Overlap_bases.V1.py Tms_1_du_cov5_CpG.txt Tms_2_du_cov5_CpG.txt Tms_1_pr_cov5_CpG.txt Tms_2_pr_cov5_CpG.txt Tdu_1_cov5_CpG.txt Tdu_2_cov5_CpG.txt Tpr_1_cov5_CpG.txt Tpr_2_cov5_CpG.txt &`. **In total, there are 1,083,850 shared bases.**

Example outputs:
```
Tdu_1_cov5_CpG_overlap.txt
Tdu_2_cov5_CpG_overlap.txt
Tms_1_du_cov5_CpG_overlap.txt
```

### 3.4 Read Bismark alignment for T. dubius, T. pratensis, T. miscellus, and two subgenomes; mincov = 1
The minimum read coverage to call a mehtylation status for a base is 1 (mincov = 1). **The purpose of this step is to collect all base information at the very beginning of the analysis, and then may apply the filter for minimum coverage in following steps. However, after identify overlapping bases across the diploid parents and the two subgenomes, if then apply mincov = 3 when using function methRead, the bases left won't be found in all samples. Therefore, it seems like the mincov = 3 filter should be applied when calling methylation from Bismark alignment (check out section 2.5).**
#### 3.4.1 CpG methylation
For T. dubius, T. pratensis, and T. miscellus, scripts `processBismarkAln_CG_minCOV-1_V1.r` and `processBismarkAln_CG_minCOV-1_V1.sh` were used. For the subgneomes of T. miscellus, scripts `processBismarkAln_subgenome_compare_CG_minCOV-1_V1.r` and `processBismarkAln_subgenome_compare_CG_minCOV-1_V1.sh` were used.

Example output files:
```
Tdu_1_cov1_CpG.txt
Tdu_2_cov1_CpG.txt
Tms_1_du_cov1_CpG.txt
Tms_2_du_cov1_CpG.txt
```
