# DMR comparison
This section will perform downstream analysis after DMR identification -- find the overlapping DMRs. In order to easily edit/test scripts and uplode/downlode files (< 1-2 Gb), [Open OnDemand](https://help.rc.ufl.edu/doc/Open_OnDemand) (a web interface to starting and attaching to HiPerGator jobs) is used. Thanks for the help from Matt!

Script [`Overlap_DMR_V1.py`](https://github.com/GatorShan/Tragopogon-Methylation-Project/blob/master/SNPsplit/scripts/Overlap_DMR_V1.py) was used for following analyses.

## 1. DMR comparision -- between parents vs between subgenomes (CpG)
<img src="https://github.com/GatorShan/Tragopogon-Methylation-Project/blob/master/SNPsplit/images/DMR_CpG_between-parents_between-subgenomes.png" width="600"/>

| Category | Number |
| - | - |
| Overlaping DMR with same direction of change (Tdu vs Tpr) | 811 |
| Overlaping DMR with different direction of change | 44 |
| DMR identified only between two parents | 9,703 |
| DMR identified only between two subgenomes | 1,035 |
### 1.1 Question: 
For those DMRs only identified between parents but not between subgenomes (e.g., 9,703), is it because there is no SNP at these regions to differentiate the reads from the polyploid? Or, within those DMRs, there is read coverage only in the parents, but not in the subgenomes? By default, the minimum base coverage is 10. Maybe only include bases with coverage across all species and subgenomes? Summary of the total number of covered bases in all files:
Tms_1_du_CpG.txt	4,588,001

Tms_2_du_CpG.txt	2,772,112

Tms_1_pr_CpG.txt	2,146,950

Tms_2_pr_CpG.txt	1,469,788

Tdu_1_CpG.txt		25,606,465

Tdu_2_CpG.txt		34,126,949

Tpr_1_CpG.txt		5,128,222

Tpr_2_CpG.txt		5,463,740

Identified 209,025 overlapping bases across all files. Maybe try to use **mincov = 5**?



## 2. DMR comparision -- between parents vs between subgenomes (CHG) -- WRONG RESULTS

| Category | Number |
| - | - |
| Overlaping DMR with same direction of change (Tdu vs Tpr) | 1,048 |
| Overlaping DMR with different direction of change | 57 |
| DMR identified only between two parents | 9,309 |
| DMR identified only between two subgenomes | 849 |

## 3. DMR comparision -- between parents vs between subgenomes (CHH) -- WRONG RESULTS

| Category | Number |
| - | - |
| Overlaping DMR with same direction of change (Tdu vs Tpr) | 189 |
| Overlaping DMR with different direction of change | 9 |
| DMR identified only between two parents | 4,241 |
| DMR identified only between two subgenomes | 557 |
