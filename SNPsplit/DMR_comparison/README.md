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

## 2. DMR comparision -- between parents vs between subgenomes (CHG)

| Category | Number |
| - | - |
| Overlaping DMR with same direction of change (Tdu vs Tpr) | 1,048 |
| Overlaping DMR with different direction of change | 57 |
| DMR identified only between two parents | 9,309 |
| DMR identified only between two subgenomes | 849 |

## 3. DMR comparision -- between parents vs between subgenomes (CHH)

| Category | Number |
| - | - |
| Overlaping DMR with same direction of change (Tdu vs Tpr) | 189 |
| Overlaping DMR with different direction of change | 9 |
| DMR identified only between two parents | 4,241 |
| DMR identified only between two subgenomes | 557 |
