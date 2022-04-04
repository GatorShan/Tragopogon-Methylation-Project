# Permutation analysis of DMR between Tdu and Tpr (CpG context)
When analyzing DMR between Tdu and Tpr, we found higher number of hyper-methylation in Tdu than that in Tpr. Is this real? We, therefore, would like to perform permutation analysis (by swithcing lables when running the scripts) to differentiate biological differences and technical differences.

| | Tdu_1 | Tdu_2 | Tpr_1 | Tpr_2 |
|-|-|-|-|-|
| Orignial | 1 | 1 | 0 | 0 |
| Permutation_1 |0 | 0|1 |1 |
| Permutation_2 | 0| 1| 0| 1|
| Permutation_3 | 0| 1| 1| 0|
| Permutation_4 | 1| 0| 0| 1|
| Permutation_5 | 1| 0| 1| 0|

Script `MethylDiff_CpG_Tdu-Tpr_permutation_V1.r` and `MethylDiff_CpG_Tdu-Tpr_permutation_V1.sh` were used.
