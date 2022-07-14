# Permutation analysis of DMR between Tdu and Tpr (CpG context)
When analyzing DMR between Tdu and Tpr, we found higher number of hyper-methylation in Tdu than that in Tpr. Is this real? We, therefore, would like to perform permutation analysis (by swithcing lables when running the scripts) to differentiate biological differences and technical differences.

Scripts `MethylDiff_CpG_Tdu-Tpr_permutation_V1.r` and `MethylDiff_CpG_Tdu-Tpr_permutation_V1.sh` were used.

Results:
| | Tdu_1 | Tdu_2 | Tpr_1 | Tpr_2 | Diff_all | Hyper (1>0) | Hypo (0>1) |
|-|-|-|-|-|-|-|-|
| Orignial | 1 | 1 | 0 | 0 | 10,558| 8,809| 1,749|
| Permutation_1 |0 | 0|1 |1 |10,558 | 1,749| 8,809|
| Permutation_2 | 0| 1| 0| 1|1,508 | 387| 1,121|
| Permutation_3 | 1| 0| 1| 0| 1,508| 1,121| 387|
| Permutation_4 | 0| 1| 1| 0| 1,901| 288| 1,613|
| Permutation_5 | 1| 0| 0| 1| 1,901| 1,613| 288|
