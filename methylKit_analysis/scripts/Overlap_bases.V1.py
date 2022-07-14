#!/usr/bin/env python3

import sys,os,re
import pandas as pd

InFileName1 = 'Tms_1_du_cov5_CpG.txt'
InFileName2 = 'Tms_2_du_cov5_CpG.txt'
InFileName3 = 'Tms_1_pr_cov5_CpG.txt'
InFileName4 = 'Tms_2_pr_cov5_CpG.txt'
InFileName5 = 'Tdu_1_cov5_CpG.txt'
InFileName6 = 'Tdu_2_cov5_CpG.txt'
InFileName7 = 'Tpr_1_cov5_CpG.txt'
InFileName8 = 'Tpr_2_cov5_CpG.txt'

OutFileName1 = os.path.splitext(InFileName1)[0] + '_overlap.txt'
OutFileName2 = os.path.splitext(InFileName2)[0] + '_overlap.txt'
OutFileName3 = os.path.splitext(InFileName3)[0] + '_overlap.txt'
OutFileName4 = os.path.splitext(InFileName4)[0] + '_overlap.txt'
OutFileName5 = os.path.splitext(InFileName5)[0] + '_overlap.txt'
OutFileName6 = os.path.splitext(InFileName6)[0] + '_overlap.txt'
OutFileName7 = os.path.splitext(InFileName7)[0] + '_overlap.txt'
OutFileName8 = os.path.splitext(InFileName8)[0] + '_overlap.txt'

OutFile1 = open(OutFileName1, 'w')
OutFile2 = open(OutFileName2, 'w')
OutFile3 = open(OutFileName3, 'w')
OutFile4 = open(OutFileName4, 'w')
OutFile5 = open(OutFileName5, 'w')
OutFile6 = open(OutFileName6, 'w')
OutFile7 = open(OutFileName7, 'w')
OutFile8 = open(OutFileName8, 'w')
Delimiter = '\t'

InFile1 = pd.read_csv(InFileName1, sep=Delimiter, header=0)
InFile2 = pd.read_csv(InFileName2, sep=Delimiter, header=0)
InFile3 = pd.read_csv(InFileName3, sep=Delimiter, header=0)
InFile4 = pd.read_csv(InFileName4, sep=Delimiter, header=0)
InFile5 = pd.read_csv(InFileName5, sep=Delimiter, header=0)
InFile6 = pd.read_csv(InFileName6, sep=Delimiter, header=0)
InFile7 = pd.read_csv(InFileName7, sep=Delimiter, header=0)
InFile8 = pd.read_csv(InFileName8, sep=Delimiter, header=0)

Overlap_1_2 = pd.merge(InFile1, InFile2, on=['chrBase','strand'], how='inner')
Overlap_1_2_3 = pd.merge(Overlap_1_2, InFile3, on=['chrBase','strand'], how='inner')
Overlap_1_2_3_4 = pd.merge(Overlap_1_2_3, InFile4, on=['chrBase','strand'], how='inner')
Overlap_1_2_3_4_5 = pd.merge(Overlap_1_2_3_4, InFile5, on=['chrBase','strand'], how='inner')
Overlap_1_2_3_4_5_6 = pd.merge(Overlap_1_2_3_4_5, InFile6, on=['chrBase','strand'], how='inner')
Overlap_1_2_3_4_5_6_7 = pd.merge(Overlap_1_2_3_4_5_6, InFile7, on=['chrBase','strand'], how='inner')
Overlap_all = pd.merge(Overlap_1_2_3_4_5_6_7, InFile8, on=['chrBase','strand'], how='inner')

InFile1_filtered = InFile1[InFile1.chrBase.isin(Overlap_all.chrBase)]
InFile2_filtered = InFile2[InFile2.chrBase.isin(Overlap_all.chrBase)]
InFile3_filtered = InFile3[InFile3.chrBase.isin(Overlap_all.chrBase)]
InFile4_filtered = InFile4[InFile4.chrBase.isin(Overlap_all.chrBase)]
InFile5_filtered = InFile5[InFile5.chrBase.isin(Overlap_all.chrBase)]
InFile6_filtered = InFile6[InFile6.chrBase.isin(Overlap_all.chrBase)]
InFile7_filtered = InFile7[InFile7.chrBase.isin(Overlap_all.chrBase)]
InFile8_filtered = InFile8[InFile8.chrBase.isin(Overlap_all.chrBase)]

InFile1_filtered.to_csv(OutFile1, sep=Delimiter, index=False)
InFile2_filtered.to_csv(OutFile2, sep=Delimiter, index=False)
InFile3_filtered.to_csv(OutFile3, sep=Delimiter, index=False)
InFile4_filtered.to_csv(OutFile4, sep=Delimiter, index=False)
InFile5_filtered.to_csv(OutFile5, sep=Delimiter, index=False)
InFile6_filtered.to_csv(OutFile6, sep=Delimiter, index=False)
InFile7_filtered.to_csv(OutFile7, sep=Delimiter, index=False)
InFile8_filtered.to_csv(OutFile8, sep=Delimiter, index=False)




