#!/usr/bin/env python3

import sys,os,re
import pandas as pd

## Identify shared cytosins between six files, only include loci with coverage >= 3
## threshold is the coverage cutoff
## There are six input files, two replicate from each species: T. dubius, T. pratensis, and T. miscellus

Delimiter = '\t'
threshold = 3

InFileName1 = sys.argv[1]
InFileName2 = sys.argv[2]
InFileName3 = sys.argv[3]
InFileName4 = sys.argv[4]
InFileName5 = sys.argv[5]
InFileName6 = sys.argv[6]

OutFileName1 = os.path.splitext(InFileName1)[0] + '_shared_filtered.cov'
OutFileName2 = os.path.splitext(InFileName2)[0] + '_shared_filtered.cov'
OutFileName3 = os.path.splitext(InFileName3)[0] + '_shared_filtered.cov'
OutFileName4 = os.path.splitext(InFileName4)[0] + '_shared_filtered.cov'
OutFileName5 = os.path.splitext(InFileName5)[0] + '_shared_filtered.cov'
OutFileName6 = os.path.splitext(InFileName6)[0] + '_shared_filtered.cov'

## read the Bismark cov files
InFile1 = pd.read_csv(InFileName1, sep=Delimiter, header=None, names=['chromosome', 'Start', 'End', 'methylation_p', 'count_methylated', 'count_unmethylated'])
InFile2 = pd.read_csv(InFileName2, sep=Delimiter, header=None, names=['chromosome', 'Start', 'End', 'methylation_p', 'count_methylated', 'count_unmethylated'])
InFile3 = pd.read_csv(InFileName3, sep=Delimiter, header=None, names=['chromosome', 'Start', 'End', 'methylation_p', 'count_methylated', 'count_unmethylated'])
InFile4 = pd.read_csv(InFileName4, sep=Delimiter, header=None, names=['chromosome', 'Start', 'End', 'methylation_p', 'count_methylated', 'count_unmethylated'])
InFile5 = pd.read_csv(InFileName5, sep=Delimiter, header=None, names=['chromosome', 'Start', 'End', 'methylation_p', 'count_methylated', 'count_unmethylated'])
InFile6 = pd.read_csv(InFileName6, sep=Delimiter, header=None, names=['chromosome', 'Start', 'End', 'methylation_p', 'count_methylated', 'count_unmethylated'])

## only remain loci with coverage >= 3 in each file
filtered_InFile1 = InFile1[InFile1['count_methylated'] + InFile1['count_unmethylated'] >= threshold]
filtered_InFile2 = InFile2[InFile2['count_methylated'] + InFile2['count_unmethylated'] >= threshold]
filtered_InFile3 = InFile3[InFile3['count_methylated'] + InFile3['count_unmethylated'] >= threshold]
filtered_InFile4 = InFile4[InFile4['count_methylated'] + InFile4['count_unmethylated'] >= threshold]
filtered_InFile5 = InFile5[InFile5['count_methylated'] + InFile5['count_unmethylated'] >= threshold]
filtered_InFile6 = InFile6[InFile6['count_methylated'] + InFile6['count_unmethylated'] >= threshold]


## identify shared loci between six files
merged_file = pd.merge(filtered_InFile1, filtered_InFile2, on=['chromosome','Start'], how='inner')
merged_file = pd.merge(merged_file, filtered_InFile3, on=['chromosome','Start'], how='inner')
merged_file = pd.merge(merged_file, filtered_InFile4, on=['chromosome','Start'], how='inner')
merged_file = pd.merge(merged_file, filtered_InFile5, on=['chromosome','Start'], how='inner')
merged_file = pd.merge(merged_file, filtered_InFile6, on=['chromosome','Start'], how='inner')


## remain only shared loci with coverage >= 3 in each file
## merge each input file with the merged_file generated from the previous step
## only retain the first six columns which contain the information from each input file
## : before the comma represents all rows in the DataFrame, while :6 after the comma indicates the range of columns from 0 to 5 (excluding column 6)
Inner_1 = pd.merge(InFile1, merged_file, on=['chromosome', 'Start'])
shared_filtered_InFile1 = Inner_1.iloc[:, :6]

Inner_2 = pd.merge(InFile2, merged_file, on=['chromosome', 'Start'])
shared_filtered_InFile2 = Inner_2.iloc[:, :6]

Inner_3 = pd.merge(InFile3, merged_file, on=['chromosome', 'Start'])
shared_filtered_InFile3 = Inner_3.iloc[:, :6]

Inner_4 = pd.merge(InFile4, merged_file, on=['chromosome', 'Start'])
shared_filtered_InFile4 = Inner_4.iloc[:, :6]

Inner_5 = pd.merge(InFile5, merged_file, on=['chromosome', 'Start'])
shared_filtered_InFile5 = Inner_5.iloc[:, :6]

Inner_6 = pd.merge(InFile6, merged_file, on=['chromosome', 'Start'])
shared_filtered_InFile6 = Inner_6.iloc[:, :6]


## export the output files
shared_filtered_InFile1.to_csv(OutFileName1, sep=Delimiter, header=False, index=False)
shared_filtered_InFile2.to_csv(OutFileName2, sep=Delimiter, header=False, index=False)
shared_filtered_InFile3.to_csv(OutFileName3, sep=Delimiter, header=False, index=False)
shared_filtered_InFile4.to_csv(OutFileName4, sep=Delimiter, header=False, index=False)
shared_filtered_InFile5.to_csv(OutFileName5, sep=Delimiter, header=False, index=False)
shared_filtered_InFile6.to_csv(OutFileName6, sep=Delimiter, header=False, index=False)

