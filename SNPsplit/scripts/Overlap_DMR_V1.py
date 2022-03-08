#!/usr/bin/env python3

Usage = """
Overlap_DMR_V1.py
Usage:
	Overlap_DMR_V1.py DMR_1.txt DMR_2.txt
"""

import sys,os,re
import pandas as pd

# Input files
InFileName1 = sys.argv[1]
InFileName2 = sys.argv[2]

# Output file 1: overlaping DMR with the same direction of change
OutFileName1 = 'Overlap_same_direction_' + os.path.splitext(InFileName1)[0] + '_' + os.path.splitext(InFileName2)[0] + '.txt'
# Output file 2: overlaping DMR with different direction of change
OutFileName2 = 'Overlap_diff_direction_' + os.path.splitext(InFileName1)[0] + '_' + os.path.splitext(InFileName2)[0] + '.txt'
# Output file 3: DMR unique to input file 1
OutFileName3 = 'Unique_' + os.path.splitext(InFileName1)[0] + '.txt'
# Output file 4: DMR unique to input file 2
OutFileName4 = 'Unique_' + os.path.splitext(InFileName2)[0] + '.txt'

OutFile1 = open(OutFileName1, 'w')
OutFile2 = open(OutFileName2, 'w')
OutFile3 = open(OutFileName3, 'w')
OutFile4 = open(OutFileName4, 'w')
Delimiter = '\t'

if len(sys.argv) < 3:
    print(Usage)

# skiprows=10: skip the first 11 rows
# names: list of column names to use
InFile1 = pd.read_csv(InFileName1, sep=Delimiter, skiprows=10, header=0, names=['Scaffold', 'Start', 'End', 'Strand', 'P', 'Q', 'Direction'])
InFile2 = pd.read_csv(InFileName2, sep=Delimiter, skiprows=10, header=0, names=['Scaffold', 'Start', 'End', 'Strand', 'P', 'Q', 'Direction'])


# Find the overlap DMR between InFile1 and InFile2 based on the scaffold ID and DMR start position
Overlap = pd.merge(InFile1, InFile2, on=['Scaffold', 'Start'], how='inner')
# Find the overlap DMR with the same direaction of change
OverlapSameDirection = Overlap[Overlap['Direction_x'] * Overlap['Direction_y'] > 0]
# Find the overlap DMR witht different direaction of change
OverlapDiffDirection = Overlap[Overlap['Direction_x'] * Overlap['Direction_y'] < 0]


# Find the full outer join between InFile1 and InFile2 based on the scaffold ID and DMR start position
Outer = pd.merge(InFile1, InFile2, on=['Scaffold', 'Start'], how='outer', indicator=True)
# Find the DMR only identified in InFile1 (left file)
UniqueInFile1 = Outer[Outer['_merge']=='left_only']
# Find the DMR only identified in InFile2 (right file)
UniqueInFile2 = Outer[Outer['_merge']=='right_only']


# Save the overlap DMR with the same direction of change
OverlapSameDirection.to_csv(OutFile1, sep=Delimiter)
# Save the overlap DMR with different direction of change
OverlapDiffDirection.to_csv(OutFile2, sep=Delimiter)
# Save the InFile1 only DMR
UniqueInFile1.to_csv(OutFile3, sep=Delimiter)
# Save the InFile2 only DMR
UniqueInFile2.to_csv(OutFile4, sep=Delimiter)


