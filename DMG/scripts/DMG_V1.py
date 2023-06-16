#!/usr/bin/env python3

Usage = """
This code is used to identify genes that are overlapped with DMRs
DMG_V1.py
Usage:
	DMG_V1.py Tdub.V1.rm.gff DMR.txt
"""

import sys,os,re
import pandas as pd

Delimiter = '\t'

# Input files
# InFileName1 should be the gff file; InFileName2 should be the DMR file
InFileName1 = sys.argv[1]
InFileName2 = sys.argv[2]

# Output file 1: identify DMRs that overlap with gene space, that is DMG
OutFileName1 = 'DMG_' + os.path.splitext(InFileName2)[0] + '.txt'
OutFile1 = open(OutFileName1, 'w')

if len(sys.argv) < 3:
    print(Usage)
    sys.exit(1)

# Read the gff file and add a header line
GFF = pd.read_csv(InFileName1, sep=Delimiter, header=None, names=['Scaffold', 'Source', 'Feature', 'Start', 'End', 'Score', 'Strand', 'Frame', 'Attribute'])
# only keep lines with a mRNA (i.e., gene) in the Feature column
GFF_gene = GFF[GFF['Feature'] == 'mRNA']

# Read the DMR file, which already contains the header
DMR = pd.read_csv(InFileName2, sep=Delimiter)
# Only keep two columns, Scaffold and Start, from the df
DMR_filtered = DMR[['Scaffold', 'Start']]

# Merge the two dfs, and only keep those lines with shared Scaffold
merged_df = pd.merge(GFF_gene, DMR_filtered, on="Scaffold")
# Only keep those genes that are overalpped with DMRs; Start_y is the start position of a DMR; Start_x is the start position
# of a gene; End is the end position of a gene; the length of each DMR is 300 bp
filtered_merged_df = merged_df[(merged_df['Start_x']-300 < merged_df['Start_y']) & (merged_df['Start_y'] <= merged_df['End'])]

# Extract the gene ID in the Attribute column; add a new column (Gene_ID)
filtered_merged_df['Gene_ID'] = filtered_merged_df['Attribute'].str.extract(r'ID=(.*?);')
# Extract the Gene_ID column and exlude the redundant IDs
Unique_gene_id = filtered_merged_df['Gene_ID'].drop_duplicates()

# Save the gene ids
Unique_gene_id.to_csv(OutFile1, sep=Delimiter, header=False, index=False)

