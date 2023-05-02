#!/usr/bin/env python3

# The script is used to calculate cytosine coverage of sample 1 combining CG, CHG, and CHH contexts
Usage = """
Cytosine_coverage_V1.py
Usage:
	Cytosine_coverage_V1.py sample1_CG_bismark.cov sample1_CHG_bismark.cov sample1_CHH_bismark.cov
"""

import sys,os,re
import pandas as pd

# Input files
InFileName1 = sys.argv[1]
InFileName2 = sys.argv[2]
InFileName3 = sys.argv[3]

Delimiter = '\t'

# header=None: the first row fo the CSV file will be treated as data rather than column
InFile1 = pd.read_csv(InFileName1, sep=Delimiter, header=None, names=['chromosome', 'Start', 'End', 'methylation_p', 'count_methylated', 'count_unmethylated'])
InFile2 = pd.read_csv(InFileName2, sep=Delimiter, header=None, names=['chromosome', 'Start', 'End', 'methylation_p', 'count_methylated', 'count_unmethylated'])
InFile3 = pd.read_csv(InFileName3, sep=Delimiter, header=None, names=['chromosome', 'Start', 'End', 'methylation_p', 'count_methylated', 'count_unmethylated'])


# Total number of reads with methylated cytocines
sum_methylated_count = InFile1["count_methylated"].sum() + InFile2["count_methylated"].sum() + InFile3["count_methylated"].sum()

# Total number of reads with unmethlated cytocines
sum_unmethylated_count = InFile1["count_unmethylated"].sum() + InFile2["count_unmethylated"].sum() + InFile3["count_unmethylated"].sum()

# Total number of sites
sum_sites = len(InFile1) + len(InFile2) + len(InFile3)

# Coverage
Coverage = (sum_methylated_count + sum_unmethylated_count)/sum_sites
print("Coverage is:", Coverage)