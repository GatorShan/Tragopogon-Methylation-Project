#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd
import os
import sys


# In[2]:


## the code within the 'if' block will be executed only when the code runs directly
if __name__ == "__main__":
    if len(sys.argv) < 3 :
        print ("Usage: Reformat_bismark_cov_V3.py <sample_mType.gz.bismark.cov> <-o=outputDirectory>")
    else:
        inPut = sys.argv[1]
        path = sys.argv[2][3:]
        name = inPut.split('_')[0]
        mType = inPut.split('_')[1].split('.')[0]
        covFile = pd.read_table(inPut, sep="\t", header = None)
        ## 0: <chromosome>; 1: <start position>; 2: <end position>; 3 <methylation percentage>;
        ## 4: <count methylated>; 5: <count unmethylated>
        
        ## Select only those scaffold with gff annotations
        gffFile = pd.read_table("/orange/soltis/Tdub_Genome_V1/Tdub.V1.rm.gff", sep="\t", header = None)
        scaffold_range = gffFile[0].unique()
        scaffold_range = scaffold_range.tolist()
        covFile = covFile[covFile[0].isin(scaffold_range)]
        
        ## Add column 4 and 5 and get column 6, which contains the number of total reads
        covFile[6] = covFile[4] + covFile[5]
        ## Drop column 2, 3, 5; axis=1 means cloumn
        covFile = covFile.drop([2,3,5], axis=1)
        ## following the format from allC file
        covFile[2] = 'NA'
        covFile[7] = 'NA'
        
        if mType == 'CpG':
            covFile[3] = 'CG'
        else:
            covFile[3] = mType
            
        covFile.sort_index(axis=1, inplace=True)
        
        ## Extract all the unique items from the column you want to split upon
        ## Convert the array to list
        data_category_range = covFile[0].unique()
        data_category_range = data_category_range.tolist()
        
        ## writeout
        for scaffold in data_category_range:
            covFile[covFile[0] == scaffold].to_csv(str(path) + 'allc_' + str(name) + '_' + str(scaffold) + '.tsv', header=False, index=False, sep="\t")


# In[ ]:




