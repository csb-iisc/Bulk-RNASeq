#!/usr/bin/env python
import argparse
import numpy as np
import pandas as pd
# Parse arguments
parser = argparse.ArgumentParser()
parser.add_argument('-f', help='Gene Length File')
args = parser.parse_args()

# Read gene lengths
gAnnot = pd.read_csv(args.f,delimiter='\t',names=['GeneID','Gene','Len'],index_col=0)
# Read raw counts
df=pd.read_csv('RawCounts_htseq.tsv',delimiter='\t',index_col=0)
# Remove last 5 lines
n=5
df.drop(df.tail(n).index,inplace=True)

# Normalize to TPM
if (len(df)!=len(gAnnot)):
    raise('Gene IDs not proper')
df = df.divide(gAnnot['Len'],axis='rows')
df = df.divide(df.sum(),axis='columns') * 1e6
df = np.log2(df+1)
# Add Gene Names
df = pd.merge(gAnnot['Gene'],df,left_index=True,right_index=True)
# Save TPM counts
df.to_csv('TPMCounts.tsv',sep='\t')
