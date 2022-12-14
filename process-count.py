#!/usr/bin/env python
import glob
import argparse
import numpy as np
import pandas as pd
# Parse arguments
parser = argparse.ArgumentParser()
parser.add_argument('-f', help='Gene Length File')
args = parser.parse_args()

# Create generator to read the htseq files
files = sorted(glob.glob("*.bam_htseq.txt"))
df_from_each_file = (pd.read_csv(f,sep='\t',names=['GeneID',f.replace('.bam_htseq.txt', '')],index_col=0) for f in files)

# Concatenate the files into a single dataframe
df = pd.concat(df_from_each_file, axis=1)
n=5
df.drop(df.tail(n).index,inplace=True)

# Save raw counts
df.to_csv('RawCount.csv')

# Read gene lengths
gAnnot = pd.read_csv(args.f,delimiter='\t',names=['GeneID','Len'],index_col=0)

# Normalize to TPM
if (len(df)!=len(gAnnot)):
    raise('Gene IDs not proper')
df = df.divide(gAnnot['Len'],axis='rows')
df = df.divide(df.sum(),axis='columns') * 1e6
df = np.log2(df+1)

# Save TPM counts
df.to_csv('TPMCounts.csv')