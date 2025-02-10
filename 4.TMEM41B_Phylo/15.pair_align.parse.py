# -*- coding: utf-8 -*-
"""
Created on Thu Jan 16 09:23:12 2025

@author: xuyu
"""
import os
import pandas as pd
import numpy as np
from itertools import chain
import matplotlib.pyplot as plt
human = "homo_sapiens"

#%% merge data
data_dir="data/sequence_identity/pair_align.csv"
df_list = list()
for filename in os.listdir(data_dir):
    if filename.endswith('.csv'):
        print(filename)
        file_path = os.path.join(data_dir, filename)
        df = pd.read_csv(file_path)
        df_list.append(df)
df_merge = pd.concat(df_list)
df_merge.to_csv("data/sequence_identity/pair_align.merge.csv", index=False)
df_merge.to_csv("tables/seq.pair_align.csv", index=False)
colnames = df_merge.columns
eg = df_merge.head(10)

#%% exploration

# df_merge.type.value_counts()
# df_merge.dn_ds.value_counts()
# df_merge.method_link_type.value_counts()
# df_merge.taxonomy_level.value_counts()

#%% get all seqs
seq_cols = ['target.species', 'target.id', 'target.protein_id', 'target.align_seq']
#df_seq = df_merge[df_merge['source.species']=="homo_sapiens"][seq_cols]
df_seq = df_merge[seq_cols]
df_seq.columns = df_seq.columns.str.replace('target.', '', regex=False)
df_seq
df_seq['seq'] = df_seq['align_seq'].apply(lambda x: x.replace('-', ''))
df_seq
df_seq_uniq = df_seq.drop('align_seq', axis=1).drop_duplicates(keep='first')
df_seq_uniq.shape
df_seq_uniq.species.value_counts()
df_seq_uniq = df_seq_uniq.rename(columns={'id':'gene_id'}).sort_values(by=['species', 'gene_id', 'protein_id'])
df_seq_uniq.to_csv("tables/seq.pair_align.seqs.csv", index=False)

#%% seqs to fa
def seq2fa(seq_file, out_file):
    seq_df = pd.read_csv(seq_file)
    with open(out_file, 'w') as ofh:
        for rec in seq_df.itertuples():
            name= "|".join((rec.protein_id, rec.species))
            ofh.write(f">{name}\n")
            ofh.write(f"{rec.seq}\n")
seq2fa("tables/seq.pair_align.seqs.csv", "tables/seq.pair_align.seqs.fa")

        
