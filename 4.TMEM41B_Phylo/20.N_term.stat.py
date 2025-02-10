# -*- coding: utf-8 -*-
"""
analysis of N-terminal length and similairty from TMEM41B pairwise alignment

Created on Thu Jan 16 09:23:12 2025

@author: xuyu
"""
import os
import pandas as pd
import numpy as np
from itertools import chain
import matplotlib.pyplot as plt
human = "homo_sapiens"

# p5 stands for 5-prime end, or N-terminal in protein

#%% get align
df_merge = pd.read_csv("tables/seq.pair_align.csv")
homo_species = pd.read_csv("data/homology/homo.species.csv")
homo_species

df_align = df_merge[(df_merge['source.species'].isin(homo_species.species) & 
                     (df_merge['target.species'].isin(homo_species.species)))]

align_cols = ['species', 'protein_id', 'align_seq']
align_colnames = list(chain(*[("source." + col, "target." + col) for col in align_cols]))
df_align = df_align[align_colnames]
df_align
df_align_raw = df_align
df_merge

#%% p5 align functions
def find_kth_non_dash(s, k):
    if k == -1:
        return -1 
    count = 0  # To count non '-' characters
    for index, char in enumerate(s):
        if char != '-':
            count += 1
            if count == k:
                return index  # Return the index of the k-th non '-' character
    return -1  # Return -1 if k is greater than the number of non '-' characters

def count_non_dash_till_index(s, i):
    if i == -1:
        return -1 
    if i < 0 or i > len(s):
        raise ValueError("Index i is out of bounds.")
    
    count = 0  # To count non '-' characters
    for char in s[:i+1]:  # Iterate until index i
        if char != '-':
            count += 1
    return count

def calc_aligned_p5len(p5len, source_align, target_align):
    end_idx = find_kth_non_dash(source_align, p5len)
    align_len = count_non_dash_till_index(target_align, end_idx)
    return end_idx, align_len

#%% p5 len for human

#human_p5_len = 60
human_p5_len = 67
p5_len = homo_species
p5_len['p5len'] = -1
p5_len.loc[p5_len.species == human,'p5len'] = human_p5_len
sum(p5_len.species == human)
sum(p5_len.p5len > 0)

#%% round 1: get p5_len for other species from human
df_align = df_align_raw.merge(p5_len, left_on = 'source.species', right_on='species', how='left')
df_align[['p5_end', 'target.p5len']] = df_align.apply(
    lambda row: calc_aligned_p5len(row['p5len'], row['source.align_seq'], row['target.align_seq']), 
    axis=1, result_type='expand')
df_align
sum(df_align['target.p5len']>0)
aligned = df_align[df_align.p5len > 0]
#tmp[tmp['target.species']==human]['target.p5len']
aligned_p5len = aligned[['target.species', 'target.protein_id', 'target.p5len']]
aligned_p5len.columns = aligned_p5len.columns.str.replace('target.', '', regex=False)
aligned_p5len = aligned_p5len.drop_duplicates()

aligned_p5len.species.value_counts().value_counts()
aligned_count = aligned_p5len.species.value_counts()
aligned_multi = aligned_count[aligned_count > 1]
aligned_p5len[aligned_p5len['species'].isin(aligned_multi.index)]
plt.hist(aligned_p5len.p5len)

p5_len = p5_len.merge(aligned_p5len, on='species', how='left')
any(p5_len['species'] == human)
p5_len['p5len'] = p5_len.apply(lambda row: row['p5len_y'] if row['p5len_x'] == -1 else row['p5len_x'], axis=1)
sum(p5_len.p5len < 0)
p5_len = p5_len.drop(columns=['p5len_x', 'p5len_y'])


#%% round 2: get 5p end by cross-species pairalign
p5_len_source = p5_len.add_prefix("source.")
align =  df_align_raw.merge(p5_len_source, on = ['source.species', 'source.protein_id'], how='left')
align[['source.p5_end', 'source.target.p5len']] = align.apply(
    lambda row: calc_aligned_p5len(row['source.p5len'], row['source.align_seq'], row['target.align_seq']), 
    axis=1, result_type='expand')

p5_len_target = p5_len.add_prefix("target.")
align = align.merge(p5_len_target, on = ['target.species', 'target.protein_id'], how='left')
align[['target.p5_end', 'target.source.p5len']] = align.apply(
    lambda row: calc_aligned_p5len(row['target.p5len'], row['target.align_seq'], row['source.align_seq'] ), 
    axis=1, result_type='expand')
#check consistence
sum(align['source.p5_end']==align['target.p5_end'])
align['p5end_delta'] =  abs(align['source.p5_end'] - align['target.p5_end'])
plt.hist(abs(align['source.p5_end'] - align['target.p5_end']))
align['p5_end'] = np.maximum(align['source.p5_end'], align['target.p5_end'])



#%% check result consistency from round 1 and round 2
align.columns
p5len_v1_a = align[['source.species', 'source.protein_id', 'source.p5len']]
p5len_v1_a.columns = p5len_v1_a.columns.str.replace('source.', '')
p5len_v1_b = align[['target.species', 'target.protein_id', 'target.p5len']]
p5len_v1_b.columns = p5len_v1_b.columns.str.replace('target.', '')
p5len_v1 = pd.concat([p5len_v1_a, p5len_v1_b])
p5len_v1 = p5len_v1.drop_duplicates()
p5len_v1.shape
plt.hist(p5len_v1.p5len)
sum(np.isnan(p5len_v1.p5len))
np.nanquantile(p5len_v1.p5len, np.linspace(0, 1, num=11))

p5len_v2_a = align[['source.species', 'source.protein_id', 'target.source.p5len']]
p5len_v2_a.columns = p5len_v2_a.columns.str.replace('target.', '').str.replace('source.', '')
p5len_v2_b = align[['target.species', 'target.protein_id', 'source.target.p5len']]
p5len_v2_b.columns = p5len_v2_b.columns.str.replace('source.', '').str.replace('target.', '')
p5len_v2 = pd.concat([p5len_v2_a, p5len_v2_b])
p5len_v2.shape

p5len_v2_q = p5len_v2.groupby(['species', 'protein_id'])['p5len'].agg(
    **{f'p5len.q{int(q*100)}': lambda x, q=q: np.quantile(x, q) for q in np.array([0,10,25,50,75,90,100])/100}
).reset_index()

p5len_merge = p5len_v2_q.merge(p5len_v1.rename(columns={'p5len': 'p5len.by_marker'}), on=['species', 'protein_id'], how='left')
(p5len_merge['p5len.q75'] - p5len_merge['p5len.q25']).value_counts()
(p5len_merge['p5len.q90'] - p5len_merge['p5len.q10']).value_counts()
(p5len_merge['p5len.by_marker'] - p5len_merge['p5len.q50']).value_counts()
sum(np.isnan(p5len_merge['p5len.q50']))
p5len_result = p5len_merge[['species', 'protein_id', 'p5len.q50']].rename(columns={'p5len.q50':'p5len'})
p5len_result.shape
#p5len_result.to_csv("tables/seq.pair_align.dm67.p5len.csv", index=False)
p5len_result.to_csv("tables/seq.pair_align.hs67.p5len.csv", index=False)


#%% distance function
def normalized_hamming(p5_end, align_a, align_b):
    d = [int(char1 != char2) for char1, char2 in zip(align_a, align_b)]
    p5_len = p5_end + 1
    if p5_len > 0:
        p5_d = sum(d[0:(p5_end+1)]) / p5_len
    else:
        p5_d = 1
    p3_len = len(align_a) - p5_len
    if p3_len > 0:
        p3_d = sum(d[(p5_end+1):]) / p3_len
    else:
        p3_d = 1
    return p5_d, p3_d

#%% calculate distance
align[['p5d', 'p3d']] = align.apply(
    lambda row: normalized_hamming(row['p5_end'], row['source.align_seq'], row['target.align_seq']),
    axis=1, result_type='expand')

align.columns
align_out_cols = ['source.species', 'target.species', 'source.protein_id','target.protein_id', 'p5_end','p5d', 'p3d' ]
align_out = align[align_out_cols]
align_out.shape
align.to_csv("tables/seq.pair_align.hs67.p5sim.ext.csv", index=False)
align_out.to_csv("tables/seq.pair_align.hs67.p5sim.csv", index=False)





