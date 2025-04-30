import os
import pandas as pd
from collections import Counter
from tqdm import tqdm
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import numpy as np
pd.set_option('display.max_columns', None)

workspace = '/home/adrian/Data/TRR319_RMaP_BaseCalling_RNA004/Adrian/DTU/output'
file_ctrl = 'CTRL_01_transcripts_table.tsv'
file_kd = 'TRUB1_KD_01_transcripts_table.tsv'

img_out = '/home/adrian/img_out/DTU_ctrl_kd'
os.makedirs(img_out, exist_ok=True)

df_ctrl = pd.read_csv((os.path.join(workspace, file_ctrl)), sep='\t')
df_ctrl = df_ctrl[(df_ctrl['cov'] > 0) * (df_ctrl['class_code'] == '=')]
valid_genes_ctrl = [gene for gene, num_transcripts in Counter(df_ctrl['ref_gene_id']).items() if num_transcripts > 1]

df_kd = pd.read_csv((os.path.join(workspace, file_kd)), sep='\t')
df_kd = df_kd[(df_kd['cov'] > 0) * (df_kd['class_code'] == '=')]
valid_genes_kd = [gene for gene, num_transcripts in Counter(df_kd['ref_gene_id']).items() if num_transcripts > 1]

common_valid_genes = list(set(valid_genes_ctrl).intersection(set(valid_genes_kd)))

keep_fields = [
    'ref_gene_id',
    'ref_id',
    'class_code',
    'num_exons',
    'cov',
    'len',
    'sample_id',
    'parent gene iso num'
]

dfs_merged = []
for this_gene in tqdm(common_valid_genes):
    sub_df_ctrl = df_ctrl[df_ctrl['ref_gene_id'] == this_gene][keep_fields]
    sub_df_kd = df_kd[df_kd['ref_gene_id'] == this_gene][keep_fields]
    sub_df_merged = pd.merge(sub_df_ctrl, sub_df_kd, on=['ref_gene_id', 'ref_id'], suffixes=['_ctrl', '_kd'])
    sub_df_merged['total_ctrl'] = sub_df_merged['cov_ctrl'].sum()
    sub_df_merged['ratio_ctrl'] = sub_df_merged['cov_ctrl'] / sub_df_merged['total_ctrl'] * 100.0
    sub_df_merged['total_kd'] = sub_df_merged['cov_kd'].sum()
    sub_df_merged['ratio_kd'] = sub_df_merged['cov_kd'] / sub_df_merged['cov_kd'].sum() * 100.0
    sub_df_merged['delta_ratio'] = sub_df_merged['ratio_kd'] - sub_df_merged['ratio_ctrl']
    sub_df_merged['max_abs_delta_ratio'] = sub_df_merged['delta_ratio'].abs().max()
    dfs_merged.append(sub_df_merged)
df_merged = pd.concat(dfs_merged)

thresh_total = 20
thresh_delta = 50.0

df_merged_filtered = df_merged[
    (df_merged['total_ctrl'] >= thresh_total)
    * (df_merged['total_kd'] >= thresh_total)
]

plt.figure(figsize=(5, 5))
plt.hist(df_merged_filtered['max_abs_delta_ratio'], range=[0, 100], bins=20, log=True)
plt.xlabel('Max. Abs. $\Delta$ isofrom %')
plt.ylabel('Counts')
plt.title(f'Gene coverage $\geq$ {thresh_total}')
plt.savefig(os.path.join(img_out, 'hist_delta_ratio.png'), bbox_inches='tight')

df_sel = df_merged_filtered[df_merged_filtered['max_abs_delta_ratio'] >= thresh_delta]
df_sel.sort_values('ref_gene_id', inplace=True)
df_sel.to_csv(os.path.join(workspace, f'transcript_table_gene_count_{thresh_total}_delta_{thresh_delta}.tsv'),
              sep='\t', index=False, float_format='%.3f')

cmap = matplotlib.colormaps['Spectral']
for this_gene in df_sel['ref_gene_id'].unique():
    sub_df_sel = df_sel[df_sel['ref_gene_id'] == this_gene].copy()
    sub_df_sel.sort_values('ref_id', inplace=True)
    transcripts = sub_df_sel['ref_id']
    num_isoforms = len(transcripts)
    shifts = np.linspace(-0.15, 0.15, num_isoforms)
    bar_width = (0.5 / num_isoforms) * 0.75
    tx_colors = cmap(np.linspace(0, 1, num_isoforms))

    plt.figure(figsize=(5, 5))
    for tx_ind, this_transcript in enumerate(transcripts):
        vec_x = np.array([1, 2]) + shifts[tx_ind]
        vec_y = sub_df_sel[sub_df_sel['ref_id'] == this_transcript][['ratio_ctrl', 'ratio_kd']].to_numpy()[0]
        plt.bar(vec_x, vec_y, width=bar_width, color=tx_colors[tx_ind], label=this_transcript)
    plt.xticks([1, 2], [f"CTRL ({sub_df_sel['total_ctrl'][0]})", f"TRUB1-KD ({sub_df_sel['total_kd'][0]})"])
    plt.xlabel('Conditions (total num. reads)')
    plt.ylabel('Isoform ratio (%)')
    plt.legend(loc='upper center')
    plt.title(this_gene)
    plt.savefig(os.path.join(img_out, f'{this_gene}_isoforms.png'), bbox_inches='tight')

