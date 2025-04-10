from Bio import SeqIO
from tqdm import tqdm
import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.use('TkAgg')
#######################################################################
cm = 1/2.54  # centimeters in inches
gr = 1.618
dpi = 1200
mpl.rcParams['figure.dpi'] = dpi
mpl.rcParams['savefig.dpi'] = dpi
mpl.rcParams['font.size'] = 8
mpl.rcParams['legend.fontsize'] = 6
mpl.rcParams['xtick.labelsize'] = 8
mpl.rcParams['ytick.labelsize'] = 8
mpl.rcParams['xtick.major.size'] = 4
mpl.rcParams['ytick.major.size'] = 4
mpl.rcParams['lines.linewidth'] = 1
mpl.rcParams['font.family'] = 'Arial'
# FMT = 'svg'
# fig_kwargs = dict(format=FMT, bbox_inches='tight', dpi=dpi, transparent=True)
FMT = 'png'
fig_kwargs = dict(format=FMT, bbox_inches='tight', dpi=dpi)
#######################################################################
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import os
import pickle
from functools import reduce


REF_FILE = '/home/adrian/Data/genomes/homo_sapiens/GRCh38_102/GRCh38_102.fa'
IMG_OUT = '/home/adrian/img_out/RNA004_psi_KD_OE_analysis'
BASE_DIR = '/home/adrian/Data/TRR319_RMaP_BaseCalling_RNA004/Adrian'
WRITER = 'TRUB1'

CONDITIONS = [f'{WRITER}-KD', 'CTRL']
# CONDITIONS = ['OE', 'WT']

def get_merged_df_with_delta(in_dfs):
    merged_fields = [
        'chrom',
        'chromStart',
        'chromEnd',
        'name',
        'strand',
        'ref_motif'
    ]

    this_mod_dfs = []
    for cond in CONDITIONS:
        df = in_dfs[cond].copy()
        this_mod_dfs.append(
            df.rename(columns={'score': f'score_{cond}', 'frequency': f'freq_{cond}'})
        )
    out_df_merged = reduce(lambda left, right: pd.merge(left, right, on=merged_fields, how='inner'), this_mod_dfs)

    mask = ((out_df_merged.loc[:, out_df_merged.columns.str.contains('freq_')] > 0).any(axis=1))
    out_df_merged = out_df_merged[mask]

    if f'{writer}-KD' in CONDITIONS:
        out_df_merged[f'delta_{writer}-KD'] = out_df_merged[f'freq_{writer}-KD'] - out_df_merged['freq_CTRL']
    elif f'{writer}-OE' in CONDITIONS:
        out_df_merged[f'delta_{writer}-OE'] = out_df_merged[f'freq_{writer}-OE'] - out_df_merged['freq_CTRL']

    return out_df_merged


img_out = IMG_OUT
os.makedirs(img_out, exist_ok=True)

base_dir = BASE_DIR
writer = WRITER

mod_names = ['17802', 'a']
dict_mod_display = {
    'a': 'm^6A',
    '17802': '\psi'
}

mod_motif_name = {
    '17802': writer,
    'a': 'DRACH'
}

cond_colors = {
    'WT': 'gray',
    'OE': 'red',
    'KD': 'blue'
}

### merged and plot delta ###
with open(os.path.join(base_dir, f'dfs_mod_filtered_{writer}.pkl'), 'rb') as pkl_in:
    dfs_mod_filtered = pickle.load(pkl_in)

thresh_delta = 5

dfs_mod_thresh_delta = {}
for this_mod in mod_names:
    this_merged_df = get_merged_df_with_delta(dfs_mod_filtered[this_mod])

    if f'{writer}-KD' in CONDITIONS:
        if this_mod == '17802':
            df_thresh_delta = this_merged_df[(this_merged_df[f'delta_{writer}-KD'] < -thresh_delta)]
        elif this_mod == 'a':
            df_thresh_delta = this_merged_df[(this_merged_df[f'delta_{writer}-KD'] >= thresh_delta)]
    elif f'{writer}-OE' in CONDITIONS:
        if this_mod == '17802':
            df_thresh_delta = this_merged_df[(this_merged_df[f'delta_{writer}-OE'] >= thresh_delta)]
        elif this_mod == 'a':
            df_thresh_delta = this_merged_df[(this_merged_df[f'delta_{writer}-OE'] < -thresh_delta)]

    dfs_mod_thresh_delta[this_mod] = df_thresh_delta

out_df = dfs_mod_thresh_delta['a'][dfs_mod_thresh_delta['a'][f'freq_{writer}-KD'] >= 80]
out_df.to_csv(os.path.join(base_dir, f'{writer}-KD_delta_sites_a.tsv'), sep='\t', index=False)

### output metaplot bed ###
out_bed = out_df.iloc[:, :6]
out_bed.to_csv(os.path.join(base_dir, 'm6a.bed'), sep='\t', index=False, header=False)
