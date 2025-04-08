from Bio import SeqIO
from tqdm import tqdm
import logomaker as lm
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


def get_ref(ref_file=REF_FILE):
    out_ref = {}
    for record in SeqIO.parse(ref_file, "fasta"):
        out_ref[record.id] = record.seq
    return out_ref


def get_central_motif(in_df, in_ref, span=2):
    print(f'\nAssiging motif to {len(in_df)} sites:')
    all_motifs = []
    for _, this_row in tqdm(in_df.iterrows()):
        this_chrom, this_chromStart, this_chromEnd, this_strand = this_row[
            ['chrom', 'chromStart', 'chromEnd', 'strand']
        ]
        this_motif = in_ref[this_chrom][(this_chromStart-span):(this_chromStart+span+1)]
        if this_strand == '-':
            this_motif = this_motif.reverse_complement()
        all_motifs.append(str(this_motif))
    in_df['ref_motif'] = all_motifs
    return in_df


def get_df_from_bed(in_base_dir, in_ds, in_sample, thresh_cov=100):
    bed_fields = [
        'chrom',
        'chromStart',
        'chromEnd',
        'name',
        'score',
        'strand',
        'frequency'
    ]
    out_df = pd.read_csv(os.path.join(in_base_dir, in_ds, f'{in_sample}.modkit042.cov10.bedmethyl'), sep='\t',
                         usecols=[0, 1, 2, 3, 4, 5, 10], names=bed_fields, dtype={'chrom': str})
    out_df = out_df[out_df['score'] >= thresh_cov]
    return out_df


def get_df_mod_filtered(in_df, in_mod, in_ref, in_motifs):
    this_df_mod = in_df[in_df['name'] == in_mod]
    this_df_mod = get_central_motif(this_df_mod, in_ref)
    return this_df_mod[this_df_mod['ref_motif'].isin(in_motifs)]


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
    for cond in conditions:
        df = in_dfs[cond].copy()
        this_mod_dfs.append(
            df.rename(columns={'score': f'score_{cond}', 'frequency': f'freq_{cond}'})
        )
    out_df_merged = reduce(lambda left, right: pd.merge(left, right, on=merged_fields, how='inner'), this_mod_dfs)
    out_df_merged = out_df_merged[
        (out_df_merged['freq_KD'] > 0) *
        (out_df_merged['freq_WT'] > 0) *
        (out_df_merged['freq_OE'] > 0)
        ]
    out_df_merged['delta_KD'] = out_df_merged['freq_KD'] - out_df_merged['freq_WT']
    out_df_merged['delta_OE'] = out_df_merged['freq_OE'] - out_df_merged['freq_WT']
    return out_df_merged


img_out = IMG_OUT
os.makedirs(img_out, exist_ok=True)

### RNA004 ###
base_dir = BASE_DIR
writer = WRITER

sel_motifs = {
    'a': {
        'DRACH': [
            'GGACT', 'GGACA', 'GAACT', 'AGACT', 'GGACC', 'TGACT',
            'AAACT', 'GAACA', 'AGACA', 'AGACC', 'GAACC', 'TGACA',
            'TAACT', 'AAACA', 'TGACC', 'TAACA', 'AAACC', 'TAACC'
        ]
    },
    '17802': {
        'TRUB1': ['GTTCA', 'GTTCC', 'GTTCG', 'GTTCT'],
        'PUS7': ['TGTAG']
    }
}

mod_names = ['17802', 'a']
dict_mod_display = {
    'a': 'm^6A',
    '17802': '\psi'
}

cond_colors = {
    'WT': 'gray',
    'OE': 'red',
    'KD': 'blue'
}

dfs = {}
dfs['WT'] = get_df_from_bed(base_dir, 'HEK293_psU-KD', 'CTRL')
dfs['KD'] = get_df_from_bed(base_dir, 'HEK293_psU-KD', f'{writer}-KD')
dfs['OE'] = get_df_from_bed(base_dir, 'HEK293_psU-OE', f'{writer}-OE')
# conditions = list(dfs.keys())
conditions = ['KD', 'WT', 'OE']

ref = get_ref()

### histogram ###
dfs_mod_filtered = {
    this_mod: {} for this_mod in mod_names
}

mod_motif_name = {
    '17802': writer,
    'a': 'DRACH'
}

### historgram of S distribution ###
num_bins = 5
for sel_mod in mod_names:
    plt.figure(figsize=(5*cm, 5*cm))
    for cond_ind, this_cond in enumerate(conditions):
        if dfs_mod_filtered[sel_mod].get(this_cond) is None:
            this_df = dfs[this_cond]
            this_df_mod_filtered = get_df_mod_filtered(this_df, sel_mod, ref, sel_motifs[sel_mod][mod_motif_name[sel_mod]])
            dfs_mod_filtered[sel_mod][this_cond] = this_df_mod_filtered
        else:
            this_df_mod_filtered = dfs_mod_filtered[sel_mod][this_cond]
        this_hist, bin_edges = np.histogram(this_df_mod_filtered['frequency'], range=[0, 100], bins=num_bins)
        this_hist_norm = this_hist / np.sum(this_hist)
        bin_centers = 0.5 * (bin_edges[1:] + bin_edges[:-1])
        plt.semilogy(bin_centers, this_hist_norm, label=f'{writer}-{this_cond}', c=cond_colors[this_cond])
        # plt.bar(bin_centers+(cond_ind-1)*5, this_hist_norm, label=f'{writer}-{this_cond}',
        #         fc=cond_colors[this_cond], width=5, log=True)
    plt.legend()
    plt.xticks(np.linspace(0, 100, num_bins+1))
    plt.xlabel('S')
    plt.ylabel('Probility (log)')
    plt.title(f'${dict_mod_display[sel_mod]}$, {mod_motif_name[sel_mod]} motifs')
    plt.savefig(os.path.join(img_out, f'histogram_mod_{sel_mod}_{WRITER}.{FMT}'), **fig_kwargs)

# with open(os.path.join(base_dir, 'dfs_mod_filtered.pkl'), 'wb') as pkl_out:
#     pickle.dump(dfs_mod_filtered, pkl_out)

### merged and plot delta ###
with open(os.path.join(base_dir, 'dfs_mod_filtered.pkl'), 'rb') as pkl_in:
    dfs_mod_filtered = pickle.load(pkl_in)

dfs_mod_merged = {}
for this_mod in mod_names:
    this_merged_df = get_merged_df_with_delta(dfs_mod_filtered[this_mod])
    dfs_mod_merged[this_mod] = this_merged_df

    this_merged_df = this_merged_df[
        (np.abs(this_merged_df['delta_KD']) >= 5) +
        (np.abs(this_merged_df['delta_OE']) >= 5)
        ]

    plt.figure(figsize=(5*cm, 5*cm))
    # plt.scatter(this_merged_df['freq_WT'], this_merged_df['delta_KD'], s=0.5, c=cond_colors['KD'])
    # plt.scatter(this_merged_df['freq_WT'], this_merged_df['delta_OE'], s=0.5, c=cond_colors['OE'])
    for _, this_row in this_merged_df.iterrows():
        plt.plot(this_row['freq_WT'], this_row['delta_OE'], '.', markersize=1, c=cond_colors['OE'])
        plt.plot(this_row['freq_WT'], this_row['delta_KD'], '.', markersize=1, c=cond_colors['KD'])
        if this_mod == '17802':
            plt.plot([this_row['freq_WT']] * 2, [this_row['delta_KD'], this_row['delta_OE']],
                     c=cond_colors['WT'], linewidth=0.5, alpha=0.5)
    plt.ylim([-100, 100])
    plt.xlabel('$S_{WT}$')
    plt.ylabel('$\Delta$ from $S_{WT}$')
    plt.title(f'${dict_mod_display[this_mod]}$, {mod_motif_name[this_mod]} motifs')

    point = Line2D([0], [0], label='OE', marker='o', markersize=1,
                   markerfacecolor=cond_colors['OE'], linestyle='')


    # handles, labels = plt.gca().get_legend_handles_labels()
    handles = [Line2D([0], [0], label=f'{writer}-{this_cond}',
                      marker='.', markersize=5, markerfacecolor=cond_colors[this_cond], markeredgecolor='None',
                      linestyle='')
               for this_cond in ['KD', 'OE']
               ]
    plt.legend(handles=handles, loc='lower left')

    plt.savefig(os.path.join(img_out, f'scatter_delta_vs_S_WT_mod_{this_mod}.{FMT}'), **fig_kwargs)


