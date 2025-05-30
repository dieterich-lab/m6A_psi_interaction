import os
import pandas as pd
import pickle
from functools import reduce
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
import pybedtools
import numpy as np
from tqdm import tqdm
from scipy import stats


def get_merged_df_with_delta(in_dfs, compare_conds):
    merged_fields = [
        'chrom',
        'chromStart',
        'chromEnd',
        'name',
        'strand',
        'ref_motif'
    ]

    this_mod_dfs = []
    for cond in compare_conds:
        df = in_dfs[cond].copy()
        this_mod_dfs.append(
            df.rename(columns={'score': f'score_{cond}', 'frequency': f'freq_{cond}'})
        )
    out_df_merged = reduce(lambda left, right: pd.merge(left, right, on=merged_fields, how='inner'), this_mod_dfs)

    mask = ((out_df_merged.loc[:, out_df_merged.columns.str.contains('freq_')] > 0).any(axis=1))
    out_df_merged = out_df_merged[mask]

    if f'{writer}-KD' in compare_conds:
        out_df_merged[f'delta_{writer}-KD'] = out_df_merged[f'freq_{writer}-KD'] - out_df_merged['freq_CTRL']
    if f'{writer}-OE' in compare_conds:
        out_df_merged[f'delta_{writer}-OE'] = out_df_merged[f'freq_{writer}-OE'] - out_df_merged['freq_CTRL']

    return out_df_merged


def get_neighboring_sites_on_exon(in_row, dict_annots, second_mod_bedtool):
    if in_row['chrom'] not in dict_annots.keys():
        return []
    in_row_bedtool = pybedtools.BedTool.from_dataframe(in_row.to_frame().T)
    this_chr_annot = dict_annots[in_row['chrom']]
    in_row_exon = this_chr_annot.intersect(in_row_bedtool, u=True)
    if len(in_row_exon) > 1:
        in_row_exon = in_row_exon.sort().merge()
    neighboring_sites = second_mod_bedtool.intersect(in_row_exon)
    return neighboring_sites.to_dataframe()


def get_vec_dist_delta_from_site_df(in_df, sec_bedtool, annots):
    dist_delta = []
    print(f'Processing {len(in_df)} sites...')
    for _, this_row in tqdm(in_df.iterrows()):
        exon_sec_sites = get_neighboring_sites_on_exon(this_row, annots, sec_bedtool)
        if len(exon_sec_sites):
            this_vec_dist = (this_row['chromStart'] - exon_sec_sites['start']).abs().to_numpy()
            this_vec_delta = exon_sec_sites['score'].to_numpy()
            dist_delta.extend(list(zip(this_vec_dist, this_vec_delta)))
    return np.vstack(dist_delta).T


def get_binned_dist_delta(in_vec_dist, in_vec_delta, in_bin_range, in_bin_width):
    bin_edges = np.arange(in_bin_range[0], in_bin_range[1]+in_bin_width, in_bin_width)
    bin_centers = 0.5 * (bin_edges[1:] + bin_edges[:-1])
    num_bins = len(bin_edges) - 1
    binned_delta = []
    for bin_i in range(num_bins):
        bin_start = bin_edges[bin_i]
        bin_end = bin_edges[bin_i+1]
        mask_dist = (in_vec_dist >= bin_start) * (in_vec_dist < bin_end)
        if mask_dist.any():
            # binned_delta.append(np.median(in_vec_delta[mask_dist]))
            # binned_delta.append(np.mean(in_vec_delta[mask_dist]))
            # binned_delta.append(np.max(in_vec_delta[mask_dist]))
            binned_delta.append(in_vec_delta[mask_dist])
        else:
            binned_delta.append(np.nan)
    return bin_centers, binned_delta


base_dir = '/home/adrian/Data/TRR319_RMaP_BaseCalling_RNA004/Adrian'
writer = 'PUS7'
cond = 'KD'
# cond = 'OE'

img_out = '/home/adrian/img_out/RNA004_psi_KD_OE_analysis'

with open(os.path.join(base_dir, f'dfs_mod_filtered_{writer}.pkl'), 'rb') as pkl_in:
    dfs_mod_cond = pickle.load(pkl_in)

annot_dir = '/home/adrian/Data/genomes/homo_sapiens/GRCh38_102/exons'
chr_annots = {
    str(this_chr): pybedtools.BedTool(os.path.join(annot_dir, f'chr{this_chr}.exons.GRCh38.102.gtf'))
    for this_chr in list(range(1, 23)) + ['X', 'Y', 'MT']
}

merged_df_psi = get_merged_df_with_delta(dfs_mod_cond['17802'], ['CTRL', f'{writer}-{cond}'])
merged_df_m6a = get_merged_df_with_delta(dfs_mod_cond['a'], ['CTRL', f'{writer}-{cond}'])

plt.figure(figsize=(10*cm, 5*cm))
plt.subplot(1, 2, 1)
plt.hist(merged_df_psi[f'delta_{writer}-{cond}'], range=[-100, 100], bins=50, log=True)
plt.xlabel('$\Delta$S($\psi$)')
plt.ylabel('Site count')
plt.axvline(x=0, c='r')
plt.subplot(1, 2, 2)
plt.hist(merged_df_m6a[f'delta_{writer}-{cond}'], range=[-100, 100], bins=50, log=True)
plt.xlabel('$\Delta$S(m6A)')
plt.ylabel('Site count')
plt.axvline(x=0, c='r')
plt.suptitle(f'{writer}-{cond} vs CTRL')
plt.tight_layout()
plt.savefig(os.path.join(img_out, f'hist_mods_{writer}-{cond}.{FMT}'), **fig_kwargs)

merged_df_m6a.iloc[:, 4] = merged_df_m6a[f'delta_{writer}-{cond}']
m6a_bedtool = pybedtools.BedTool.from_dataframe(merged_df_m6a.iloc[:, :6])

thresh_delta = 5

# df_psi_same = merged_df_psi[
#     (merged_df_psi[f'delta_{writer}-{cond}'] >= -thresh_delta)
#     * (merged_df_psi[f'delta_{writer}-{cond}'] < thresh_delta)
#     ]

bin_range = [0, 2000]
bin_width = 200

if cond == 'KD':
    df_psi_down = merged_df_psi[merged_df_psi[f'delta_{writer}-{cond}'] < -thresh_delta]
    num_sites = len(df_psi_down)
    psi_down_vec_dist, psi_down_vec_delta = get_vec_dist_delta_from_site_df(df_psi_down, m6a_bedtool, chr_annots)
    psi_down_binned_dist, psi_down_binned_delta = get_binned_dist_delta(
        psi_down_vec_dist, psi_down_vec_delta, bin_range, bin_width
    )
    psi_down_binned_delta_lower = np.array([np.quantile(this_bin, 0.01) for this_bin in psi_down_binned_delta])
    psi_down_binned_delta_upper = np.array([np.quantile(this_bin, 0.99) for this_bin in psi_down_binned_delta])
    psi_down_binned_delta_mean = 0.5 * (psi_down_binned_delta_upper + psi_down_binned_delta_lower)
elif cond == 'OE':
    df_psi_up = merged_df_psi[merged_df_psi[f'delta_{writer}-{cond}'] >= thresh_delta]
    num_sites = len(df_psi_up)
    psi_up_vec_dist, psi_up_vec_delta = get_vec_dist_delta_from_site_df(df_psi_up, m6a_bedtool, chr_annots)
    psi_up_binned_dist, psi_up_binned_delta = get_binned_dist_delta(
        psi_up_vec_dist, psi_up_vec_delta, bin_range, bin_width
    )
    psi_up_binned_delta_lower = np.array([np.quantile(this_bin, 0.01) for this_bin in psi_up_binned_delta])
    psi_up_binned_delta_upper = np.array([np.quantile(this_bin, 0.99) for this_bin in psi_up_binned_delta])
    psi_up_binned_delta_mean = 0.5 * (psi_up_binned_delta_upper + psi_up_binned_delta_lower)
# psi_same_vec_dist, psi_same_vec_delta = get_vec_dist_delta_from_site_df(df_psi_same, m6a_bedtool, chr_annots)

# psi_same_binned_dist, psi_same_binned_delta = get_binned_dist_delta(
#     psi_same_vec_dist, psi_same_vec_delta, bin_range, bin_width
# )
# psi_down_binned_delta_mean = [stats.trim_mean(this_bin, 0.25) for this_bin in psi_down_binned_delta]
# psi_up_binned_delta_mean = [stats.trim_mean(this_bin, 0.25) for this_bin in psi_up_binned_delta]
# psi_up_binned_delta_mean = psi_up_binned_delta_lower - psi_up_binned_delta_upper

xticks = np.int64(np.linspace(*bin_range, 5))

plt.figure(figsize=(5*cm, 5*cm))
# plt.scatter(psi_up_vec_dist, psi_up_vec_delta, s=1, c='gray', label=f'$\Delta$S($\psi$)$\geq${thresh_delta}')
if cond == 'OE':
    plt.scatter(psi_up_vec_dist, psi_up_vec_delta, s=1, c='gray')
    plt.axhline(y=0, c='g', ls='--')
    plt.plot(psi_up_binned_dist, psi_up_binned_delta_mean, c='r', label='Trimmed mean')
    plt.plot(psi_up_binned_dist, psi_up_binned_delta_mean, 'r.')
elif cond == 'KD':
    plt.scatter(psi_down_vec_dist, psi_down_vec_delta, s=1, c='gray')
    plt.axhline(y=0, c='g', ls='--')
    plt.plot(psi_down_binned_dist, psi_down_binned_delta_mean, c='r', label='Trimmed mean')
    plt.plot(psi_down_binned_dist, psi_down_binned_delta_mean, 'r.')
# plt.scatter(psi_down_vec_dist, psi_down_vec_delta, s=1, c='b', label=f'$\Delta$S($\psi$) < -{thresh_delta}')
# plt.plot(psi_down_binned_dist, psi_down_binned_delta, 'b', label=f'$\Delta$S($\psi$) < -{thresh_delta}')
# plt.plot(psi_same_binned_dist, psi_same_binned_delta, 'gray', label=f'-{thresh_delta} $\leq$ $\Delta$S($\psi$) < {thresh_delta}')
# plt.plot(psi_up_binned_dist, psi_up_binned_delta, 'r', label=f'$\Delta$S($\psi$)$\geq${thresh_delta}')
# plt.boxplot(psi_down_binned_delta, positions=psi_down_binned_dist, showfliers=False, label=f'$\Delta$S($\psi$) < -{thresh_delta}')
# plt.boxplot(psi_up_binned_delta, positions=psi_up_binned_dist, widths=1, label=f'$\Delta$S($\psi$)$\geq${thresh_delta}')
plt.xlim(bin_range)
plt.ylim([-25, 25])
plt.xticks(xticks)
plt.xlabel('Distance from $\psi$ site on same exon (nts)')
plt.ylabel('$\Delta$S(m6A)')
# plt.plot(psi_up_binned_dist, psi_up_binned_delta_lower, c='b')
if cond == 'OE':
    plt.title(f'{writer}-{cond} vs CTRL\n{num_sites} sites with $\Delta$S($\psi$)$\geq${thresh_delta}')
elif cond == 'KD':
    plt.title(f'{writer}-{cond} vs CTRL\n{num_sites} sites with $\Delta$S($\psi$)<-{thresh_delta}')
plt.legend()
plt.savefig(os.path.join(img_out, f'dist_corr_m6a_from_psi_{writer}-{cond}_thresh{thresh_delta}.{FMT}'), **fig_kwargs)
# plt.savefig(os.path.join(img_out, f'boxplot_m6a_from_psi_{writer}-{cond}_thresh{thresh_delta}.{FMT}'), **fig_kwargs)
