import os
import numpy as np
import matplotlib as mpl
#######################################################################
cm = 1/2.54  # centimeters in inches
gr = 1.618
dpi = 1200
mpl.rcParams['figure.dpi'] = dpi
mpl.rcParams['savefig.dpi'] = dpi
mpl.rcParams['font.size'] = 6
mpl.rcParams['legend.fontsize'] = 4
mpl.rcParams['xtick.labelsize'] = 5
mpl.rcParams['ytick.labelsize'] = 5
mpl.rcParams['xtick.major.size'] = 1.5
mpl.rcParams['ytick.major.size'] = 1.5
mpl.rcParams['lines.linewidth'] = 0.5
# mpl.rcParams['font.family'] = 'Arial'
FMT = 'svg'
fig_kwargs = dict(format=FMT, bbox_inches='tight', dpi=dpi, transparent=True)
# FMT = 'png'
# fig_kwargs = dict(format=FMT, bbox_inches='tight', dpi=dpi)
######################################################################
import matplotlib.pyplot as plt
from scipy.stats import pearsonr, spearmanr, trim_mean


img_out = '/home/achan/img_out/m6A_psi_cross_talk'
os.makedirs(img_out, exist_ok=True)

dict_mod_display = {
    'm6A': 'm^6A',
    'psi': '\psi'
}

m6a_motifs = '_withoutT'

base_dir = '/home/achan/prj/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA_v1'
in_file = os.path.join(base_dir, f'single_read_frac_m6A{m6a_motifs}_psi_HEK293_WT_R002.npz')
temp = np.load(in_file)
vec_m6A = temp['vec_m6A']
vec_psi = temp['vec_psi']

binned_psi = []
binned_m6A = []
bin_edges = np.round(np.linspace(0, 1, 6), 1)
for bin_i in range(len(bin_edges)-1):
    bin_start = bin_edges[bin_i]
    bin_end = bin_edges[bin_i+1]
    if bin_end == 1.0:
        bin_end += 0.001
    mask_m6A = (vec_m6A >= bin_start) * (vec_m6A < bin_end)
    binned_psi.append(vec_psi[mask_m6A])
    mask_psi = (vec_psi >= bin_start) * (vec_psi < bin_end)
    binned_m6A.append(vec_m6A[mask_psi])

flierprops = dict(marker='o', markerfacecolor='none', markersize=2, markeredgecolor='gray',
                  alpha=0.5, rasterized=True)
xy_ticks = np.int32(bin_edges * 100)

# trim_q = 0.0
# trimmed_mean_m6A = [trim_mean(x, trim_q) for x in binned_m6A]
trimmed_mean_m6A = [np.median(x) for x in binned_m6A]
binned_centers = 0.5 * (bin_edges[1:] + bin_edges[:-1])
# spearman_rho, spearman_pval = spearmanr(trimmed_mean_m6A, binned_centers)
# pearson_rho, pearson_pval = pearsonr(trimmed_mean_m6A, binned_centers)

### boxplot ###
plt.figure(figsize=(5*cm, 5*cm))
plt.boxplot(binned_m6A, flierprops=flierprops)
plt.ylim([-0.01, 1.05])
plt.xticks(np.arange(len(bin_edges)) + 0.5, xy_ticks)
plt.yticks(bin_edges, xy_ticks)
plt.xlabel(f"occ(${dict_mod_display['psi']}$) per read")
plt.ylabel(f"occ(${dict_mod_display['m6A']}$) per read")
#plt.title(fr'$\rho$={pearson_rho:.3f}, p-value={pearson_pval:.3f}'+f'(\n{trim_q*100:.0f}% outlier trimmed)')
# plt.title(fr'$\rho$={pearson_rho:.3f}, p-value={pearson_pval:.3f}')
plt.savefig(os.path.join(img_out, f'boxplot{m6a_motifs}.{FMT}'), **fig_kwargs)

### outlier ###
# top_reads = 100
arr_top_reads = [100, 1000, 10000, 100000, 1000000]
arr_colors = [str(this_float) for this_float in np.linspace(0, 0.5, len(arr_top_reads))]
plt.figure(figsize=(5*cm, 5*cm))
for top_reads, top_color in zip(arr_top_reads, arr_colors):
    label = f'Top {top_reads:.0E} reads'
    top_psi = [np.median(np.sort(this_bin)[-top_reads:]) for this_bin in binned_psi]
    top_m6A = [np.median(np.sort(this_bin)[-top_reads:]) for this_bin in binned_m6A]
    plt.plot(np.arange(1, len(bin_edges)), top_m6A, c=top_color, ls='-', label=label)
    plt.plot(np.arange(1, len(bin_edges)), top_m6A, c=top_color, marker='o', markersize=2)
plt.legend(title='Median', loc='lower right')
plt.xticks(np.arange(len(bin_edges)) + 0.5, xy_ticks)
plt.ylim([-0.01, 1.05])
plt.yticks(bin_edges, xy_ticks)
plt.xlabel(f"occ(${dict_mod_display['psi']}$) per read")
plt.ylabel(f"occ(${dict_mod_display['m6A']}$) per read")
plt.tight_layout()
plt.savefig(os.path.join(img_out, f'outlier{m6a_motifs}.{FMT}'), **fig_kwargs)
