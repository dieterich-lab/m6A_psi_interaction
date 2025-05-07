import os
import pandas as pd
import pysam
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
mpl.rcParams['font.family'] = 'Arial'
# FMT = 'svg'
# fig_kwargs = dict(format=FMT, bbox_inches='tight', dpi=dpi, transparent=True)
FMT = 'png'
fig_kwargs = dict(format=FMT, bbox_inches='tight', dpi=dpi)
######################################################################
mpl.use('TkAgg')
import matplotlib.pyplot as plt
from tqdm import tqdm

img_out = '/home/adrian/img_out/RNA004_psi_KD_OE_analysis'
os.makedirs(img_out, exist_ok=True)

dict_mod_display = {
    'm6A': 'm^6A',
    'psi': '\psi'
}

thresh_valid_reads = 1000
num_top_locs = 5

def get_mean_logit(in_probs):
    if len(in_probs) == 0:
        return np.nan
    top_probs = np.sort(in_probs)[-num_top_locs:]
    rescaled_probs = np.clip(np.array(top_probs) / 255.0, a_max=0.999, a_min=0.001)
    logits = np.log2(rescaled_probs / (1-rescaled_probs))
    return np.mean(logits)


def get_mean_logit_mod_level(in_read):
    mod_mean_logit = {}
    for this_mod, this_tag in mod_tags.items():
        this_mod_probs = [this_tup[1] for this_tup in in_read.modified_bases.get(this_tag, [])]
        mod_mean_logit[this_mod] = get_mean_logit(this_mod_probs)
    return mod_mean_logit

thresh_min_locs = 10
def get_mod_mean_occupancy(in_read, min_locs=thresh_min_locs):
    mod_mean_occupancy = {}
    for this_mod, this_tag in mod_tags.items():
        this_mod_probs = np.array([this_tup[1] for this_tup in in_read.modified_bases.get(this_tag, [])]) / 255.0
        if len(this_mod_probs) >= min_locs:
            mod_mean_occupancy[this_mod] = np.mean(this_mod_probs >= 0.5)
        else:
            mod_mean_occupancy[this_mod] = np.nan
    return mod_mean_occupancy

########################################################################################################################
### R002 ###############################################################################################################
########################################################################################################################
# mod_tags = {
#     'm6A': ('N', 0, 21891),
#     'psi': ('N', 0, 17802)
# }
#
# base_dir = '/home/adrian/Data/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA_v1'
#
# ds = 'HEK293_WT_R002'
# bam_file = os.path.join(base_dir, 'HEK293/WT_P2/chrALL.mAFiA.reads.bam')

# ds = 'M3KO'
# bam_file = os.path.join(base_dir, 'HEK293T_Mettl3_KO/merged/chrALL.mAFiA.reads.bam')

# ds = 'M3KD'
# bam_file = os.path.join(base_dir, 'NanoSPA/HEK_siMETTL3_input_merged/chrALL.mAFiA.reads.bam')

# ds = 'TRUB1KD'
# bam_file = os.path.join(base_dir, 'NanoSPA/HEK_siTRUB1_input_merged/chrALL.mAFiA.reads.bam')

# ds = 'TRUB1OE'
# bam_file = os.path.join(base_dir, 'HEK293_TRUB1_OE/merged/chrALL.mAFiA.reads.bam')

########################################################################################################################
### R004 ###############################################################################################################
########################################################################################################################
mod_tags = {
    'm6A': ('A', 0, 'a'),
    'psi': ('T', 0, 17802)
}

base_dir = '/home/adrian/Data/TRR319_RMaP_BaseCalling_RNA004/Isabel/20250224_HEK293_psU_kds_RTA/Dorado_082'

ds = 'HEK293_ctrl_R004'
bam_file = os.path.join(base_dir, 'HEK293_ctrl_RTA/calls_2025-02-26_T06-44-51.bam')

# ds = 'HEK293_TRUB1_kd'
# bam_file = os.path.join(base_dir, 'HEK293_TRUB1_kd_RTA/calls_2025-02-26_T06-43-59.bam')

########################################################################################################################

single_read_mean_occupancy = []
with pysam.AlignmentFile(bam_file, 'rb', check_sq=False) as bam:
    for this_read in tqdm(bam.fetch(until_eof=True)):
        if this_read.modified_bases is not None:
            # single_read_mean_logit.append(get_mean_logit_mod_level(this_read))
            single_read_mean_occupancy.append(get_mod_mean_occupancy(this_read))

vec_m6A, vec_psi = np.vstack([
    (this_read_mean_occupancy['m6A'],  this_read_mean_occupancy['psi'])
    for this_read_mean_occupancy in single_read_mean_occupancy
    if ~np.isnan(this_read_mean_occupancy['m6A']) and ~np.isnan(this_read_mean_occupancy['psi'])
]).T
num_valid_reads = len(vec_m6A)

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

outfile_name = os.path.join(img_out, f'boxplot_mean_occupancy_per_read_{ds}.{FMT}')
# outfile_name = os.path.join(img_out, f'figureS4.{FMT}')

plt.figure(figsize=(8*cm, 7*cm))
plt.subplot(2, 2, 1)
plt.boxplot(binned_psi, flierprops=flierprops)
# plt.violinplot(binned_psi, quantiles=[[0.5]]*len(binned_psi))
# plt.plot(np.arange(1, len(bin_edges)), top_psi, c='r', ls='-')
# plt.plot(np.arange(1, len(bin_edges)), top_psi, 'r+', markersize=2, label=f'N$\geq${thresh_top_reads}')
# bin_sizes = [len(this_bin) for this_bin in binned_psi]
# for bin_ind, this_bin_size in enumerate(bin_sizes):
#     plt.text(bin_ind+0.5, 1.05, this_bin_size)
# plt.legend(loc='upper right')
plt.ylim([-0.01, 1.05])
plt.xticks(np.arange(len(bin_edges)) + 0.5, xy_ticks)
plt.yticks(bin_edges, xy_ticks)
plt.xlabel(f"frac(${dict_mod_display['m6A']}$) per read")
plt.ylabel(f"frac(${dict_mod_display['psi']}$) per read")
plt.subplot(2, 2, 2)
plt.boxplot(binned_m6A, flierprops=flierprops)
# plt.violinplot(binned_m6A, quantiles=[[0.5]]*len(binned_m6A))
# plt.plot(np.arange(1, len(bin_edges)), top_m6A, c='r', ls='-')
# plt.plot(np.arange(1, len(bin_edges)), top_m6A, 'r+', markersize=2, label=f'N$\geq${thresh_top_reads}')
# plt.legend(loc='upper right')
plt.ylim([-0.01, 1.05])
plt.xticks(np.arange(len(bin_edges)) + 0.5, xy_ticks)
plt.yticks(bin_edges, xy_ticks)
plt.xlabel(f"frac(${dict_mod_display['psi']}$) per read")
plt.ylabel(f"frac(${dict_mod_display['m6A']}$) per read")
# plt.savefig(os.path.join(img_out, f'boxplot_mean_occupancy_per_read_top{num_top_locs}_{ds}.{FMT}'), **fig_kwargs)
# plt.suptitle(f'{ds}\nMin. {thresh_min_locs} locs per read')

# thresh_top_reads = 0.75
# top_psi = [np.mean(this_bin[this_bin >= thresh_top_reads]) for this_bin in binned_psi]
# top_m6A = [np.mean(this_bin[this_bin >= thresh_top_reads]) for this_bin in binned_m6A]
top_reads = 100
label = f'Median, top {top_reads}'
top_psi = [np.median(np.sort(this_bin)[-top_reads:]) for this_bin in binned_psi]
top_m6A = [np.median(np.sort(this_bin)[-top_reads:]) for this_bin in binned_m6A]

# ylim = [0.79, 0.85]
# yticks = np.linspace(*ylim, 3)
top_color = 'k'
# plt.figure(figsize=(8*cm, 4*cm))
plt.subplot(2, 2, 3)
plt.plot(np.arange(1, len(bin_edges)), top_psi, c=top_color, ls='-', label=label)
plt.plot(np.arange(1, len(bin_edges)), top_psi, f'{top_color}o', markersize=2)
plt.legend(loc='lower left')
plt.xticks(np.arange(len(bin_edges)) + 0.5, xy_ticks)
plt.ylim([-0.01, 1.05])
plt.yticks(bin_edges, xy_ticks)
plt.xlabel(f"frac(${dict_mod_display['m6A']}$) per read")
plt.ylabel(f"frac(${dict_mod_display['psi']}$) per read")
# plt.ylim(ylim)
plt.subplot(2, 2, 4)
plt.plot(np.arange(1, len(bin_edges)), top_m6A, c=top_color, ls='-', label=label)
plt.plot(np.arange(1, len(bin_edges)), top_m6A, f'{top_color}o', markersize=2)
plt.legend(loc='lower left')
plt.xticks(np.arange(len(bin_edges)) + 0.5, xy_ticks)
plt.ylim([-0.01, 1.05])
plt.yticks(bin_edges, xy_ticks)
plt.xlabel(f"frac(${dict_mod_display['psi']}$) per read")
plt.ylabel(f"frac(${dict_mod_display['m6A']}$) per read")
# plt.ylim(ylim)
# plt.suptitle(f'N$\geq${int(thresh_top_reads*100)}%')
plt.tight_layout()
# plt.savefig(os.path.join(img_out, f'boxplot_mean_occupancy_per_read_{ds}_above{thresh_top_reads}.{FMT}'), **fig_kwargs)
plt.savefig(outfile_name, **fig_kwargs)
