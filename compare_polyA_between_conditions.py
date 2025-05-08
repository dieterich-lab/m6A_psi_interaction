import os
import pandas as pd
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)
import numpy as np
import pysam
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from collections import Counter
from tqdm import tqdm

########################################################################################################################
### define datasets ####################################################################################################
########################################################################################################################
conditions = ['HEK293_CTRL', 'HEK293_TRUB1_OE']
cond_colors = {this_cond: this_color for this_cond, this_color in zip(conditions, ['b', 'r'])}

polyA_dir = '/home/adrian/Data/TRR319_RMaP_BaseCalling_RNA004/Adrian/polyA/'
polyA_files = {this_cond: os.path.join(polyA_dir, f'read_polyA_{this_cond}.tsv') for this_cond in conditions}

home = os.environ['HOME']
img_out = f'{home}/img_out/RNA004_psi_KD_OE_analysis'
os.makedirs(img_out, exist_ok=True)

########################################################################################################################
### overall polyA distribution #########################################################################################
########################################################################################################################
df_polyA = {this_cond: pd.read_csv(polyA_files[this_cond], sep='\t', names=['read_id', 'polyA_len'])
            for this_cond in conditions}
dict_polyA = {this_cond: {k: v for k, v in df_polyA[this_cond][['read_id', 'polyA_len']].values}
              for this_cond in conditions}

xlim = [0, 300]
xticks = np.linspace(*xlim, 4)
num_bins = 100
bin_edges = np.linspace(*xlim, num_bins+1)
bin_centers = 0.5 * (bin_edges[1:] + bin_edges[:-1])

plt.figure(figsize=(4, 4))
for cond_ind, this_cond in enumerate(conditions):
    this_hist, _ = np.histogram(list(dict_polyA[this_cond].values()), bins=bin_edges)
    num_reads = this_hist.sum()
    norm_hist = this_hist / np.sum(this_hist)
    plt.plot(bin_centers, norm_hist, c=cond_colors[this_cond], label=f'{this_cond} ({num_reads})')

# ymax = (np.max(norm_hist) // 0.05 + 1) * 0.05
# ylim = [0, ymax]
# yticks = np.linspace(*ylim, 3)
# plt.ylim(ylim)
# plt.yticks(yticks)

plt.legend(fontsize=10)
plt.xlim(xlim)
plt.xticks(xticks)
plt.xlabel('polyA length (bps)', fontsize=12)
plt.ylabel('Density', fontsize=12)
plt.savefig(
    os.path.join(img_out, f'hist_{conditions[0]}_vs_{conditions[1]}_polyA_length.png'),
    bbox_inches='tight'
)
