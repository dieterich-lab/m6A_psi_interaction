import os
import pandas as pd
pd.set_option('display.max_columns', None)
import numpy as np
import matplotlib as mpl
mpl.use('TkAgg')
import matplotlib.pyplot as plt
from tqdm import tqdm


def get_vec_change_neg_log_pval(in_dmr_file):
    change_pval = []
    for this_chunk in tqdm(pd.read_csv(in_dmr_file, sep='\t', dtype={'chrom': str}, iterator=True, chunksize=10000)):
        this_chunk_valid = this_chunk[this_chunk['p_value'] < 1]
        change_pval.append(this_chunk_valid[['effect_size', 'p_value']].values)
    vec_effect_size, vec_pval = np.vstack(change_pval).T
    vec_change = -vec_effect_size * 100
    vec_neg_log_pval = -np.log10(vec_pval)
    return vec_change, vec_neg_log_pval


def get_mask(in_vec_change, in_vec_neg_log_pval, in_thresh_neg_log_pval):
    mask_no_change = in_vec_neg_log_pval < in_thresh_neg_log_pval
    mask_pos_change = (in_vec_neg_log_pval >= in_thresh_neg_log_pval) * (in_vec_change >= 0)
    mask_neg_change = (in_vec_neg_log_pval >= in_thresh_neg_log_pval) * (in_vec_change < 0)

    num_pos_change = np.sum(mask_pos_change)
    num_neg_change = np.sum(mask_neg_change)

    return mask_no_change, mask_neg_change, mask_pos_change, num_neg_change, num_pos_change


home = os.environ['HOME']
data_dir = f'{home}/Data/TRR319_RMaP_BaseCalling_RNA004/Adrian'
img_out = f'{home}/img_out/RNA004_psi_KD_OE_analysis'
os.makedirs(img_out, exist_ok=True)

dict_display_mod = {
    'm6A': 'm^6A',
    'psi': '\Psi'
}

# exp = 'HEK293_psU-KD'
exp = 'HEK293_psU-OE'
cond0 = 'CTRL-KD'
cond1 = 'PUS7-OE'

thresh_neg_log_pval = 2.0

plt.figure(figsize=(10, 5))
for subplot_ind, mod in enumerate(['psi', 'm6A']):
    this_dmr_file = os.path.join(data_dir, exp, f'{cond0}_{cond1}.cov10.{mod}.diff_sites.dmr')
    this_vec_change, this_vec_neg_log_pval = get_vec_change_neg_log_pval(this_dmr_file)
    this_mask_no_change, this_mask_neg_change, this_mask_pos_change, this_num_neg_change, this_num_pos_change = get_mask(this_vec_change, this_vec_neg_log_pval, thresh_neg_log_pval)

    plt.subplot(1, 2, subplot_ind+1)
    plt.scatter(this_vec_change[this_mask_no_change], this_vec_neg_log_pval[this_mask_no_change], s=1, c='gray')
    plt.scatter(this_vec_change[this_mask_pos_change], this_vec_neg_log_pval[this_mask_pos_change], s=3, c='red')
    plt.scatter(this_vec_change[this_mask_neg_change], this_vec_neg_log_pval[this_mask_neg_change], s=3, c='blue')
    plt.axhline(y=thresh_neg_log_pval, c='gray', ls='--')
    plt.xlabel(f'% Mod. level change')
    plt.ylabel('$-log_{10}$ p-val')
    plt.title(f'${dict_display_mod[mod]}$')
    plt.xlim([-101, 101])
    plt.ylim([0, plt.gca().get_ylim()[1]])
    # plt.ylim([0, 20])
    plt.text(0.01, 1.01, f'{this_num_neg_change} down', c='blue', ha='left', transform=plt.gca().transAxes)
    plt.text(0.99, 1.01, f'{this_num_pos_change} up', c='red', ha='right', transform=plt.gca().transAxes)
plt.suptitle(f'{cond1} versus {cond0}')
plt.savefig(os.path.join(img_out, f'volcano_{exp}_{cond0}_{cond1}_diff_sites.png'), bbox_inches='tight')
