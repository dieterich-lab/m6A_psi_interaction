import os
import pandas as pd
pd.set_option('display.max_columns', None)
import numpy as np
import matplotlib as mpl
mpl.use('TkAgg')
import matplotlib.pyplot as plt
from tqdm import tqdm
from Bio import SeqIO
import logomaker


def get_filtered_df_by_motifs(in_df, in_motifs):
    df_out = []
    for _, this_row in in_df.iterrows():
        this_motif = ref[str(this_row['chrom'])][(this_row['chromStart'] - 2):(this_row['chromStart'] + 3)]
        if this_row['strand'] == '-':
            this_motif = this_motif.reverse_complement()
        if this_motif in in_motifs:
            df_out.append(this_row)
    return pd.DataFrame(df_out)

def get_vec_change_neg_log_pval(in_dmr_file, in_thresh_count, in_motifs):
    change_pval_score = []
    sites = []
    for this_chunk in tqdm(pd.read_csv(in_dmr_file, sep='\t', dtype={'chrom': str}, iterator=True, chunksize=10000)):
        if replicates:
            this_chunk_valid = this_chunk[
                (this_chunk['balanced_p_value'] < 1)
                * (this_chunk['sample_a_total'] >= in_thresh_count)
                * (this_chunk['sample_b_total'] >= in_thresh_count)
                ]
            this_chunk_valid = get_filtered_df_by_motifs(this_chunk_valid, in_motifs)
            if len(this_chunk_valid):
                change_pval_score.append(this_chunk_valid[['balanced_effect_size', 'balanced_p_value', 'score']].values)
        else:
            this_chunk_valid = this_chunk[
                (this_chunk['p_value'] < 1)
                * (this_chunk['sample_a_total'] >= in_thresh_count)
                * (this_chunk['sample_b_total'] >= in_thresh_count)
                ]
            this_chunk_valid = get_filtered_df_by_motifs(this_chunk_valid, in_motifs)
            if len(this_chunk_valid):
                change_pval_score.append(this_chunk_valid[['effect_size', 'p_value', 'score']].values)
        if len(this_chunk_valid):
            sites.append(this_chunk_valid[['chrom', 'chromStart', 'strand']])
    vec_effect_size, vec_pval, vec_score = np.vstack(change_pval_score).T
    vec_change = -vec_effect_size * 100
    vec_neg_log_pval = -np.log10(vec_pval+10**(-20))
    df_sites = pd.concat(sites)
    return vec_change, vec_neg_log_pval, vec_score, df_sites


def get_mask(in_vec_change, in_vec_neg_log_pval, in_thresh_change, in_thresh_neg_log_pval):
    # mask_no_change = in_vec_neg_log_pval < in_thresh_neg_log_pval
    mask_no_change = (in_vec_change < in_thresh_change) * (in_vec_change >= -in_thresh_change)
    mask_pos_change = (in_vec_neg_log_pval >= in_thresh_neg_log_pval) * (in_vec_change >= in_thresh_change)
    mask_neg_change = (in_vec_neg_log_pval >= in_thresh_neg_log_pval) * (in_vec_change < -in_thresh_change)
    num_pos_change = np.sum(mask_pos_change)
    num_neg_change = np.sum(mask_neg_change)
    return mask_no_change, mask_neg_change, mask_pos_change, num_neg_change, num_pos_change


REF_FILE = '/home/adrian/Data/genomes/homo_sapiens/GRCh38_102/GRCh38_102.fa'
def get_ref(ref_file=REF_FILE):
    out_ref = {}
    for record in SeqIO.parse(ref_file, "fasta"):
        out_ref[record.id] = record.seq
    return out_ref


def plot_logo(in_df_sites, in_mask, in_mod, in_change, span=3):
    sites_change = in_df_sites[in_mask]
    sites_change.to_csv(os.path.join(img_out, f'sites_{cond1}_{in_change}_{in_mod}.tsv'), sep='\t', index=False)
    all_refseq = []
    for _, row in sites_change.iterrows():
        refseq = ref[row['chrom']][(row['chromStart'] - span):(row['chromStart'] + span + 1)]
        if row['strand'] == '-':
            refseq = refseq.reverse_complement()
        all_refseq.append(refseq)
    df_logo = logomaker.alignment_to_matrix([str(seq) for seq in all_refseq])

    fig = plt.figure(figsize=(10, 5))
    ax = fig.subplots()
    logomaker.Logo(df_logo, ax=ax)
    ax.set_xticks(np.arange(0, 2 * span + 1, span), np.arange(0, 2 * span + 1, span) - span)
    ax.set_title(f'{in_mod}, {in_change}')
    fig.savefig(os.path.join(img_out, f'logo_{cond1}_{in_change}_{in_mod}_site_motifs.png'), bbox_inches='tight')
    plt.close(fig)


ref = get_ref()

home = os.environ['HOME']
data_dir = f'{home}/Data/TRR319_RMaP_BaseCalling_RNA004/Adrian'
img_out = f'{home}/img_out/RNA004_psi_KD_OE_analysis'
os.makedirs(img_out, exist_ok=True)

dict_display_mod = {
    'm6A': 'm^6A',
    'psi': '\Psi'
}

exp = 'HEK293_M3I'
# exp = 'HEK293_psU-OE'
# exp = 'HEK293_psU-KD'
# exp = 'HEK293_M3I'
cond0 = 'DMSO_rep1'
cond1 = 'M3I_24h_rep1'
# cond0 = 'CTRL'
# cond1 = 'PUS7_OE'
# cond1 = 'TRUB1-KD'
# cond0 = 'CTRL_OE'
# cond1 = 'CTRL_KD'
replicates = False

thresh_neg_log_pval = 2.0
thresh_score = 2.0
thresh_count = 50
thresh_change = 25.0

mod_motifs = {
    'm6A': [
            'GGACT', 'GGACA', 'GAACT', 'AGACT', 'GGACC', 'TGACT',
            'AAACT', 'GAACA', 'AGACA', 'AGACC', 'GAACC', 'TGACA',
            'TAACT', 'AAACA', 'TGACC', 'TAACA', 'AAACC', 'TAACC'
    ],
    'psi': ['GTTCA', 'GTTCC', 'GTTCG', 'GTTCT'] + ['TGTAG'] +
           [f'{pos1}{pos2}T{pos4}{pos5}'
            for pos1 in ['A', 'C', 'G', 'T']
            for pos2 in ['A', 'G']
            for pos4 in ['A', 'G']
            for pos5 in ['A', 'C', 'G', 'T']
            ]
}

plt.figure(figsize=(10, 5))
for subplot_ind, mod in enumerate(['psi', 'm6A']):
    if replicates:
        this_dmr_file = os.path.join(data_dir, exp, 'dmr', f'{cond0}_{cond1}_replicates.cov10.{mod}.diff_sites.dmr')
    else:
        this_dmr_file = os.path.join(data_dir, exp, 'dmr', f'{cond0}_{cond1}.cov10.{mod}.diff_sites.dmr')
    this_vec_change, this_vec_neg_log_pval, this_vec_score, this_df_sites = get_vec_change_neg_log_pval(this_dmr_file, thresh_count, mod_motifs[mod])
    this_mask_no_change, this_mask_neg_change, this_mask_pos_change, this_num_neg_change, this_num_pos_change = get_mask(this_vec_change, this_vec_neg_log_pval, thresh_change, thresh_neg_log_pval)
    for this_change, this_mask in zip(['up', 'down'], [this_mask_pos_change, this_mask_neg_change]):
        if this_mask.any():
            plot_logo(this_df_sites, this_mask, mod, this_change)

    plt.subplot(1, 2, subplot_ind+1)
    plt.scatter(this_vec_change[this_mask_no_change], this_vec_neg_log_pval[this_mask_no_change], s=1, c='gray')
    plt.scatter(this_vec_change[this_mask_pos_change], this_vec_neg_log_pval[this_mask_pos_change], s=3, c='red')
    plt.scatter(this_vec_change[this_mask_neg_change], this_vec_neg_log_pval[this_mask_neg_change], s=3, c='blue')
    plt.axhline(y=thresh_neg_log_pval, c='gray', ls='--')
    plt.axvline(x=-thresh_change, c='gray', ls='--')
    plt.axvline(x=thresh_change, c='gray', ls='--')
    plt.xlabel(f'% Mod. level change')
    plt.ylabel('$-log_{10}$ p-val')
    plt.title(f'${dict_display_mod[mod]}$')
    plt.xlim([-101, 101])
    plt.ylim([0, plt.gca().get_ylim()[1]])
    # plt.ylim([0, 20])
    plt.text(0.01, 1.01, f'{this_num_neg_change} down', c='blue', ha='left', transform=plt.gca().transAxes)
    plt.text(0.99, 1.01, f'{this_num_pos_change} up', c='red', ha='right', transform=plt.gca().transAxes)
if replicates:
    plt.suptitle(f'{cond1} versus {cond0}, count$\geq${thresh_count}, with replicates')
    plt.savefig(os.path.join(img_out, f'volcano_{exp}_{cond0}_{cond1}_replicates_diff_sites.png'), bbox_inches='tight')
else:
    plt.suptitle(f'{cond1} versus {cond0}, count$\geq${thresh_count}')
    plt.savefig(os.path.join(img_out, f'volcano_{exp}_{cond0}_{cond1}_count{thresh_count}_diff_sites.png'), bbox_inches='tight')
