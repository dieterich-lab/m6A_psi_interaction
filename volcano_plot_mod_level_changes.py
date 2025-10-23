import os
import sys
import re
import pandas as pd
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from tqdm import tqdm
from Bio import SeqIO
import logomaker
from argparse import ArgumentParser


def get_filtered_df_by_motifs(in_df, in_ref, in_motifs):
    df_out = []
    for _, this_row in in_df.iterrows():
        this_motif = in_ref[str(this_row['chrom'])][(this_row['start'] - 2):(this_row['start'] + 3)]
        if this_row['strand'] == '-':
            this_motif = this_motif.reverse_complement()
        if this_motif in in_motifs:
            df_out.append(this_row)
    return pd.DataFrame(df_out)

def get_mod_effect_size(in_mods, this_chunk_valid):
    num_mods = len(in_mods)
    acols = [f'a_{mod}' for mod in in_mods]
    bcols = [f'b_{mod}' for mod in in_mods]
    this_chunk_valid[acols] = this_chunk_valid['a_mod_percentages'].str.split(',', n=num_mods, expand=True).replace(r'.*\:', '', regex=True).astype(float)
    this_chunk_valid[bcols] = this_chunk_valid['b_mod_percentages'].str.split(',', n=num_mods, expand=True).replace(r'.*\:', '', regex=True).astype(float)
    cols = []
    for mod in in_mods:
        this_chunk_valid[mod] = np.abs(this_chunk_valid[f'a_{mod}'] - this_chunk_valid[f'b_{mod}'])
        cols.append(mod)
    this_chunk_valid['mod'] = this_chunk_valid[cols].idxmax(axis=1)
    return this_chunk_valid
    
def get_vec_change_neg_log_pval(in_dmr_file, in_args, in_ref, in_motifs=None):
    change_pval_score_mod = []
    sites = []
    print(f'Parsing {in_dmr_file}')

    with open(in_dmr_file) as f:
        m = next(re.finditer('#([^\n]+)\n(?!#)', f.readline()), None)
        names = m.group(1).split() if m else None
        in_mods = f.readline().split('\t')[6].split(',')
        in_mods = [mod.split(':')[0] for mod in in_mods]
    csv_args = {'sep': '\t', 'dtype': {'chrom': str}}
    if names is not None:
        csv_args.update({'names': names, 'comment': '#'})

    pval_col = 'balanced_map_pvalue'
    final_cols = ['balanced_effect_size', 'balanced_map_pvalue', 'score', 'mod']
    if not in_args.balanced:
        pval_col = 'p_value'
        final_cols = ['effect_size', 'p_value', 'score', 'mod']
    
    for this_chunk in tqdm(pd.read_csv(in_dmr_file, **csv_args, iterator=True, chunksize=10000)):
        if in_args.balanced and this_chunk[['balanced_map_pvalue', 'balanced_effect_size']].isna().all(axis=None):
            print('Balanced values not found. Check that DMR input contains replicates!')
            sys.exit(0)
        this_chunk_valid = this_chunk[
            (this_chunk[pval_col] < in_args.thresh_pval)
            * (this_chunk['a_total'] >= in_args.thresh_count)
            * (this_chunk['b_total'] >= in_args.thresh_count)
            * (1 - np.exp(-1*this_chunk['score']) > in_args.thresh_score)
            ].copy()
        if in_motifs:
            this_chunk_valid = get_filtered_df_by_motifs(this_chunk_valid, in_ref, in_motifs)
        if len(this_chunk_valid):
            this_chunk_valid = get_mod_effect_size(in_mods, this_chunk_valid)
            change_pval_score_mod.append(this_chunk_valid[final_cols].values)
            sites.append(this_chunk_valid[['chrom', 'start', 'end', 'mod', 'score', 'strand']])
            
    vec_effect_size, vec_pval, vec_score, vec_mod = np.vstack(change_pval_score_mod).T
    vec_effect_size= vec_effect_size.astype(float)
    vec_pval = vec_pval.astype(float)
    vec_score = vec_score.astype(float)
    vec_change = -vec_effect_size * 100
    vec_neg_log_pval = -np.log10(vec_pval+10**(-20))
    df_sites = pd.concat(sites)
    return vec_change, vec_neg_log_pval, vec_score, vec_mod, df_sites


def get_mask(in_vec_change, in_vec_neg_log_pval, in_vec_mod, in_args):
    mask_no_change = (in_vec_change < in_args.thresh_change) * (in_vec_change >= -in_args.thresh_change) \
                     + (in_vec_neg_log_pval < in_args.thresh_neg_log_pval)
    mask_pos_change = (in_vec_neg_log_pval >= in_args.thresh_neg_log_pval) * (in_vec_change >= in_args.thresh_change)
    mask_neg_change = (in_vec_neg_log_pval >= in_args.thresh_neg_log_pval) * (in_vec_change < -in_args.thresh_change)
    num_pos_change = np.sum(mask_pos_change)
    num_neg_change = np.sum(mask_neg_change)
    mask_change = dict()
    mask_change['pos_change'] = mask_pos_change
    mask_change['neg_change'] = mask_neg_change
    for mod in np.unique(in_vec_mod):
        mask_change[f'{mod}_pos_change'] = (mask_pos_change) * (in_vec_mod == mod)
        mask_change[f'{mod}_neg_change'] = (mask_neg_change) * (in_vec_mod == mod)
    return mask_no_change, mask_change, num_neg_change, num_pos_change


def get_ref(ref_file):
    if ref_file is None:
        print('Need to suppy reference!')
        sys.exit(0)
    print(f'Loading reference {ref_file}...')
    out_ref = {}
    for record in SeqIO.parse(ref_file, "fasta"):
        out_ref[record.id] = record.seq
    return out_ref


def plot_logo(in_df_sites, in_mask, in_mod, in_change, in_ref, in_args, span=3):
    sites_change = in_df_sites[in_mask]
    sites_change.to_csv(os.path.join(in_args.out_dir, f'sites_{in_change}_{in_mod}.bed'), sep='\t', index=False, header=False)
    all_refseq = []
    for _, row in sites_change.iterrows():
        refseq = in_ref[row['chrom']][(row['start'] - span):(row['start'] + span + 1)]
        if row['strand'] == '-':
            refseq = refseq.reverse_complement()
        all_refseq.append(refseq)
    df_logo = logomaker.alignment_to_matrix([str(seq) for seq in all_refseq])

    fig = plt.figure(figsize=(10, 5))
    ax = fig.subplots()
    logomaker.Logo(df_logo, ax=ax)
    ax.set_xticks(np.arange(0, 2 * span + 1, span), np.arange(0, 2 * span + 1, span) - span)
    ax.set_title(f'{in_mod}, {in_change}')
    fig.savefig(os.path.join(in_args.out_dir, f'logo_{in_change}_{in_mod}_site_motifs.png'), bbox_inches='tight')
    plt.close(fig)


def main():
    home = os.environ['HOME']
    parser = ArgumentParser()
    parser.add_argument('--dmr_base_A', type=str, required=True,
                        help='modkit DMR output for base A')
    parser.add_argument('--dmr_base_U', type=str, required=True,
                        help='modkit DMR output for base U')
    parser.add_argument('--reference', type=str, required=True,
                        help='reference fasta, must be the same one used for read alignment')
    parser.add_argument('--out_dir', type=str, default=home,
                        help='output directory')
    parser.add_argument('--plot_name', type=str, default='volcano_plot.png',
                        help='name of volcano plot')
    parser.add_argument('--thresh_neg_log_pval', type=float, default=2.0,
                        help='threshold for negative log p-value')
    parser.add_argument('--thresh_pval', type=float, default=1.0,
                        help='MAP p-value cutoff')
    parser.add_argument('--thresh_score', type=float, default=0.5,
                        help='Score cutoff as P(different)')
    parser.add_argument('--thresh_count', type=int, default=50,
                        help='threshold for site coverage')
    parser.add_argument('--thresh_change', type=float, default=25.0,
                        help='threshold for effect size')
    parser.add_argument('--balanced', action='store_true',
                        help='use balanced effect size and p-value, only available with replicates')
    parser.add_argument('--filter_by_motifs', action='store_true',
                        help='filter sites by pre-defined m6A and psi motifs')
    parser.add_argument('--xlim', type=int, default=101,
                        help='x-axis limit for display')
    
    args = parser.parse_args()

    ref = get_ref(args.reference)
    if args.filter_by_motifs:
        mod_motifs = {
            'A': [
                    'GGACT', 'GGACA', 'GAACT', 'AGACT', 'GGACC', 'TGACT',
                    'AAACT', 'GAACA', 'AGACA', 'AGACC', 'GAACC', 'TGACA',
                    'TAACT', 'AAACA', 'TGACC', 'TAACA', 'AAACC', 'TAACC'
            ],
            'U': ['GTTCA', 'GTTCC', 'GTTCG', 'GTTCT'] + ['TGTAG'] +
                   [f'{pos1}{pos2}T{pos4}{pos5}'
                    for pos1 in ['A', 'C', 'G', 'T']
                    for pos2 in ['A', 'G']
                    for pos4 in ['A', 'G']
                    for pos5 in ['A', 'C', 'G', 'T']
                    ]
        }
    else:
        mod_motifs = {
            'A': None,
            'U': None
        }

    os.makedirs(args.out_dir, exist_ok=True)

    dict_display_mod = {
        'A': r'Base\ A',
        'U': r'Base\ U'
    }
    dict_display_labels = {
        'a': 'm6A',
        '69426': 'Am',
        '17596': 'I',
        '17802': 'Y',
        '19227': 'Um',
        '19228': 'Cm',
        'm': 'm5C',
        '19229': 'Gm',
    }
    markers = ['o', 's', 'D']

    plt.figure(figsize=(10, 5))
    for subplot_ind, base in enumerate(dict_display_mod.keys()):
        this_dmr_file = args.__dict__.get(f'dmr_base_{base}')
        this_vec_change, this_vec_neg_log_pval, this_vec_score, this_vec_mod, this_df_sites = get_vec_change_neg_log_pval(this_dmr_file, args, ref, mod_motifs[base])
        this_df_sites['mod'] = this_df_sites['mod'].map(dict_display_labels)
        this_mask_no_change, this_mask_change, this_num_neg_change, this_num_pos_change = get_mask(this_vec_change, this_vec_neg_log_pval, this_vec_mod, args)
        # for this_change, this_mask in zip(['up', 'down'], [this_mask_change['pos_change'], this_mask_change['neg_change']]):
        #     if this_mask.any():
        #         plot_logo(this_df_sites, this_mask, base, this_change, ref, args)

        plt.subplot(1, 2, subplot_ind+1)
        plt.scatter(this_vec_change[this_mask_no_change], this_vec_neg_log_pval[this_mask_no_change], s=1, c='gray')
        for idx, mod in enumerate(np.unique(this_vec_mod)):
            plt.scatter(
                this_vec_change[this_mask_change[f'{mod}_pos_change']], 
                this_vec_neg_log_pval[this_mask_change[f'{mod}_pos_change']], 
                marker=markers[idx], 
                s=8, 
                facecolors='none', 
                edgecolors='red', 
                label=dict_display_labels[mod]
            )
            plt.scatter(
                this_vec_change[this_mask_change[f'{mod}_neg_change']], 
                this_vec_neg_log_pval[this_mask_change[f'{mod}_neg_change']], 
                marker=markers[idx], 
                s=8, 
                facecolors='none', 
                edgecolors='blue'
            )
        legend = plt.legend(frameon=False)
        handles = legend.legend_handles
        for i, handle in enumerate(handles):
            handle.set_edgecolor("#000000")
            handle.set_facecolor("#FFFFFF")
        plt.axhline(y=args.thresh_neg_log_pval, c='gray', ls='--')
        plt.axvline(x=-args.thresh_change, c='gray', ls='--')
        plt.axvline(x=args.thresh_change, c='gray', ls='--')
        plt.xlabel(f'% Mod. level change')
        plt.ylabel('$-log_{10}$ p-val')
        plt.title(f'${dict_display_mod[base]}$')
        plt.xlim([-args.xlim, args.xlim])
        plt.ylim([0, plt.gca().get_ylim()[1]])
        plt.text(0.01, 1.01, f'{this_num_neg_change} down', c='blue', ha='left', transform=plt.gca().transAxes)
        plt.text(0.99, 1.01, f'{this_num_pos_change} up', c='red', ha='right', transform=plt.gca().transAxes)
        # plt.savefig(os.path.join(args.out_dir, args.plot_name), bbox_inches='tight')
        plt.savefig(os.path.join(args.out_dir, f"{args.plot_name}.pdf"), format="pdf", bbox_inches='tight')
        plt.savefig(os.path.join(args.out_dir, f"{args.plot_name}.svg"), format="svg", bbox_inches='tight')

if __name__ == '__main__':
    main()
    print('Finished')
