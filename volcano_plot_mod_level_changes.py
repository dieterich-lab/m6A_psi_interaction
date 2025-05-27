import os
import sys
import pandas as pd
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from tqdm import tqdm
from Bio import SeqIO
import logomaker
from argparse import ArgumentParser
mpl.use('TkAgg')


def get_filtered_df_by_motifs(in_df, in_ref, in_motifs):
    df_out = []
    for _, this_row in in_df.iterrows():
        this_motif = in_ref[str(this_row['chrom'])][(this_row['chromStart'] - 2):(this_row['chromStart'] + 3)]
        if this_row['strand'] == '-':
            this_motif = this_motif.reverse_complement()
        if this_motif in in_motifs:
            df_out.append(this_row)
    return pd.DataFrame(df_out)


def get_vec_change_neg_log_pval(in_dmr_file, in_args, in_ref, in_motifs=None):
    change_pval_score = []
    sites = []
    print(f'Parsing {in_dmr_file}')
    for this_chunk in tqdm(pd.read_csv(in_dmr_file, sep='\t', dtype={'chrom': str}, iterator=True, chunksize=10000)):
        if in_args.balanced:
            if this_chunk[['balanced_p_value', 'balanced_effect_size']].isna().all(axis=None):
                print('Balanced values not found. Check that DMR input contains replicates!')
                sys.exit(0)
            this_chunk_valid = this_chunk[
                (this_chunk['balanced_p_value'] < 1)
                * (this_chunk['sample_a_total'] >= in_args.thresh_count)
                * (this_chunk['sample_b_total'] >= in_args.thresh_count)
                ]
            if in_motifs:
                this_chunk_valid = get_filtered_df_by_motifs(this_chunk_valid, in_ref, in_motifs)
            if len(this_chunk_valid):
                change_pval_score.append(this_chunk_valid[['balanced_effect_size', 'balanced_p_value', 'score']].values)
        else:
            this_chunk_valid = this_chunk[
                (this_chunk['p_value'] < 1)
                * (this_chunk['sample_a_total'] >= in_args.thresh_count)
                * (this_chunk['sample_b_total'] >= in_args.thresh_count)
                ]
            if in_motifs:
                this_chunk_valid = get_filtered_df_by_motifs(this_chunk_valid, in_ref, in_motifs)
            if len(this_chunk_valid):
                change_pval_score.append(this_chunk_valid[['effect_size', 'p_value', 'score']].values)
        if len(this_chunk_valid):
            sites.append(this_chunk_valid[['chrom', 'chromStart', 'strand']])
    vec_effect_size, vec_pval, vec_score = np.vstack(change_pval_score).T
    vec_change = -vec_effect_size * 100
    vec_neg_log_pval = -np.log10(vec_pval+10**(-20))
    df_sites = pd.concat(sites)
    return vec_change, vec_neg_log_pval, vec_score, df_sites


def get_mask(in_vec_change, in_vec_neg_log_pval, in_args):
    mask_no_change = (in_vec_change < in_args.thresh_change) * (in_vec_change >= -in_args.thresh_change) \
                     + (in_vec_neg_log_pval < in_args.thresh_neg_log_pval)
    mask_pos_change = (in_vec_neg_log_pval >= in_args.thresh_neg_log_pval) * (in_vec_change >= in_args.thresh_change)
    mask_neg_change = (in_vec_neg_log_pval >= in_args.thresh_neg_log_pval) * (in_vec_change < -in_args.thresh_change)
    num_pos_change = np.sum(mask_pos_change)
    num_neg_change = np.sum(mask_neg_change)
    return mask_no_change, mask_neg_change, mask_pos_change, num_neg_change, num_pos_change


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
    sites_change.to_csv(os.path.join(in_args.out_dir, f'sites_{in_change}_{in_mod}.tsv'), sep='\t', index=False)
    all_refseq = []
    for _, row in sites_change.iterrows():
        refseq = in_ref[row['chrom']][(row['chromStart'] - span):(row['chromStart'] + span + 1)]
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
    parser.add_argument('--out_dir', type=str, default=home,
                        help='output directory')
    parser.add_argument('--plot_name', type=str, default='volcano_plot.png',
                        help='name of volcano plot')
    parser.add_argument('--thresh_neg_log_pval', type=float, default=2.0,
                        help='threshold for negative log p-value')
    parser.add_argument('--thresh_count', type=int, default=50,
                        help='threshold for site coverage')
    parser.add_argument('--thresh_change', type=float, default=25.0,
                        help='threshold for effect size')
    parser.add_argument('--reference', type=str, default=None,
                        help='reference fasta, must be the same one used for read alignment')
    parser.add_argument('--balanced', action='store_true',
                        help='use balanced effect size and p-value, only available with replicates')
    parser.add_argument('--filter_by_motifs', action='store_true',
                        help='filter sites by pre-defined m6A and psi motifs')
    parser.add_argument('--output_logo', action='store_true',
                        help='output sequence context of up / down sites and their locations')
    args = parser.parse_args()

    if args.filter_by_motifs or args.output_logo:
        ref = get_ref(args.reference)
    else:
        ref = None

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

    # dict_display_mod = {
    #     'A': 'm^6A',
    #     'U': '\Psi'
    # }
    dict_display_mod = {
        'A': r'Base\ A',
        'U': r'Base\ U'
    }

    plt.figure(figsize=(10, 5))
    for subplot_ind, base in enumerate(dict_display_mod.keys()):
        this_dmr_file = args.__dict__.get(f'dmr_base_{base}')
        this_vec_change, this_vec_neg_log_pval, this_vec_score, this_df_sites = get_vec_change_neg_log_pval(this_dmr_file, args, ref, mod_motifs[base])
        this_mask_no_change, this_mask_neg_change, this_mask_pos_change, this_num_neg_change, this_num_pos_change = get_mask(this_vec_change, this_vec_neg_log_pval, args)
        if args.output_logo:
            for this_change, this_mask in zip(['up', 'down'], [this_mask_pos_change, this_mask_neg_change]):
                if this_mask.any():
                    plot_logo(this_df_sites, this_mask, base, this_change, ref, args)

        plt.subplot(1, 2, subplot_ind+1)
        plt.scatter(this_vec_change[this_mask_no_change], this_vec_neg_log_pval[this_mask_no_change], s=1, c='gray')
        plt.scatter(this_vec_change[this_mask_pos_change], this_vec_neg_log_pval[this_mask_pos_change], s=3, c='red')
        plt.scatter(this_vec_change[this_mask_neg_change], this_vec_neg_log_pval[this_mask_neg_change], s=3, c='blue')
        plt.axhline(y=args.thresh_neg_log_pval, c='gray', ls='--')
        plt.axvline(x=-args.thresh_change, c='gray', ls='--')
        plt.axvline(x=args.thresh_change, c='gray', ls='--')
        plt.xlabel(f'% Mod. level change')
        plt.ylabel('$-log_{10}$ p-val')
        plt.title(f'${dict_display_mod[base]}$')
        plt.xlim([-101, 101])
        plt.ylim([0, plt.gca().get_ylim()[1]])
        plt.text(0.01, 1.01, f'{this_num_neg_change} down', c='blue', ha='left', transform=plt.gca().transAxes)
        plt.text(0.99, 1.01, f'{this_num_pos_change} up', c='red', ha='right', transform=plt.gca().transAxes)
        plt.savefig(os.path.join(args.out_dir, args.plot_name), bbox_inches='tight')


if __name__ == '__main__':
    main()
    print('Finished')