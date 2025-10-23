from Bio import SeqIO
from tqdm import tqdm
import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import os
import pickle
from argparse import ArgumentParser


def configure_matplotlib(dpi=1200, font_size=8):
    mpl.rcParams['figure.dpi'] = dpi
    mpl.rcParams['savefig.dpi'] = dpi
    mpl.rcParams['font.size'] = font_size
    mpl.rcParams['legend.fontsize'] = font_size - 2
    mpl.rcParams['xtick.labelsize'] = font_size
    mpl.rcParams['ytick.labelsize'] = font_size
    mpl.rcParams['xtick.major.size'] = 4
    mpl.rcParams['ytick.major.size'] = 4
    mpl.rcParams['lines.linewidth'] = 1
    mpl.rcParams['font.family'] = 'Arial'


def get_ref(ref_file):
    print(f'Loading reference genome from: {ref_file}')
    out_ref = {}
    for record in SeqIO.parse(ref_file, "fasta"):
        out_ref[record.id] = record.seq
    return out_ref


def get_central_motif(in_df, in_ref, span=2):
    print(f'\nAssigning motif to {len(in_df)} sites:')
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
    bed_file = os.path.join(in_base_dir, in_ds, f'{in_sample}.modkit042.cov10.bedmethyl')
    print(f'Reading: {bed_file}')
    out_df = pd.read_csv(bed_file, sep='\t',
                         usecols=[0, 1, 2, 3, 4, 5, 10], names=bed_fields, dtype={'chrom': str})
    out_df = out_df[out_df['score'] >= thresh_cov]
    print(f'  Sites with coverage >= {thresh_cov}: {len(out_df)}')
    return out_df


def get_df_mod_filtered(in_df, in_mod, in_ref, in_motifs):
    this_df_mod = in_df[in_df['name'] == in_mod]
    this_df_mod = get_central_motif(this_df_mod, in_ref)
    return this_df_mod[this_df_mod['ref_motif'].isin(in_motifs)]


def get_motif_definitions(writer='TRUB1'):
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
            'PUS7': ['TGTAG'],
            'PUS1': [f'{pos1}{pos2}T{pos4}{pos5}'
                     for pos1 in ['A', 'C', 'G', 'T']
                     for pos2 in ['A', 'G']
                     for pos4 in ['A', 'G']
                     for pos5 in ['A', 'C', 'G', 'T']]
        }
    }
    
    mod_motif_name = {
        '17802': writer,
        'a': 'DRACH'
    }
    
    return sel_motifs, mod_motif_name


def plot_histograms(dfs_mod_filtered, conditions, writer, output_dir, 
                   fmt='png', dpi=1200, num_bins=5, transparent=False):
    cm = 1/2.54  # centimeters in inches
    fig_kwargs = dict(format=fmt, bbox_inches='tight', dpi=dpi, transparent=transparent)
    
    mod_names = ['17802', 'a']
    dict_mod_display = {
        'a': 'm^6A',
        '17802': '\psi'
    }
    
    cond_colors = {
        'CTRL': 'gray',
        f'{writer}-OE': 'red',
        f'{writer}-KD': 'blue'
    }
    
    _, mod_motif_name = get_motif_definitions(writer)
    
    for sel_mod in mod_names:
        plt.figure(figsize=(5*cm, 5*cm))
        for this_cond in conditions:
            this_df_mod_filtered = dfs_mod_filtered[sel_mod][this_cond]
            this_hist, bin_edges = np.histogram(this_df_mod_filtered['frequency'], range=[0, 100], bins=num_bins)
            this_hist_norm = this_hist / np.sum(this_hist)
            bin_centers = 0.5 * (bin_edges[1:] + bin_edges[:-1])
            plt.semilogy(bin_centers, this_hist_norm, label=f'{this_cond}', c=cond_colors[this_cond])
        
        plt.legend()
        plt.xticks(np.linspace(0, 100, num_bins+1))
        plt.xlabel('S')
        plt.ylabel('Probability (log)')
        plt.title(f'${dict_mod_display[sel_mod]}$, {mod_motif_name[sel_mod]} motifs')
        
        output_file = os.path.join(output_dir, f'histogram_mod_{sel_mod}_{writer}.{fmt}')
        plt.savefig(output_file, **fig_kwargs)
        plt.close()
        print(f'Saved: {output_file}')


def main():
    parser = ArgumentParser(description='Generate histograms of modification frequency (S) distribution')
    
    parser.add_argument('--base_dir', type=str, required=True,
                        help='Base directory containing bedmethyl files')
    parser.add_argument('--reference', type=str, required=True,
                        help='Reference genome FASTA file')
    parser.add_argument('--output_dir', '-o', type=str, required=True,
                        help='Output directory for plots')
    parser.add_argument('--writer', type=str, default='TRUB1',
                        help='Writer enzyme name (default: TRUB1)')
    parser.add_argument('--ctrl_dataset', type=str, default='HEK293_psU-KD',
                        help='Dataset name for control (default: HEK293_psU-KD)')
    parser.add_argument('--kd_dataset', type=str, default='HEK293_psU-KD',
                        help='Dataset name for knockdown (default: HEK293_psU-KD)')
    parser.add_argument('--oe_dataset', type=str, default='HEK293_psU-OE',
                        help='Dataset name for overexpression (default: HEK293_psU-OE)')
    parser.add_argument('--ctrl_sample', type=str, default='CTRL',
                        help='Sample name for control (default: CTRL)')
    parser.add_argument('--kd_sample', type=str, default=None,
                        help='Sample name for knockdown (default: {writer}-KD)')
    parser.add_argument('--oe_sample', type=str, default=None,
                        help='Sample name for overexpression (default: {writer}-OE)')
    parser.add_argument('--thresh_cov', type=int, default=100,
                        help='Coverage threshold for filtering sites (default: 100)')
    parser.add_argument('--num_bins', type=int, default=5,
                        help='Number of histogram bins (default: 5)')
    parser.add_argument('--format', type=str, default='png', choices=['png', 'pdf', 'svg'],
                        help='Output format (default: png)')
    parser.add_argument('--dpi', type=int, default=1200,
                        help='Resolution for output (default: 1200)')
    parser.add_argument('--font_size', type=int, default=8,
                        help='Font size for labels (default: 8)')
    parser.add_argument('--transparent', action='store_true',
                        help='Save with transparent background')
    parser.add_argument('--cache_file', type=str, default=None,
                        help='Pickle file to cache/load filtered dataframes')
    parser.add_argument('--force_recompute', action='store_true',
                        help='Force recomputation even if cache exists')
    
    args = parser.parse_args()
    
    configure_matplotlib(dpi=args.dpi, font_size=args.font_size)
    
    kd_sample = args.kd_sample if args.kd_sample else f'{args.writer}-KD'
    oe_sample = args.oe_sample if args.oe_sample else f'{args.writer}-OE'
    
    conditions = ['CTRL', kd_sample, oe_sample]
    
    os.makedirs(args.output_dir, exist_ok=True)
    
    if args.cache_file is None:
        cache_file = os.path.join(args.base_dir, f'dfs_mod_filtered_{args.writer}.pkl')
    else:
        cache_file = args.cache_file
    
    if os.path.exists(cache_file) and not args.force_recompute:
        print(f'Loading cached data from: {cache_file}')
        with open(cache_file, 'rb') as pkl_in:
            dfs_mod_filtered = pickle.load(pkl_in)
    else:
        print('Computing filtered dataframes...')
        
        dfs = {}
        dfs['CTRL'] = get_df_from_bed(args.base_dir, args.ctrl_dataset, args.ctrl_sample, args.thresh_cov)
        dfs[kd_sample] = get_df_from_bed(args.base_dir, args.kd_dataset, kd_sample, args.thresh_cov)
        dfs[oe_sample] = get_df_from_bed(args.base_dir, args.oe_dataset, oe_sample, args.thresh_cov)
        
        ref = get_ref(args.reference)
        
        sel_motifs, mod_motif_name = get_motif_definitions(args.writer)
        
        mod_names = ['17802', 'a']
        dfs_mod_filtered = {this_mod: {} for this_mod in mod_names}
        
        for sel_mod in mod_names:
            print(f'\nProcessing modification: {sel_mod}')
            for this_cond in conditions:
                print(f'  Condition: {this_cond}')
                this_df = dfs[this_cond]
                this_df_mod_filtered = get_df_mod_filtered(
                    this_df, sel_mod, ref, 
                    sel_motifs[sel_mod][mod_motif_name[sel_mod]]
                )
                dfs_mod_filtered[sel_mod][this_cond] = this_df_mod_filtered
                print(f'    Filtered sites: {len(this_df_mod_filtered)}')
        
        print(f'\nSaving cache to: {cache_file}')
        with open(cache_file, 'wb') as pkl_out:
            pickle.dump(dfs_mod_filtered, pkl_out)
    
    print('\nGenerating histograms...')
    plot_histograms(
        dfs_mod_filtered, 
        conditions, 
        args.writer, 
        args.output_dir,
        fmt=args.format,
        dpi=args.dpi,
        num_bins=args.num_bins,
        transparent=args.transparent
    )
    
    print('\nFinished')


if __name__ == '__main__':
    main()