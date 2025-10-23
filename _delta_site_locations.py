from Bio import SeqIO
from tqdm import tqdm
import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import os
import pickle
from functools import reduce
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


def get_merged_df_with_delta(in_dfs, conditions, writer):
  
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

    mask = ((out_df_merged.loc[:, out_df_merged.columns.str.contains('freq_')] > 0).any(axis=1))
    out_df_merged = out_df_merged[mask]

    if f'{writer}-KD' in conditions:
        out_df_merged[f'delta_{writer}-KD'] = out_df_merged[f'freq_{writer}-KD'] - out_df_merged['freq_CTRL']
    if f'{writer}-OE' in conditions:
        out_df_merged[f'delta_{writer}-OE'] = out_df_merged[f'freq_{writer}-OE'] - out_df_merged['freq_CTRL']

    return out_df_merged


def filter_by_delta_threshold(merged_df, conditions, writer, mod_type, thresh_delta):
    
    if f'{writer}-KD' in conditions:
        if mod_type == '17802':  # psi should decrease in KD
            df_filtered = merged_df[merged_df[f'delta_{writer}-KD'] < -thresh_delta]
        elif mod_type == 'a':  # m6A should increase in KD
            df_filtered = merged_df[merged_df[f'delta_{writer}-KD'] >= thresh_delta]
    elif f'{writer}-OE' in conditions:
        if mod_type == '17802':  # psi should increase in OE
            df_filtered = merged_df[merged_df[f'delta_{writer}-OE'] >= thresh_delta]
        elif mod_type == 'a':  # m6A should decrease in OE
            df_filtered = merged_df[merged_df[f'delta_{writer}-OE'] < -thresh_delta]
    else:
        print(f'Warning: No KD or OE condition found for {writer}')
        df_filtered = merged_df
    
    return df_filtered


def main():
    parser = ArgumentParser(description='Identify modification sites with significant delta changes between conditions')
    
    parser.add_argument('--pickle_file', type=str, required=True,
                        help='Pickle file with filtered modification dataframes')
    parser.add_argument('--output_dir', '-o', type=str, required=True,
                        help='Output directory for results')
    parser.add_argument('--writer', type=str, default='TRUB1',
                        help='Writer enzyme name (default: TRUB1)')
    parser.add_argument('--condition', type=str, required=True, choices=['KD', 'OE'],
                        help='Condition to analyze (KD or OE)')
    parser.add_argument('--thresh_delta', type=float, default=5.0,
                        help='Threshold for delta S filtering (default: 5.0)')
    parser.add_argument('--thresh_freq', type=float, default=None,
                        help='Frequency threshold for final m6A filtering (optional)')
    parser.add_argument('--output_prefix', type=str, default=None,
                        help='Prefix for output files (default: {writer}-{condition})')
    parser.add_argument('--save_bed', action='store_true',
                        help='Save BED file for m6A sites (first 6 columns)')
    parser.add_argument('--dpi', type=int, default=1200,
                        help='Resolution for matplotlib (default: 1200)')
    parser.add_argument('--font_size', type=int, default=8,
                        help='Font size for labels (default: 8)')
    
    args = parser.parse_args()
    
    configure_matplotlib(dpi=args.dpi, font_size=args.font_size)
    
    os.makedirs(args.output_dir, exist_ok=True)
    
    if args.condition == 'KD':
        conditions = [f'{args.writer}-KD', 'CTRL']
    else:  # OE
        conditions = [f'{args.writer}-OE', 'CTRL']
    
    if args.output_prefix is None:
        output_prefix = f'{args.writer}-{args.condition}'
    else:
        output_prefix = args.output_prefix
    
    print(f'Loading data from: {args.pickle_file}')
    with open(args.pickle_file, 'rb') as pkl_in:
        dfs_mod_filtered = pickle.load(pkl_in)
    
    mod_names = ['17802', 'a']  # psi and m6A
    dict_mod_display = {
        'a': 'm6A',
        '17802': 'ψ'
    }
    
    dfs_mod_thresh_delta = {}
    
    for this_mod in mod_names:
        print(f'\nProcessing modification: {dict_mod_display[this_mod]}')
        
        this_merged_df = get_merged_df_with_delta(
            dfs_mod_filtered[this_mod], 
            conditions, 
            args.writer
        )
        print(f'  Total merged sites: {len(this_merged_df)}')
        
        df_thresh_delta = filter_by_delta_threshold(
            this_merged_df, 
            conditions, 
            args.writer, 
            this_mod, 
            args.thresh_delta
        )
        print(f'  Sites passing delta threshold (|Δ| > {args.thresh_delta}): {len(df_thresh_delta)}')
        
        dfs_mod_thresh_delta[this_mod] = df_thresh_delta
    
    if args.thresh_freq is not None:
        print(f'\nApplying frequency threshold for m6A: >= {args.thresh_freq}')
        freq_col = f'freq_{args.writer}-{args.condition}'
        out_df = dfs_mod_thresh_delta['a'][dfs_mod_thresh_delta['a'][freq_col] >= args.thresh_freq]
        print(f'  m6A sites after frequency filtering: {len(out_df)}')
    else:
        out_df = dfs_mod_thresh_delta['a']
    
    output_file_m6a = os.path.join(args.output_dir, f'{output_prefix}_delta_sites_a.tsv')
    out_df.to_csv(output_file_m6a, sep='\t', index=False)
    print(f'\nSaved m6A delta sites: {output_file_m6a}')
    
    output_file_psi = os.path.join(args.output_dir, f'{output_prefix}_delta_sites_psi.tsv')
    dfs_mod_thresh_delta['17802'].to_csv(output_file_psi, sep='\t', index=False)
    print(f'Saved ψ delta sites: {output_file_psi}')
    
    if args.save_bed:
        out_bed = out_df.iloc[:, :6]
        bed_file = os.path.join(args.output_dir, f'{output_prefix}_m6a.bed')
        out_bed.to_csv(bed_file, sep='\t', index=False, header=False)
        print(f'Saved m6A BED file: {bed_file}')
    
    print('\n' + '='*60)
    print('SUMMARY')
    print('='*60)
    print(f'Writer: {args.writer}')
    print(f'Condition: {args.condition}')
    print(f'Delta threshold: {args.thresh_delta}')
    if args.thresh_freq is not None:
        print(f'Frequency threshold: {args.thresh_freq}')
    print(f'\nψ sites with significant change: {len(dfs_mod_thresh_delta["17802"])}')
    print(f'm6A sites with significant change: {len(out_df)}')
    print('='*60)
    
    print('\nFinished')


if __name__ == '__main__':
    main()