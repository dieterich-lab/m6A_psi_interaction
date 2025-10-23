import os
import pandas as pd
from collections import Counter
from tqdm import tqdm
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from argparse import ArgumentParser
pd.set_option('display.max_columns', None)


def configure_matplotlib(dpi=300, font_size=10):
    matplotlib.rcParams['figure.dpi'] = dpi
    matplotlib.rcParams['savefig.dpi'] = dpi
    matplotlib.rcParams['font.size'] = font_size
    matplotlib.rcParams['legend.fontsize'] = font_size - 2
    matplotlib.rcParams['xtick.labelsize'] = font_size
    matplotlib.rcParams['ytick.labelsize'] = font_size
    matplotlib.rcParams['xtick.major.size'] = 4
    matplotlib.rcParams['ytick.major.size'] = 4
    matplotlib.rcParams['lines.linewidth'] = 1


def load_and_filter_transcripts(file_path, min_cov=0, class_code='='):
    print(f'Loading: {file_path}')
    df = pd.read_csv(file_path, sep='\t')
    df = df[(df['cov'] > min_cov) & (df['class_code'] == class_code)]
    print(f'  Transcripts after filtering: {len(df)}')
    
    valid_genes = [gene for gene, num_transcripts in Counter(df['ref_gene_id']).items() 
                   if num_transcripts > 1]
    print(f'  Genes with multiple transcripts: {len(valid_genes)}')
    
    return df, valid_genes


def merge_and_compute_ratios(df_ctrl, df_kd, common_genes, keep_fields):
    print(f'\nMerging data for {len(common_genes)} genes...')
    dfs_merged = []
    
    for this_gene in tqdm(common_genes):
        sub_df_ctrl = df_ctrl[df_ctrl['ref_gene_id'] == this_gene][keep_fields]
        sub_df_kd = df_kd[df_kd['ref_gene_id'] == this_gene][keep_fields]
        sub_df_merged = pd.merge(sub_df_ctrl, sub_df_kd, on=['ref_gene_id', 'ref_id'], 
                                 suffixes=['_ctrl', '_kd'])
        sub_df_merged['total_ctrl'] = sub_df_merged['cov_ctrl'].sum()
        sub_df_merged['ratio_ctrl'] = sub_df_merged['cov_ctrl'] / sub_df_merged['total_ctrl'] * 100.0
        sub_df_merged['total_kd'] = sub_df_merged['cov_kd'].sum()
        sub_df_merged['ratio_kd'] = sub_df_merged['cov_kd'] / sub_df_merged['cov_kd'].sum() * 100.0
        sub_df_merged['delta_ratio'] = sub_df_merged['ratio_kd'] - sub_df_merged['ratio_ctrl']
        sub_df_merged['max_abs_delta_ratio'] = sub_df_merged['delta_ratio'].abs().max()
        dfs_merged.append(sub_df_merged)
    
    df_merged = pd.concat(dfs_merged)
    print(f'Total merged entries: {len(df_merged)}')
    
    return df_merged


def plot_histogram(df_merged_filtered, thresh_total, output_dir, fmt='png', dpi=300):
    plt.figure(figsize=(5, 5))
    plt.hist(df_merged_filtered['max_abs_delta_ratio'], range=[0, 100], bins=20, log=True)
    plt.xlabel('Max. Abs. $\Delta$ isoform %')
    plt.ylabel('Counts')
    plt.title(f'Gene coverage $\geq$ {thresh_total}')
    
    output_file = os.path.join(output_dir, f'hist_delta_ratio.{fmt}')
    plt.savefig(output_file, bbox_inches='tight', dpi=dpi)
    plt.close()
    print(f'Saved histogram: {output_file}')


def plot_isoform_bars(df_sel, output_dir, ctrl_label='CTRL', kd_label='KD', 
                     fmt='png', dpi=300, figsize=(5, 5)):
    cmap = matplotlib.colormaps['Spectral']
    unique_genes = df_sel['ref_gene_id'].unique()
    
    print(f'\nGenerating isoform plots for {len(unique_genes)} genes...')
    
    for this_gene in tqdm(unique_genes):
        sub_df_sel = df_sel[df_sel['ref_gene_id'] == this_gene].copy()
        sub_df_sel.sort_values('ref_id', inplace=True)
        transcripts = sub_df_sel['ref_id']
        num_isoforms = len(transcripts)
        shifts = np.linspace(-0.15, 0.15, num_isoforms)
        bar_width = (0.5 / num_isoforms) * 0.75
        tx_colors = cmap(np.linspace(0, 1, num_isoforms))

        plt.figure(figsize=figsize)
        for tx_ind, this_transcript in enumerate(transcripts):
            vec_x = np.array([1, 2]) + shifts[tx_ind]
            vec_y = sub_df_sel[sub_df_sel['ref_id'] == this_transcript][['ratio_ctrl', 'ratio_kd']].to_numpy()[0]
            plt.bar(vec_x, vec_y, width=bar_width, color=tx_colors[tx_ind], label=this_transcript)
        
        total_ctrl = sub_df_sel['total_ctrl'].iloc[0]
        total_kd = sub_df_sel['total_kd'].iloc[0]
        plt.xticks([1, 2], [f"{ctrl_label} ({total_ctrl})", f"{kd_label} ({total_kd})"])
        plt.xlabel('Conditions (total num. reads)')
        plt.ylabel('Isoform ratio (%)')
        plt.legend(loc='upper center')
        plt.title(this_gene)
        
        output_file = os.path.join(output_dir, f'{this_gene}_isoforms.{fmt}')
        plt.savefig(output_file, bbox_inches='tight', dpi=dpi)
        plt.close()


def main():
    parser = ArgumentParser(description='Differential Transcript Usage (DTU) Analysis')
    
    parser.add_argument('--ctrl_file', type=str, required=True,
                        help='Path to control condition transcript table (TSV)')
    parser.add_argument('--kd_file', type=str, required=True,
                        help='Path to knockdown condition transcript table (TSV)')
    parser.add_argument('--output_dir', '-o', type=str, required=True,
                        help='Output directory for results and plots')
    parser.add_argument('--ctrl_label', type=str, default='CTRL',
                        help='Label for control condition (default: CTRL)')
    parser.add_argument('--kd_label', type=str, default='KD',
                        help='Label for knockdown condition (default: KD)')
    parser.add_argument('--thresh_total', type=int, default=20,
                        help='Minimum total gene coverage threshold (default: 20)')
    parser.add_argument('--thresh_delta', type=float, default=50.0,
                        help='Minimum absolute delta isoform ratio threshold (default: 50.0)')
    parser.add_argument('--min_cov', type=int, default=0,
                        help='Minimum coverage for initial filtering (default: 0)')
    parser.add_argument('--class_code', type=str, default='=',
                        help='Class code filter (default: =)')
    parser.add_argument('--format', type=str, default='png', choices=['png', 'pdf', 'svg'],
                        help='Output format for plots (default: png)')
    parser.add_argument('--dpi', type=int, default=300,
                        help='Resolution for output plots (default: 300)')
    parser.add_argument('--font_size', type=int, default=10,
                        help='Font size for labels (default: 10)')
    parser.add_argument('--figsize', type=float, nargs=2, default=[5.0, 5.0],
                        help='Figure size (width height) (default: 5 5)')
    
    args = parser.parse_args()
    
    configure_matplotlib(dpi=args.dpi, font_size=args.font_size)
    
    os.makedirs(args.output_dir, exist_ok=True)
    
    df_ctrl, valid_genes_ctrl = load_and_filter_transcripts(
        args.ctrl_file, args.min_cov, args.class_code
    )
    df_kd, valid_genes_kd = load_and_filter_transcripts(
        args.kd_file, args.min_cov, args.class_code
    )
    
    common_valid_genes = list(set(valid_genes_ctrl).intersection(set(valid_genes_kd)))
    print(f'\nCommon genes with multiple transcripts: {len(common_valid_genes)}')
    
    keep_fields = [
        'ref_gene_id',
        'ref_id',
        'class_code',
        'num_exons',
        'cov',
        'len',
        'sample_id',
        'parent gene iso num'
    ]
    
    df_merged = merge_and_compute_ratios(df_ctrl, df_kd, common_valid_genes, keep_fields)
    
    df_merged_filtered = df_merged[
        (df_merged['total_ctrl'] >= args.thresh_total)
        & (df_merged['total_kd'] >= args.thresh_total)
    ]
    print(f'\nGenes passing coverage threshold (>= {args.thresh_total}): '
          f'{len(df_merged_filtered["ref_gene_id"].unique())}')
    
    plot_histogram(df_merged_filtered, args.thresh_total, args.output_dir, 
                   args.format, args.dpi)
    
    df_sel = df_merged_filtered[df_merged_filtered['max_abs_delta_ratio'] >= args.thresh_delta]
    df_sel.sort_values('ref_gene_id', inplace=True)
    
    print(f'\nGenes passing delta threshold (>= {args.thresh_delta}): '
          f'{len(df_sel["ref_gene_id"].unique())}')
    
    output_table = os.path.join(
        args.output_dir, 
        f'transcript_table_gene_count_{args.thresh_total}_delta_{args.thresh_delta}.tsv'
    )
    df_sel.to_csv(output_table, sep='\t', index=False, float_format='%.3f')
    print(f'Saved filtered results: {output_table}')
    
    if len(df_sel) > 0:
        plot_isoform_bars(
            df_sel, 
            args.output_dir, 
            ctrl_label=args.ctrl_label,
            kd_label=args.kd_label,
            fmt=args.format,
            dpi=args.dpi,
            figsize=tuple(args.figsize)
        )
    else:
        print('No genes passed the thresholds for plotting.')
    
    print('\nFinished')


if __name__ == '__main__':
    main()

