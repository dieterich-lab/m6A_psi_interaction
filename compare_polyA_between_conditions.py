import os
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from argparse import ArgumentParser
from tqdm import tqdm


def configure_matplotlib(dpi=300, font_size=12):
    matplotlib.rcParams['figure.dpi'] = dpi
    matplotlib.rcParams['savefig.dpi'] = dpi
    matplotlib.rcParams['font.size'] = font_size
    matplotlib.rcParams['legend.fontsize'] = font_size - 2
    matplotlib.rcParams['xtick.labelsize'] = font_size
    matplotlib.rcParams['ytick.labelsize'] = font_size
    matplotlib.rcParams['xtick.major.size'] = 4
    matplotlib.rcParams['ytick.major.size'] = 4
    matplotlib.rcParams['lines.linewidth'] = 1


def load_polyA_data(file_path):
  
    print(f'Loading poly(A) data from: {file_path}')
    df = pd.read_csv(file_path, sep='\t', names=['read_id', 'polyA_len'])
    print(f'  Loaded {len(df)} reads')
    polyA_dict = {k: v for k, v in df[['read_id', 'polyA_len']].values}
    return polyA_dict


def plot_polyA_distribution(dict_polyA, conditions, condition_labels, condition_colors,
                            xlim, num_bins, output_path, figsize=(4, 4),
                            fmt='png', dpi=300, transparent=False):
    
    fig_kwargs = dict(format=fmt, bbox_inches='tight', dpi=dpi, transparent=transparent)
    
    bin_edges = np.linspace(*xlim, num_bins+1)
    bin_centers = 0.5 * (bin_edges[1:] + bin_edges[:-1])
    xticks = np.linspace(*xlim, 4)
    
    plt.figure(figsize=figsize)
    
    for this_cond in conditions:
        this_hist, _ = np.histogram(list(dict_polyA[this_cond].values()), bins=bin_edges)
        num_reads = this_hist.sum()
        norm_hist = this_hist / np.sum(this_hist)
        label = condition_labels.get(this_cond, this_cond)
        plt.plot(bin_centers, norm_hist, c=condition_colors[this_cond], 
                label=f'{label} ({num_reads})')
    
    plt.legend(fontsize=10)
    plt.xlim(xlim)
    plt.xticks(xticks)
    plt.xlabel('poly(A) length (bp)', fontsize=12)
    plt.ylabel('Density', fontsize=12)
    
    plt.savefig(output_path, **fig_kwargs)
    plt.close()
    print(f'Saved plot: {output_path}')


def compute_statistics(dict_polyA, conditions):
    
    print('\nPoly(A) length statistics:')
    for this_cond in conditions:
        lengths = list(dict_polyA[this_cond].values())
        print(f'\n  {this_cond}:')
        print(f'    Count: {len(lengths)}')
        print(f'    Mean: {np.mean(lengths):.2f} bp')
        print(f'    Median: {np.median(lengths):.2f} bp')
        print(f'    Std: {np.std(lengths):.2f} bp')
        print(f'    Min: {np.min(lengths)} bp')
        print(f'    Max: {np.max(lengths)} bp')


def main():
    parser = ArgumentParser(description='Compare poly(A) tail length distributions between conditions')
    
    parser.add_argument('--polyA_files', type=str, nargs='+', required=True,
                        help='Paths to poly(A) length files (TSV format with read_id and polyA_len columns)')
    parser.add_argument('--condition_names', type=str, nargs='+', required=True,
                        help='Names for each condition (must match number of files)')
    parser.add_argument('--output_file', '-o', type=str, required=True,
                        help='Output file path for histogram plot')
    parser.add_argument('--condition_labels', type=str, nargs='*', default=None,
                        help='Display labels for conditions (optional, defaults to condition_names)')
    parser.add_argument('--condition_colors', type=str, nargs='*', default=None,
                        help='Colors for each condition (optional, defaults to automatic colors)')
    parser.add_argument('--xlim', type=float, nargs=2, default=[0.0, 300.0],
                        help='X-axis limits (min max) (default: 0 300)')
    parser.add_argument('--num_bins', type=int, default=100,
                        help='Number of histogram bins (default: 100)')
    parser.add_argument('--figsize', type=float, nargs=2, default=[4.0, 4.0],
                        help='Figure size (width height) (default: 4 4)')
    parser.add_argument('--format', type=str, default='png', choices=['png', 'pdf', 'svg'],
                        help='Output format (default: png)')
    parser.add_argument('--dpi', type=int, default=300,
                        help='Resolution for output (default: 300)')
    parser.add_argument('--font_size', type=int, default=12,
                        help='Font size for labels (default: 12)')
    parser.add_argument('--transparent', action='store_true',
                        help='Save with transparent background')
    parser.add_argument('--show_stats', action='store_true',
                        help='Print poly(A) length statistics')
    
    args = parser.parse_args()
    
    if len(args.polyA_files) != len(args.condition_names):
        print('Error: Number of files must match number of condition names')
        return
    
    conditions = args.condition_names
    
    if args.condition_labels is None:
        condition_labels = {cond: cond for cond in conditions}
    else:
        if len(args.condition_labels) != len(conditions):
            print('Error: Number of labels must match number of conditions')
            return
        condition_labels = {cond: label for cond, label in zip(conditions, args.condition_labels)}
    
    if args.condition_colors is None:
        default_colors = ['blue', 'red', 'green', 'orange', 'purple', 'brown', 'pink', 'gray']
        condition_colors = {cond: default_colors[i % len(default_colors)] 
                          for i, cond in enumerate(conditions)}
    else:
        if len(args.condition_colors) != len(conditions):
            print('Error: Number of colors must match number of conditions')
            return
        condition_colors = {cond: color for cond, color in zip(conditions, args.condition_colors)}
    
    configure_matplotlib(dpi=args.dpi, font_size=args.font_size)
    
    print('Loading poly(A) data...')
    dict_polyA = {}
    for this_cond, this_file in zip(conditions, args.polyA_files):
        dict_polyA[this_cond] = load_polyA_data(this_file)
    
    if args.show_stats:
        compute_statistics(dict_polyA, conditions)
    
    print('\nGenerating plot...')
    plot_polyA_distribution(
        dict_polyA,
        conditions,
        condition_labels,
        condition_colors,
        args.xlim,
        args.num_bins,
        args.output_file,
        figsize=tuple(args.figsize),
        fmt=args.format,
        dpi=args.dpi,
        transparent=args.transparent
    )
    
    print('\nFinished')


if __name__ == '__main__':
    main()