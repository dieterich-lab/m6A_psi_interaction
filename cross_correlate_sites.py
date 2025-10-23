import os
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from argparse import ArgumentParser
from tqdm import tqdm


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


def load_sites(file_path):
    print(f'Loading sites from: {file_path}')
    df = pd.read_csv(file_path, sep='\t', dtype={'chrom': str})
    print(f'  Loaded {len(df)} sites')
    return df


def calculate_nearest_distances(df_sites0, df_sites1, label0='sites0', label1='sites1'):

    print(f'\nCalculating nearest distances from {label0} to {label1}...')
    min_dist = []
    
    for _, this_row in tqdm(df_sites0.iterrows(), total=len(df_sites0)):
        chrom0, chromStart0, strand0 = this_row[['chrom', 'chromStart', 'strand']]
        sub_df_sites1 = df_sites1[
            (df_sites1['chrom'] == chrom0)
            & (df_sites1['strand'] == strand0)
        ]

        if len(sub_df_sites1) > 0:
            distances = (sub_df_sites1['chromStart'] - chromStart0).values
            if strand0 == '-':
                distances = -distances
            min_dist.append(distances[np.argmin(np.abs(distances))])
    
    print(f'  Found nearest neighbors for {len(min_dist)} / {len(df_sites0)} sites')
    return min_dist


def plot_distance_histogram(min_dist, cond0, cond1, dataset_label, output_path,
                            xmax=10, figsize=(4, 4), fmt='png', dpi=300, transparent=False):

    fig_kwargs = dict(format=fmt, bbox_inches='tight', dpi=dpi, transparent=transparent)
    
    plt.figure(figsize=figsize)
    counts, _, _ = plt.hist(min_dist, range=[-xmax, xmax], bins=2*xmax)
    plt.xlabel(f'Nearest distance (nt)\nPos({cond1}) - Pos({cond0})')
    plt.ylabel('Site count')
    plt.title(f'{int(np.sum(counts))} / {len(min_dist)} {cond0} sites\n{dataset_label}')
    
    plt.savefig(output_path, **fig_kwargs)
    plt.close()
    print(f'Saved plot: {output_path}')


def main():
    parser = ArgumentParser(description='Cross-correlate two sets of modification sites by calculating nearest distances')
    
    parser.add_argument('--sites_file0', type=str, required=True,
                        help='Path to first sites file (TSV format)')
    parser.add_argument('--sites_file1', type=str, required=True,
                        help='Path to second sites file (TSV format)')
    parser.add_argument('--output_file', '-o', type=str, required=True,
                        help='Output file path for histogram plot')
    parser.add_argument('--label0', type=str, default='sites0',
                        help='Label for first set of sites (default: sites0)')
    parser.add_argument('--label1', type=str, default='sites1',
                        help='Label for second set of sites (default: sites1)')
    parser.add_argument('--dataset_label', type=str, default='',
                        help='Label for dataset (default: empty)')
    parser.add_argument('--xmax', type=int, default=10,
                        help='Maximum distance for histogram range in nucleotides (default: 10)')
    parser.add_argument('--figsize', type=float, nargs=2, default=[4.0, 4.0],
                        help='Figure size (width height) (default: 4 4)')
    parser.add_argument('--format', type=str, default='png', choices=['png', 'pdf', 'svg'],
                        help='Output format (default: png)')
    parser.add_argument('--dpi', type=int, default=300,
                        help='Resolution for output (default: 300)')
    parser.add_argument('--font_size', type=int, default=10,
                        help='Font size for labels (default: 10)')
    parser.add_argument('--transparent', action='store_true',
                        help='Save with transparent background')
    
    args = parser.parse_args()

    configure_matplotlib(dpi=args.dpi, font_size=args.font_size)

    df_sites0 = load_sites(args.sites_file0)
    df_sites1 = load_sites(args.sites_file1)

    min_dist = calculate_nearest_distances(df_sites0, df_sites1, args.label0, args.label1)
    
    if len(min_dist) == 0:
        print('\nNo overlapping sites found. Cannot generate plot.')
        return
    
    print(f'\nGenerating histogram...')
    plot_distance_histogram(
        min_dist,
        args.label0,
        args.label1,
        args.dataset_label,
        args.output_file,
        xmax=args.xmax,
        figsize=tuple(args.figsize),
        fmt=args.format,
        dpi=args.dpi,
        transparent=args.transparent
    )

    print(f'\nDistance statistics:')
    print(f'  Mean: {np.mean(min_dist):.2f} nt')
    print(f'  Median: {np.median(min_dist):.2f} nt')
    print(f'  Std: {np.std(min_dist):.2f} nt')
    print(f'  Min: {np.min(min_dist)} nt')
    print(f'  Max: {np.max(min_dist)} nt')
    
    print('\nFinished')


if __name__ == '__main__':
    main()