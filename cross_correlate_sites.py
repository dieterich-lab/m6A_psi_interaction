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


def load_sites(file_path, mod_filter=None):
    """Loads a BED file, optionally filtering by the name column."""
    print(f'Loading sites from: {file_path}')
    try:
        # Standard BED files don't have headers. We need chrom, start, name, strand.
        bed_cols = ['chrom', 'chromStart', 'chromEnd', 'name', 'score', 'strand']
        df = pd.read_csv(
            file_path,
            sep='\t',
            header=None,
            usecols=[0, 1, 3, 5],
            names=[bed_cols[i] for i in [0, 1, 3, 5]],
            dtype={'chrom': str}
        )
    except (FileNotFoundError, pd.errors.EmptyDataError) as e:
        print(f"Error loading {file_path}: {e}")
        return pd.DataFrame()

    if mod_filter:
        print(f"  Filtering by modification: '{mod_filter}'")
        df = df[df['name'] == mod_filter]

    print(f'  Loaded {len(df)} sites')
    return df


def calculate_nearest_distances(df_sites1, df_sites2, scaling):
    """Calculates nearest distances from df_sites1 to df_sites2."""
    print(f'\nCalculating nearest distances...')
    min_dist = []
    
    # Group by chromosome and strand for efficiency
    df2_grouped = df_sites2.groupby(['chrom', 'strand'])

    for _, site1 in tqdm(df_sites1.iterrows(), total=len(df_sites1)):
        chrom, start1, strand = site1['chrom'], site1['chromStart'], site1['strand']
        
        try:
            # Find the corresponding group in the second dataframe
            group_key = (chrom, strand)
            sub_df_sites2 = df2_grouped.get_group(group_key)
            
            distances = (sub_df_sites2['chromStart'] - start1).values
            if strand == '-':
                distances = -distances
            if len(distances):
                # Find the distance with the minimum absolute value
                min_dist.append(distances[np.argmin(np.abs(distances))]/scaling)
        except KeyError:
            # No matching chrom/strand in df2
            continue
    
    print(f'  Found nearest neighbors for {len(min_dist)} / {len(df_sites1)} sites')
    return min_dist


def plot_distance_distribution(min_dist, args, xmax=10000):
    """Plots the cumulative distribution of nearest-neighbor distances."""
    fig_kwargs = dict(bbox_inches='tight', dpi=300)
    plt.figure(figsize=(4, 4))
    if args.xmax:
        counts, _, _ = plt.hist(min_dist, range=[-args.xmax, args.xmax], bins=2*args.xmax, color="gray")
        label = "(nt)"
    else:
        counts, _, _ = plt.hist(min_dist, color="gray")
        label = "(kb)"
    label1 = f"{args.mod_filter1} ({args.name1})" if args.mod_filter1 else f"All ({args.name1})"
    label2 = f"{args.mod_filter2} ({args.name2})" if args.mod_filter2 else f"All ({args.name2})"
    plt.xlabel(f'Nearest distance {label}\n{label2} - {label1}')
    plt.ylabel('Site count')
    plt.title(f'{int(np.sum(counts))} / {len(min_dist)} {label1} sites')
    plt.tight_layout()
    # Save the primary output (e.g., PNG)
    plt.savefig(args.output_file, format='png', **fig_kwargs)
    # Save the PDF version
    pdf_output_path = os.path.splitext(args.output_file)[0] + '.pdf'
    plt.savefig(pdf_output_path, format="pdf", **fig_kwargs)
    
    plt.close()
    print(f'Saved plot: {args.output_file}')
    print(f'Saved plot: {pdf_output_path}')


def main():
    parser = ArgumentParser(description='Cross-correlate two sets of modification sites from BED files.')
    
    parser.add_argument('--bed1', type=str, required=True,
                        help='Path to first BED file.')
    parser.add_argument('--name1', type=str, required=True,
                        help='BED file 1 description e.g. down (for dmr down sites).')
    parser.add_argument('--bed2', type=str, required=True,
                        help='Path to second BED file.')
    parser.add_argument('--name2', type=str, required=True,
                        help='BED file 2 description e.g. up (for dmr up sites).')
    parser.add_argument('--output_file', '-o', type=str, required=True,
                        help='Output file path for the plot (e.g., plot.png).')
    parser.add_argument('--mod_filter1', type=str,
                        help='Optional: filter sites in --bed1 by the "name" column (short name).')
    parser.add_argument('--mod_filter2', type=str,
                        help='Optional: filter sites in --bed2 by the "name" column  (short name).')
    parser.add_argument('--xmax', type=int, default=None,
                        help='Distance range in nt (if not given, all is shown in kb)')
    # xmax is no longer used for plotting range but could be kept for other purposes if needed.
    # For now, it's removed to avoid confusion.
    
    args = parser.parse_args()
    scaling = 1000
    label = "(kb)"
    if args.xmax:
        scaling = 1
        label = "(nt)"
        
    configure_matplotlib()

    df_sites1 = load_sites(args.bed1, mod_filter=args.mod_filter1)
    df_sites2 = load_sites(args.bed2, mod_filter=args.mod_filter2)

    if df_sites1.empty or df_sites2.empty:
        print('\nOne or both input files are empty or could not be loaded. Exiting.')
        return
    
    min_dist = calculate_nearest_distances(df_sites1, df_sites2, scaling)
    
    if not min_dist:
        print('\nNo overlapping sites found (check chromosome and strand). Cannot generate plot.')
        return
    
    print(f'\nGenerating distribution plot...')
    plot_distance_distribution(
        min_dist,
        args
    )

    print(f'\nDistance statistics:')
    print(f'  Mean: {np.mean(min_dist):.2f} {label}')
    print(f'  Median: {np.median(min_dist):.2f} {label}')
    print(f'  Std Dev: {np.std(min_dist):.2f} {label}')
    print(f'  Min: {np.min(min_dist)} {label}')
    print(f'  Max: {np.max(min_dist)} {label}')
    
    print('\nFinished.')


if __name__ == '__main__':
    main()
