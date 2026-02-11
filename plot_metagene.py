import os
import pandas as pd
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
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


def get_longest_isoform(in_df):
    # Ensure consistent types for merging/grouping
    in_df['chr'] = in_df['chr'].astype(str)
    
    # Calculate transcript length
    in_df['tx_len'] = in_df['utr5_size'] + in_df['cds_size'] + in_df['utr3_size']
    
    # Sort by tx_len descending
    in_df_sorted = in_df.sort_values('tx_len', ascending=False)
    
    # Drop duplicates by chr/coord, keeping first (longest)
    out_df = in_df_sorted.drop_duplicates(subset=['chr', 'coord'])
    
    # Restore sorting by coordinate
    out_df = out_df.sort_values(['chr', 'coord'])
    return out_df


def plot_metagene(df, output_path, plot_name='metagene_plot', fmt='png', 
                  bins=30, figsize=(5, 5), dpi=1200, transparent=False):
    
    cm = 1/2.54  # centimeters in inches
    fig_kwargs = dict(format=fmt, bbox_inches='tight', dpi=dpi, transparent=transparent)
    
    plt.figure(figsize=(figsize[0]*cm, figsize[1]*cm))
    plt.hist(df['rel_location'], range=[0, 3], bins=bins)
    plt.xlabel('Gene region')
    plt.ylabel('Site count')
    plt.xticks([])
    plt.xlim([0, 3])
    plt.axvline(x=1, c='gray', ls='--')
    plt.axvline(x=2, c='gray', ls='--')
    
    # Add labels for gene regions
    ymin, ymax = plt.ylim()
    plt.text(0.5, ymax * 0.95, "5'UTR", ha='center', va='top')
    plt.text(1.5, ymax * 0.95, "CDS", ha='center', va='top')
    plt.text(2.5, ymax * 0.95, "3'UTR", ha='center', va='top')

    os.makedirs(output_path, exist_ok=True)
    output_file = os.path.join(output_path, f'{plot_name}.{fmt}')
    plt.savefig(output_file, **fig_kwargs)
    plt.close()
    print(f'Plot saved to: {output_file}')


def main():
    parser = ArgumentParser(description='Plot metagene distribution of modification sites')
    
    parser.add_argument('--input', '-i', type=str, required=True,
                        help='Input file with metagene distribution data (TSV format)')
    parser.add_argument('--output_dir', '-o', type=str, required=True,
                        help='Output directory for plots')
    parser.add_argument('--plot_name', type=str, default='metagene_plot',
                        help='Base name for output plot (default: metagene_plot)')
    parser.add_argument('--format', type=str, default='png', choices=['png', 'pdf', 'svg'],
                        help='Output format (default: png)')
    parser.add_argument('--bins', type=int, default=30,
                        help='Number of histogram bins (default: 30)')
    parser.add_argument('--figsize', type=float, nargs=2, default=[5.0, 5.0],
                        help='Figure size in cm (width height) (default: 5 5)')
    parser.add_argument('--dpi', type=int, default=1200,
                        help='Resolution for output (default: 1200)')
    parser.add_argument('--font_size', type=int, default=8,
                        help='Font size for labels (default: 8)')
    parser.add_argument('--transparent', action='store_true',
                        help='Save with transparent background')
    parser.add_argument('--use_longest_isoform', action='store_true',
                        help='Filter to use only the longest isoform per location')
    
    args = parser.parse_args()
    
    configure_matplotlib(dpi=args.dpi, font_size=args.font_size)
    
    print(f'Reading data from: {args.input}')
    df_metagene = pd.read_csv(args.input, sep='\t')
    
    if args.use_longest_isoform:
        print('Filtering to longest isoforms...')
        df_metagene = get_longest_isoform(df_metagene)
    
    plot_metagene(
        df_metagene,
        args.output_dir,
        plot_name=args.plot_name,
        fmt=args.format,
        bins=args.bins,
        figsize=args.figsize,
        dpi=args.dpi,
        transparent=args.transparent
    )
    
    print('Finished')


if __name__ == '__main__':
    main()