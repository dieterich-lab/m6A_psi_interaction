import os
import pandas as pd
import pickle
from functools import reduce
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import pybedtools
import numpy as np
from tqdm import tqdm
from scipy import stats
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


def get_merged_df_with_delta(in_dfs, compare_conds, writer):

    merged_fields = [
        'chrom',
        'chromStart',
        'chromEnd',
        'name',
        'strand',
        'ref_motif'
    ]

    this_mod_dfs = []
    for cond in compare_conds:
        df = in_dfs[cond].copy()
        this_mod_dfs.append(
            df.rename(columns={'score': f'score_{cond}', 'frequency': f'freq_{cond}'})
        )
    out_df_merged = reduce(lambda left, right: pd.merge(left, right, on=merged_fields, how='inner'), this_mod_dfs)

    mask = ((out_df_merged.loc[:, out_df_merged.columns.str.contains('freq_')] > 0).any(axis=1))
    out_df_merged = out_df_merged[mask]

    if f'{writer}-KD' in compare_conds:
        out_df_merged[f'delta_{writer}-KD'] = out_df_merged[f'freq_{writer}-KD'] - out_df_merged['freq_CTRL']
    if f'{writer}-OE' in compare_conds:
        out_df_merged[f'delta_{writer}-OE'] = out_df_merged[f'freq_{writer}-OE'] - out_df_merged['freq_CTRL']

    return out_df_merged


def get_neighboring_sites_on_exon(in_row, dict_annots, second_mod_bedtool):
    
    if in_row['chrom'] not in dict_annots.keys():
        return []
    in_row_bedtool = pybedtools.BedTool.from_dataframe(in_row.to_frame().T)
    this_chr_annot = dict_annots[in_row['chrom']]
    in_row_exon = this_chr_annot.intersect(in_row_bedtool, u=True)
    if len(in_row_exon) > 1:
        in_row_exon = in_row_exon.sort().merge()
    neighboring_sites = second_mod_bedtool.intersect(in_row_exon)
    return neighboring_sites.to_dataframe()


def get_vec_dist_delta_from_site_df(in_df, sec_bedtool, annots):
    
    dist_delta = []
    print(f'Processing {len(in_df)} sites...')
    for _, this_row in tqdm(in_df.iterrows()):
        exon_sec_sites = get_neighboring_sites_on_exon(this_row, annots, sec_bedtool)
        if len(exon_sec_sites):
            this_vec_dist = (this_row['chromStart'] - exon_sec_sites['start']).abs().to_numpy()
            this_vec_delta = exon_sec_sites['score'].to_numpy()
            dist_delta.extend(list(zip(this_vec_dist, this_vec_delta)))
    return np.vstack(dist_delta).T


def get_binned_dist_delta(in_vec_dist, in_vec_delta, in_bin_range, in_bin_width):
    
    bin_edges = np.arange(in_bin_range[0], in_bin_range[1]+in_bin_width, in_bin_width)
    bin_centers = 0.5 * (bin_edges[1:] + bin_edges[:-1])
    num_bins = len(bin_edges) - 1
    binned_delta = []
    for bin_i in range(num_bins):
        bin_start = bin_edges[bin_i]
        bin_end = bin_edges[bin_i+1]
        mask_dist = (in_vec_dist >= bin_start) * (in_vec_dist < bin_end)
        if mask_dist.any():
            binned_delta.append(in_vec_delta[mask_dist])
        else:
            binned_delta.append(np.nan)
    return bin_centers, binned_delta


def load_chromosome_annotations(annot_dir):
    
    print(f'Loading chromosome annotations from: {annot_dir}')
    chr_annots = {}
    chromosomes = list(range(1, 23)) + ['X', 'Y', 'MT']
    
    for this_chr in tqdm(chromosomes, desc='Loading annotations'):
        annot_file = os.path.join(annot_dir, f'chr{this_chr}.exons.GRCh38.102.gtf')
        if os.path.exists(annot_file):
            chr_annots[str(this_chr)] = pybedtools.BedTool(annot_file)
        else:
            print(f'  Warning: {annot_file} not found')
    
    return chr_annots


def plot_delta_histograms(merged_df_psi, merged_df_m6a, writer, cond, output_dir, 
                          fmt='png', dpi=1200, transparent=False):
    
    cm = 1/2.54
    fig_kwargs = dict(format=fmt, bbox_inches='tight', dpi=dpi, transparent=transparent)
    
    plt.figure(figsize=(10*cm, 5*cm))
    plt.subplot(1, 2, 1)
    plt.hist(merged_df_psi[f'delta_{writer}-{cond}'], range=[-100, 100], bins=50, log=True)
    plt.xlabel('$\Delta$S($\psi$)')
    plt.ylabel('Site count')
    plt.axvline(x=0, c='r')
    plt.subplot(1, 2, 2)
    plt.hist(merged_df_m6a[f'delta_{writer}-{cond}'], range=[-100, 100], bins=50, log=True)
    plt.xlabel('$\Delta$S(m6A)')
    plt.ylabel('Site count')
    plt.axvline(x=0, c='r')
    plt.suptitle(f'{writer}-{cond} vs CTRL')
    plt.tight_layout()
    
    output_file = os.path.join(output_dir, f'hist_mods_{writer}-{cond}.{fmt}')
    plt.savefig(output_file, **fig_kwargs)
    plt.close()
    print(f'Saved histogram: {output_file}')


def plot_distance_correlation(vec_dist, vec_delta, binned_dist, binned_delta_mean,
                              writer, cond, thresh_delta, num_sites, bin_range,
                              output_dir, fmt='png', dpi=1200, transparent=False):
    
    cm = 1/2.54
    fig_kwargs = dict(format=fmt, bbox_inches='tight', dpi=dpi, transparent=transparent)
    xticks = np.int64(np.linspace(*bin_range, 5))
    
    plt.figure(figsize=(5*cm, 5*cm))
    plt.scatter(vec_dist, vec_delta, s=1, c='gray')
    plt.axhline(y=0, c='g', ls='--')
    plt.plot(binned_dist, binned_delta_mean, c='r', label='Trimmed mean')
    plt.plot(binned_dist, binned_delta_mean, 'r.')
    plt.xlim(bin_range)
    plt.ylim([-25, 25])
    plt.xticks(xticks)
    plt.xlabel('Distance from $\psi$ site on same exon (nts)')
    plt.ylabel('$\Delta$S(m6A)')
    
    if cond == 'OE':
        plt.title(f'{writer}-{cond} vs CTRL\n{num_sites} sites with $\Delta$S($\psi$)$\geq${thresh_delta}')
    elif cond == 'KD':
        plt.title(f'{writer}-{cond} vs CTRL\n{num_sites} sites with $\Delta$S($\psi$)<-{thresh_delta}')
    
    plt.legend()
    output_file = os.path.join(output_dir, f'dist_corr_m6a_from_psi_{writer}-{cond}_thresh{thresh_delta}.{fmt}')
    plt.savefig(output_file, **fig_kwargs)
    plt.close()
    print(f'Saved correlation plot: {output_file}')


def main():
    parser = ArgumentParser(description='Analyze correlation between psi and m6A modifications at different distances')
    
    parser.add_argument('--pickle_file', type=str, required=True,
                        help='Pickle file with filtered modification dataframes')
    parser.add_argument('--annot_dir', type=str, required=True,
                        help='Directory containing exon annotation GTF files')
    parser.add_argument('--output_dir', '-o', type=str, required=True,
                        help='Output directory for plots')
    parser.add_argument('--writer', type=str, default='TRUB1',
                        help='Writer enzyme name (default: TRUB1)')
    parser.add_argument('--condition', type=str, required=True, choices=['KD', 'OE'],
                        help='Condition to analyze (KD or OE)')
    parser.add_argument('--thresh_delta', type=float, default=5.0,
                        help='Threshold for delta S filtering (default: 5.0)')
    parser.add_argument('--bin_range', type=int, nargs=2, default=[0, 2000],
                        help='Distance bin range (min max) (default: 0 2000)')
    parser.add_argument('--bin_width', type=int, default=200,
                        help='Width of distance bins (default: 200)')
    parser.add_argument('--format', type=str, default='png', choices=['png', 'pdf', 'svg'],
                        help='Output format (default: png)')
    parser.add_argument('--dpi', type=int, default=1200,
                        help='Resolution for output (default: 1200)')
    parser.add_argument('--font_size', type=int, default=8,
                        help='Font size for labels (default: 8)')
    parser.add_argument('--transparent', action='store_true',
                        help='Save with transparent background')
    
    args = parser.parse_args()
    
    configure_matplotlib(dpi=args.dpi, font_size=args.font_size)
    
    os.makedirs(args.output_dir, exist_ok=True)
    
    print(f'Loading data from: {args.pickle_file}')
    with open(args.pickle_file, 'rb') as pkl_in:
        dfs_mod_cond = pickle.load(pkl_in)
    
    chr_annots = load_chromosome_annotations(args.annot_dir)
    
    print('\nMerging psi dataframes...')
    merged_df_psi = get_merged_df_with_delta(
        dfs_mod_cond['17802'], 
        ['CTRL', f'{args.writer}-{args.condition}'],
        args.writer
    )
    print(f'  Merged psi sites: {len(merged_df_psi)}')
    
    print('\nMerging m6A dataframes...')
    merged_df_m6a = get_merged_df_with_delta(
        dfs_mod_cond['a'], 
        ['CTRL', f'{args.writer}-{args.condition}'],
        args.writer
    )
    print(f'  Merged m6A sites: {len(merged_df_m6a)}')
    
    plot_delta_histograms(
        merged_df_psi, merged_df_m6a, args.writer, args.condition, 
        args.output_dir, args.format, args.dpi, args.transparent
    )
    
    merged_df_m6a_bed = merged_df_m6a.copy()
    merged_df_m6a_bed.iloc[:, 4] = merged_df_m6a_bed[f'delta_{args.writer}-{args.condition}']
    m6a_bedtool = pybedtools.BedTool.from_dataframe(merged_df_m6a_bed.iloc[:, :6])
    
    if args.condition == 'KD':
        df_psi_filtered = merged_df_psi[
            merged_df_psi[f'delta_{args.writer}-{args.condition}'] < -args.thresh_delta
        ]
        print(f'\nPsi sites with delta < -{args.thresh_delta}: {len(df_psi_filtered)}')
    elif args.condition == 'OE':
        df_psi_filtered = merged_df_psi[
            merged_df_psi[f'delta_{args.writer}-{args.condition}'] >= args.thresh_delta
        ]
        print(f'\nPsi sites with delta >= {args.thresh_delta}: {len(df_psi_filtered)}')
    
    if len(df_psi_filtered) == 0:
        print('No sites passed the threshold. Exiting.')
        return
    
    vec_dist, vec_delta = get_vec_dist_delta_from_site_df(df_psi_filtered, m6a_bedtool, chr_annots)
    
    if len(vec_dist) == 0:
        print('No neighboring sites found. Exiting.')
        return
    
    print(f'Found {len(vec_dist)} distance-delta pairs')
    
    binned_dist, binned_delta = get_binned_dist_delta(
        vec_dist, vec_delta, args.bin_range, args.bin_width
    )
    
    binned_delta_lower = np.array([np.quantile(this_bin, 0.01) if not isinstance(this_bin, float) 
                                   else np.nan for this_bin in binned_delta])
    binned_delta_upper = np.array([np.quantile(this_bin, 0.99) if not isinstance(this_bin, float) 
                                   else np.nan for this_bin in binned_delta])
    binned_delta_mean = 0.5 * (binned_delta_upper + binned_delta_lower)
    
    plot_distance_correlation(
        vec_dist, vec_delta, binned_dist, binned_delta_mean,
        args.writer, args.condition, args.thresh_delta, len(df_psi_filtered),
        args.bin_range, args.output_dir, args.format, args.dpi, args.transparent
    )
    
    print('\nFinished')


if __name__ == '__main__':
    main()