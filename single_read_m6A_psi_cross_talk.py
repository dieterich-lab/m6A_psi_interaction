import os
import pandas as pd
import pysam
import numpy as np
from Bio.Seq import Seq
from scipy.stats import pearsonr
from argparse import ArgumentParser
import matplotlib as mpl
#######################################################################
cm = 1/2.54  # centimeters in inches
gr = 1.618
dpi = 1200
mpl.rcParams['figure.dpi'] = dpi
mpl.rcParams['savefig.dpi'] = dpi
mpl.rcParams['font.size'] = 6
mpl.rcParams['legend.fontsize'] = 4
mpl.rcParams['xtick.labelsize'] = 5
mpl.rcParams['ytick.labelsize'] = 5
mpl.rcParams['xtick.major.size'] = 1.5
mpl.rcParams['ytick.major.size'] = 1.5
mpl.rcParams['lines.linewidth'] = 0.5
mpl.rcParams['font.family'] = 'DejaVu Sans'
FMT = 'png'
fig_kwargs = dict(bbox_inches='tight', dpi=dpi)
######################################################################
mpl.use('Agg')
import matplotlib.pyplot as plt
from tqdm import tqdm


dict_display_labels = {
    'a': 'm^6A',
    '69426': 'Am',
    '17596': 'I',
    '17802': '\psi',
    '19227': 'Um',
    '19228': 'Cm',
    'm': 'm^5C',
    '19229': 'Gm',
}

MOD_INFO = {
    # mod_name: (canonical_base, mod_code, display_name, motifs)
    'm6A': ('A', 'a', 'm^6A', [
        'GGACT', 'GGACA', 'GAACT', 'AGACT', 'GGACC', 'TGACT',
        'AAACT', 'GAACA', 'AGACA', 'AGACC', 'GAACC', 'TGACA',
        'TAACT', 'AAACA', 'TGACC', 'TAACA', 'AAACC', 'TAACC'
    ]),
    'psi': ('T', 17802, '\psi', ['GTTCA', 'GTTCC', 'GTTCG', 'GTTCT'] + ['TGTAG'] +
        [f'{pos1}{pos2}T{pos4}{pos5}'
            for pos1 in ['A', 'C', 'G', 'T']
            for pos2 in ['A', 'G']
            for pos4 in ['A', 'G']
            for pos5 in ['A', 'C', 'G', 'T']
        ]
    ),
    'Am': ('A', 69426, 'Am', []),
    'I': ('A', 17596, 'I', []),
    'Um': ('T', 19227, 'Um', []),
    'Cm': ('C', 19228, 'Cm', []),
    'm5C': ('C', 'm', 'm^5C', []),
    'Gm': ('G', 19229, 'Gm', []),
}

# Pre-compute reverse-complemented motifs
for mod_name, info in MOD_INFO.items():
    motifs = info[3]
    rc_motifs = {str(Seq(m).reverse_complement()) for m in motifs}
    MOD_INFO[mod_name] = info + (rc_motifs,)


            
def get_mean_logit(in_probs, num_top_locs=5):
    if len(in_probs) == 0:
        return np.nan
    top_probs = np.sort(in_probs)[-num_top_locs:]
    rescaled_probs = np.clip(np.array(top_probs) / 255.0, a_max=0.999, a_min=0.001)
    logits = np.log2(rescaled_probs / (1-rescaled_probs))
    return np.mean(logits)


def get_mean_logit_mod_level(in_read, mod_names):
    mod_mean_logit = {}
    for mod_name in mod_names:
        base, _, mod_code, _, _ = MOD_INFO[mod_name]
        # This logic for getting the tag might need adjustment based on strand-specificity
        tag = (base, 0, mod_code) 
        this_mod_probs = [this_tup[1] for this_tup in in_read.modified_bases.get(tag, [])]
        mod_mean_logit[mod_name] = get_mean_logit(this_mod_probs)
    return mod_mean_logit


def get_mod_mean_occupancy(in_read, mod_names, min_locs=10, filter_by_motifs=False, strand_specific=True):
    mod_mean_occupancy = {}
    n_sites = {}
    for mod_name in mod_names:
        base, mod_code, _, motifs, rc_motifs = MOD_INFO[mod_name]
        
        strand_idx = 1 if strand_specific and in_read.is_reverse else 0
        in_tag = (base, strand_idx, mod_code)
        
        tups = in_read.modified_bases.get(in_tag, [])
        this_mod_locs = [this_tup[0] for this_tup in tups]

        if len(this_mod_locs) < min_locs:
            mod_mean_occupancy[mod_name] = np.nan
            n_sites[mod_name] = np.nan
            continue
            
        this_mod_probs = np.array([this_tup[1] for this_tup in tups]) / 255.0
        
        if filter_by_motifs and motifs:
            loc_motifs = [in_read.query_sequence[this_loc - 2:this_loc + 3] for this_loc in this_mod_locs]
            current_motifs = rc_motifs if in_read.is_reverse else motifs
            motif_included = np.array([this_motif in current_motifs for this_motif in loc_motifs])
            filtered_mod_probs = this_mod_probs[motif_included]
        else:
            filtered_mod_probs = this_mod_probs

        if len(filtered_mod_probs) >= min_locs:
            mod_mean_occupancy[mod_name] = np.mean(filtered_mod_probs >= 0.5)
            n_sites[mod_name] = np.sum(filtered_mod_probs >= 0.5)
        else:
            mod_mean_occupancy[mod_name] = np.nan
            n_sites[mod_name] = np.nan
    return mod_mean_occupancy, n_sites


def main():
    home = os.environ['HOME']
    parser = ArgumentParser()
    parser.add_argument('--img_out', type=str, default=home,
                        help='output directory')
    parser.add_argument('--bam_file', type=str, required=True,
                        help='BAM file')
    parser.add_argument('--ds', type=str, required=True,
                        help='Name')
    parser.add_argument('--mod_names', type=str, nargs=2, default=['m6A', 'psi'],
                        help='Two modification names to compare from: ' + ', '.join(MOD_INFO.keys()))
    parser.add_argument('--filter_by_motifs', action='store_true',
                        help='Filter modification sites by known motifs.')
    parser.add_argument('--no_strand_specific', dest='strand_specific', action='store_false',
                        help='Do not treat strands separately for modification detection.')
    
    args = parser.parse_args()
    
    os.makedirs(args.img_out, exist_ok=True)
    
    mod1_name, mod2_name = args.mod_names
    
    print(f'Processing {args.bam_file}...')
    print(f'Comparing {mod1_name} and {mod2_name}')
    if args.filter_by_motifs:
        print('Filtering by motifs.')
    if not args.strand_specific:
        print('Not using strand-specific detection.')

    single_read_mean_occupancy = []
    single_read_n_sites = []
    with pysam.AlignmentFile(args.bam_file, 'rb', check_sq=False) as bam:
        for this_read in tqdm(bam.fetch(until_eof=True)):
            if this_read.modified_bases is not None:
                occupancy, sites = get_mod_mean_occupancy(
                    this_read, 
                    args.mod_names,
                    filter_by_motifs=args.filter_by_motifs, 
                    strand_specific=args.strand_specific
                )
                single_read_mean_occupancy.append(occupancy)
                single_read_n_sites.append(sites)

    valid_occupancies = [
        (occ[mod1_name], occ[mod2_name])
        for occ in single_read_mean_occupancy
        if not np.isnan(occ[mod1_name]) and not np.isnan(occ[mod2_name])
    ]
    if not valid_occupancies:
        print("No reads found with both modifications. Exiting.")
        return
        
    vec_mod1, vec_mod2 = np.vstack(valid_occupancies).T

    valid_sites = [
        (sites[mod1_name], sites[mod2_name])
        for sites in single_read_n_sites
        if not np.isnan(sites[mod1_name]) and not np.isnan(sites[mod2_name])
    ]
    vec_mod1_sites, vec_mod2_sites = np.vstack(valid_sites).T
    
    print(f'... Done! Writing to {args.img_out}...')

    out_file = os.path.join(args.img_out, f'single_read_occ_{mod1_name}_{mod2_name}_{args.ds}.npz')
    np.savez(out_file, **{f'vec_{mod1_name}': vec_mod1, f'vec_{mod2_name}': vec_mod2})
    
    out_file = os.path.join(args.img_out, f'single_read_sites_{mod1_name}_{mod2_name}_{args.ds}.npz')
    np.savez(out_file, **{f'vec_{mod1_name}_sites': vec_mod1_sites, f'vec_{mod2_name}_sites': vec_mod2_sites})
    
    corr, pval = pearsonr(vec_mod2, vec_mod1)
    corr = "{:.2}".format(corr)
    pval = "{:.2E}".format(pval)

    top_reads = [10**i for i in range(2, 7)]
    top_labels = [f"Top $10^{i}$ reads" for i in range(2, 7)]
    cmap = plt.cm.Greys
    colors = cmap(np.linspace(0.5, 1, len(top_reads)))

    binned_mod2 = []
    binned_mod1 = []
    bin_edges = np.round(np.linspace(0, 1, 6), 1)
    for bin_i in range(len(bin_edges)-1):
        bin_start = bin_edges[bin_i]
        bin_end = bin_edges[bin_i+1]
        if bin_end == 1.0:
            bin_end += 0.001
        mask_mod1 = (vec_mod1 >= bin_start) * (vec_mod1 < bin_end)
        binned_mod2.append(vec_mod2[mask_mod1])
        mask_mod2 = (vec_mod2 >= bin_start) * (vec_mod2 < bin_end)
        binned_mod1.append(vec_mod1[mask_mod2])

    flierprops = dict(marker='o', markerfacecolor='none', markersize=2, markeredgecolor='gray',
                    alpha=0.5, rasterized=True)

    xy_ticks = np.int32(bin_edges * 100)

    outfile_name = os.path.join(args.img_out, f'boxplot_mean_occupancy_per_read_{args.ds}.{FMT}')
    
    mod1_display = MOD_INFO[mod1_name][2]
    mod2_display = MOD_INFO[mod2_name][2]

    plt.figure(figsize=(8*cm, 7*cm))
    plt.subplot(2, 2, 1)
    bp1 = plt.boxplot(binned_mod2, flierprops=flierprops)
    # Add counts to boxplot
    counts = [len(bin_data) for bin_data in binned_mod2]
    for i, count in enumerate(counts):
        y = np.median(bp1['medians'][i].get_ydata())
        plt.text(i + 1, 1.06, f'n={count}', ha='center', va='bottom', fontsize='x-small')
    plt.ylim([-0.01, 1.15])
    plt.xticks(np.arange(len(bin_edges)) + 0.5, xy_ticks)
    plt.yticks(bin_edges, xy_ticks)
    plt.xlabel(f"occ(${mod1_display}$) per read")
    plt.ylabel(f"occ(${mod2_display}$) per read")
    plt.text(1, .95, rf"$\rho$={corr}, p-value={pval}", fontsize="xx-small")
    plt.subplot(2, 2, 2)
    bp2 = plt.boxplot(binned_mod1, flierprops=flierprops)
    # Add counts to boxplot
    counts = [len(bin_data) for bin_data in binned_mod1]
    for i, count in enumerate(counts):
        y = np.median(bp2['medians'][i].get_ydata())
        plt.text(i + 1, 1.06, f'n={count}', ha='center', va='bottom', fontsize='x-small')
    plt.ylim([-0.01, 1.15])
    plt.xticks(np.arange(len(bin_edges)) + 0.5, xy_ticks)
    plt.yticks(bin_edges, xy_ticks)
    plt.xlabel(f"occ(${mod2_display}$) per read")
    plt.ylabel(f"occ(${mod1_display}$) per read")
    plt.subplot(2, 2, 3)
    for top_read, top_label, color in zip(top_reads, top_labels, colors):
        top_mod2 = [np.median(np.sort(this_bin)[-top_read:]) if len(this_bin) > 0 else np.nan for this_bin in binned_mod2]
        plt.plot(np.arange(1, len(bin_edges)), top_mod2, color=color, ls='-', label=top_label)
        plt.plot(np.arange(1, len(bin_edges)), top_mod2, 'o', color=color, markersize=2)
    # plt.legend(loc='upper right', title="Median", title_fontsize="medium",fontsize="x-small", frameon=False)
    plt.legend(loc='upper right',fontsize="x-small", frameon=False)
    plt.xticks(np.arange(len(bin_edges)) + 0.5, xy_ticks)
    plt.ylim([-0.01, 1.05])
    plt.yticks(bin_edges, xy_ticks)
    plt.xlabel(f"occ(${mod1_display}$) per read")
    plt.ylabel(f"Mdn. occ(${mod2_display}$) per read")
    plt.subplot(2, 2, 4)
    for top_read, top_label, color in zip(top_reads, top_labels, colors):
        top_mod1 = [np.median(np.sort(this_bin)[-top_read:]) for this_bin in binned_mod1]
        plt.plot(np.arange(1, len(bin_edges)), top_mod1, color=color, ls='-', label=top_label)
        plt.plot(np.arange(1, len(bin_edges)), top_mod1, 'o', color=color, markersize=2)
    plt.xticks(np.arange(len(bin_edges)) + 0.5, xy_ticks)
    plt.ylim([-0.01, 1.05])
    plt.yticks(bin_edges, xy_ticks)
    plt.xlabel(f"occ(${mod2_display}$) per read")
    plt.ylabel(f"Mdn. occ(${mod1_display}$) per read")
    plt.tight_layout()
    plt.savefig(outfile_name, format=FMT, **fig_kwargs)
    plt.savefig(f"{outfile_name}.pdf", format="pdf", **fig_kwargs)
    plt.savefig(f"{outfile_name}.svg", format="svg", **fig_kwargs)


if __name__ == '__main__':
    main()
    print('... finished!')