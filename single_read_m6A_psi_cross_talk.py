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


dict_mod_display = {
    'm6A': 'm^6A',
    'psi': '\psi'
}


mod_tags = {
    'm6A': ('A', 0, 'a'),
    'psi': ('T', 0, 17802)
}
# mod_tags = {
#     'm6A': ('A', 1, 'a'),
#     'psi': ('T', 1, 17802)
# }


mod_motifs = {
    'm6A': [
            'GGACT', 'GGACA', 'GAACT', 'AGACT', 'GGACC', 'TGACT',
            'AAACT', 'GAACA', 'AGACA', 'AGACC', 'GAACC', 'TGACA',
            'TAACT', 'AAACA', 'TGACC', 'TAACA', 'AAACC', 'TAACC'
    ],
    'psi': ['GTTCA', 'GTTCC', 'GTTCG', 'GTTCT'] + ['TGTAG'] +
        [f'{pos1}{pos2}T{pos4}{pos5}'
            for pos1 in ['A', 'C', 'G', 'T']
            for pos2 in ['A', 'G']
            for pos4 in ['A', 'G']
            for pos5 in ['A', 'C', 'G', 'T']
            ]
}
        
rc_mod_motifs = {
    "m6A" : {str(Seq(m).reverse_complement()) for m in mod_motifs["m6A"]},
    "psi" : {str(Seq(m).reverse_complement()) for m in mod_motifs["psi"]},
}

            
def get_mean_logit(in_probs, num_top_locs=5):
    if len(in_probs) == 0:
        return np.nan
    top_probs = np.sort(in_probs)[-num_top_locs:]
    rescaled_probs = np.clip(np.array(top_probs) / 255.0, a_max=0.999, a_min=0.001)
    logits = np.log2(rescaled_probs / (1-rescaled_probs))
    return np.mean(logits)


def get_mean_logit_mod_level(in_read):
    mod_mean_logit = {}
    for this_mod, this_tag in mod_tags.items():
        this_mod_probs = [this_tup[1] for this_tup in in_read.modified_bases.get(this_tag, [])]
        mod_mean_logit[this_mod] = get_mean_logit(this_mod_probs)
    return mod_mean_logit


def get_mod_mean_occupancy(in_read, min_locs=10):
    mod_mean_occupancy = {}
    n_sites = {}
    for this_mod, this_tag in mod_tags.items():
        
        # if reverse 1, else 0
        in_tag = (this_tag[0], in_read.is_reverse, this_tag[2])
        
        tups = in_read.modified_bases.get(in_tag, [])
        this_mod_locs = [this_tup[0] for this_tup in tups]
        if len(this_mod_locs) < min_locs:
            mod_mean_occupancy[this_mod] = np.nan
            n_sites[this_mod] = np.nan
            continue
        loc_motifs = [in_read.query_sequence[this_loc - 2:this_loc + 3] for this_loc in this_mod_locs]
        if in_read.is_reverse:
            motif_included = np.array([this_motif in rc_mod_motifs[this_mod] for this_motif in loc_motifs])
        else:
            motif_included = np.array([this_motif in mod_motifs[this_mod] for this_motif in loc_motifs])
        this_mod_probs = np.array([this_tup[1] for this_tup in tups]) / 255.0
        filtered_mod_probs = this_mod_probs[motif_included]
        if len(filtered_mod_probs) >= min_locs:
            mod_mean_occupancy[this_mod] = np.mean(filtered_mod_probs >= 0.5)
            n_sites[this_mod] = np.sum(filtered_mod_probs >= 0.5)
        else:
            mod_mean_occupancy[this_mod] = np.nan
            n_sites[this_mod] = np.nan
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
    
    args = parser.parse_args()
    
    os.makedirs(args.img_out, exist_ok=True)
    
    print(f'Processing {args.bam_file}...')

    single_read_mean_occupancy = []
    single_read_n_sites = []
    with pysam.AlignmentFile(args.bam_file, 'rb', check_sq=False) as bam:
        # for this_read in tqdm(bam.fetch(until_eof=True)):
        for this_read in bam.fetch(until_eof=True):
            if this_read.modified_bases is not None:
                # single_read_mean_logit.append(get_mean_logit_mod_level(this_read))
                # single_read_mean_occupancy.append(get_mod_mean_occupancy(this_read))
                occupancy, sites = get_mod_mean_occupancy(this_read)
                single_read_mean_occupancy.append(occupancy)
                single_read_n_sites.append(sites)

    vec_m6A, vec_psi = np.vstack([
        (this_read_mean_occupancy['m6A'],  this_read_mean_occupancy['psi'])
        for this_read_mean_occupancy in single_read_mean_occupancy
        if ~np.isnan(this_read_mean_occupancy['m6A']) and ~np.isnan(this_read_mean_occupancy['psi'])
    ]).T
    vec_m6A_sites, vec_psi_sites = np.vstack([
        (this_read_n_sites['m6A'],  this_read_n_sites['psi'])
        for this_read_n_sites in single_read_n_sites
        if ~np.isnan(this_read_n_sites['m6A']) and ~np.isnan(this_read_n_sites['psi'])
    ]).T
    # num_valid_reads = len(vec_m6A)
    
    print(f'... Done! Writing to {args.img_out}...')

    out_file = os.path.join(args.img_out, f'single_read_occ_m6A_psi_{args.ds}.npz')
    np.savez(out_file, vec_m6A=vec_m6A, vec_psi=vec_psi)
    
    out_file = os.path.join(args.img_out, f'single_read_sites_m6A_psi_{args.ds}.npz')
    np.savez(out_file, vec_m6A=vec_m6A_sites, vec_psi=vec_psi_sites)
    
    corr, pval = pearsonr(vec_psi, vec_m6A)
    corr = "{:.2}".format(corr)
    pval = "{:.2E}".format(pval)

    top_reads = [10**i for i in range(2, 7)]
    top_labels = [f"Top $10^{i}$ reads" for i in range(2, 7)]
    cmap = plt.cm.Greys
    colors = cmap(np.linspace(0.5, 1, len(top_reads)))

    binned_psi = []
    binned_m6A = []
    bin_edges = np.round(np.linspace(0, 1, 6), 1)
    for bin_i in range(len(bin_edges)-1):
        bin_start = bin_edges[bin_i]
        bin_end = bin_edges[bin_i+1]
        if bin_end == 1.0:
            bin_end += 0.001
        mask_m6A = (vec_m6A >= bin_start) * (vec_m6A < bin_end)
        binned_psi.append(vec_psi[mask_m6A])
        mask_psi = (vec_psi >= bin_start) * (vec_psi < bin_end)
        binned_m6A.append(vec_m6A[mask_psi])

    flierprops = dict(marker='o', markerfacecolor='none', markersize=2, markeredgecolor='gray',
                    alpha=0.5, rasterized=True)

    xy_ticks = np.int32(bin_edges * 100)

    outfile_name = os.path.join(args.img_out, f'boxplot_mean_occupancy_per_read_{args.ds}.{FMT}')

    plt.figure(figsize=(8*cm, 7*cm))
    plt.subplot(2, 2, 1)
    plt.boxplot(binned_psi, flierprops=flierprops)
    plt.ylim([-0.01, 1.05])
    plt.xticks(np.arange(len(bin_edges)) + 0.5, xy_ticks)
    plt.yticks(bin_edges, xy_ticks)
    plt.xlabel(f"occ(${dict_mod_display['m6A']}$) per read")
    plt.ylabel(f"occ(${dict_mod_display['psi']}$) per read")
    plt.text(1, .95, rf"$\rho$={corr}, p-value={pval}", fontsize="xx-small")
    plt.subplot(2, 2, 2)
    plt.boxplot(binned_m6A, flierprops=flierprops)
    plt.ylim([-0.01, 1.05])
    plt.xticks(np.arange(len(bin_edges)) + 0.5, xy_ticks)
    plt.yticks(bin_edges, xy_ticks)
    plt.xlabel(f"occ(${dict_mod_display['psi']}$) per read")
    plt.ylabel(f"occ(${dict_mod_display['m6A']}$) per read")
    plt.subplot(2, 2, 3)
    for top_read, top_label, color in zip(top_reads, top_labels, colors):
        top_psi = [np.median(np.sort(this_bin)[-top_read:]) for this_bin in binned_psi]
        plt.plot(np.arange(1, len(bin_edges)), top_psi, color=color, ls='-', label=top_label)
        plt.plot(np.arange(1, len(bin_edges)), top_psi, 'o', color=color, markersize=2)
    # plt.legend(loc='upper right', title="Median", title_fontsize="medium",fontsize="x-small", frameon=False)
    plt.legend(loc='upper right',fontsize="x-small", frameon=False)
    plt.xticks(np.arange(len(bin_edges)) + 0.5, xy_ticks)
    plt.ylim([-0.01, 1.05])
    plt.yticks(bin_edges, xy_ticks)
    plt.xlabel(f"occ(${dict_mod_display['m6A']}$) per read")
    plt.ylabel(f"Mdn. occ(${dict_mod_display['psi']}$) per read")
    plt.subplot(2, 2, 4)
    for top_read, top_label, color in zip(top_reads, top_labels, colors):
        top_m6A = [np.median(np.sort(this_bin)[-top_read:]) for this_bin in binned_m6A]
        plt.plot(np.arange(1, len(bin_edges)), top_m6A, color=color, ls='-', label=top_label)
        plt.plot(np.arange(1, len(bin_edges)), top_m6A, 'o', color=color, markersize=2)
    plt.xticks(np.arange(len(bin_edges)) + 0.5, xy_ticks)
    plt.ylim([-0.01, 1.05])
    plt.yticks(bin_edges, xy_ticks)
    plt.xlabel(f"occ(${dict_mod_display['psi']}$) per read")
    plt.ylabel(f"Mdn. occ(${dict_mod_display['m6A']}$) per read")
    plt.tight_layout()
    plt.savefig(outfile_name, format=FMT, **fig_kwargs)
    plt.savefig(f"{outfile_name}.pdf", format="pdf", **fig_kwargs)
    plt.savefig(f"{outfile_name}.svg", format="svg", **fig_kwargs)


if __name__ == '__main__':
    main()
    print('... finished!')