import pandas as pd
import os
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from tqdm import tqdm
from argparse import ArgumentParser


dict_mod_display = {
    'a': 'm^6A',
    '17802': '\psi'
}


def get_df_gtf(in_args):
    print(f'Loading gtf {in_args.gtf_exon}')
    gtf_fields = [
        'chrom',
        'source',
        'feature',
        'start',
        'end',
        'score',
        'strand',
        'frame',
        'attribute'
    ]

    out_df_gtf = pd.read_csv(in_args.gtf_exon, sep='\t', names=gtf_fields, dtype={'chrom': str})
    out_df_gtf = out_df_gtf[out_df_gtf['source'] == 'ensembl_havana']
    return out_df_gtf


def get_df_bed_mod(in_args, mod_code):
    print(f'Loading bedmethyl(s) {in_args.bedmethyl}')

    bedmethyl_fields = [
        'chrom',
        'chromStart',
        'chromEnd',
        'name',
        'score',
        'strand',
        'thickStart',
        'thickEnd',
        'itemRgb',
        'coverage',
        'frequency'
    ]
    sel_bedmethyl_fields = [bedmethyl_fields[i] for i in list(range(6)) + [9,10]]
    out_df_bed_mod_strand_list = []
    for bedmethyl in in_args.bedmethyl:
        out_df_bed = pd.read_csv(bedmethyl, sep='\t', usecols=list(range(6)) + [9,10],
                                names=sel_bedmethyl_fields, dtype={'chrom': str})
        out_df_bed_mod = out_df_bed[(out_df_bed['name'] == mod_code)&(out_df_bed['coverage']>=in_args.thresh_count)]
        out_df_bed_mod_strand = {
            '+': out_df_bed_mod[out_df_bed_mod['strand'] == '+'],
            '-': out_df_bed_mod[out_df_bed_mod['strand'] == '-']
        }
        out_df_bed_mod_strand_list.append(out_df_bed_mod_strand)
    return out_df_bed_mod_strand_list


def get_exon_junction_freq(in_df_gtf, in_df_bed_list, in_args):
    print(f'Parsing sites in {len(in_df_gtf)} exons...')

    out_junction_freq = []
    out_non_junction_freq = []
    # for _, this_row in tqdm(in_df_gtf.iterrows()):
    for _, this_row in in_df_gtf.iterrows():
        this_chrom_start = this_row['start'] - 1
        this_chrom_end = this_row['end']
        this_strand = str(this_row['strand'])
        
        for in_df_bed in in_df_bed_list:
            this_strand_bed = in_df_bed[this_strand]
            sub_df_bed = this_strand_bed[
                (this_strand_bed['chromStart'] >= this_chrom_start)
                * (this_strand_bed['chromEnd'] <= this_chrom_end)
                ]

            if len(sub_df_bed):
                mask = (sub_df_bed['chromStart'] <= (this_chrom_start + in_args.exon_junction_margin)) \
                    + (sub_df_bed['chromStart'] >= (this_chrom_end - in_args.exon_junction_margin))
                out_junction_freq.extend(sub_df_bed[mask]['frequency'].values)
                out_non_junction_freq.extend(sub_df_bed[~mask]['frequency'].values)

    return out_junction_freq, out_non_junction_freq


def make_plot_by_2(in_junction_freq, in_non_junction_freq, in_args):

    bin_range = [0, 100]
    num_bins = 10
    
    # len(in_args.mod_code) == 2!

    plt.figure(figsize=(8, 4))
    for idx, mod_code in enumerate(in_args.mod_code):
        plt.subplot(1,2,idx+1)
        plt.hist(in_junction_freq[mod_code], density=True, log=True, range=bin_range, bins=num_bins,
                alpha=0.5, color='r', histtype='step', label=f'$\leq${in_args.exon_junction_margin} nt of exon edge')
        plt.hist(in_non_junction_freq[mod_code], density=True, log=True, range=bin_range, bins=num_bins,
                alpha=0.5, color='b', histtype='step', label=f'>{in_args.exon_junction_margin} nt')
        if idx == 0:
            plt.ylabel('Density')
        plt.xlabel(f'Stoichiometry (${dict_mod_display[mod_code]}$)')
        plt.xlim(bin_range)
    plt.legend()
    out_file = os.path.join(in_args.img_out, f'exon_junction_mod_freq_{in_args.ds}.png')
    plt.savefig(out_file, bbox_inches='tight')
    out_file = os.path.join(in_args.img_out, f'exon_junction_mod_freq_{in_args.ds}.pdf')
    plt.savefig(out_file, format="pdf", bbox_inches='tight')
    out_file = os.path.join(in_args.img_out, f'exon_junction_mod_freq_{in_args.ds}.svg')
    plt.savefig(out_file, format="svg", bbox_inches='tight')

def main():
    home = os.environ['HOME']
    parser = ArgumentParser()
    parser.add_argument('--gtf_exon', type=str, required=True,
                        help='annotation file for exons')
    parser.add_argument('--bedmethyl', type=str, required=True, nargs='+',
                        help='bedmethyl output from modkit pileup')
    parser.add_argument('--mod_code', type=str, required=True, nargs='+',
                        help='modification code in bedmethyl, eg, "a" for m6A, "17802" for pseudouridine')
    parser.add_argument('--exon_junction_margin', type=int, default=20,
                        help='cut-off for classification of sites on exon junction')
    parser.add_argument('--thresh_count', type=int, default=10,
                        help='threshold for site coverage')
    parser.add_argument('--img_out', type=str, default=home,
                        help='output directory')
    parser.add_argument('--ds', type=str, required=True,
                        help='Name')
    args = parser.parse_args()
    
    os.makedirs(args.img_out, exist_ok=True)
    
    df_gtf = get_df_gtf(args)
    junction_freq = dict()
    non_junction_freq = dict()
    for mod_code in args.mod_code:
        print(f'Processing {mod_code} for {args.bedmethyl}...')
        df_bed = get_df_bed_mod(args, mod_code)
        junction_freq[mod_code], non_junction_freq[mod_code] = get_exon_junction_freq(df_gtf, df_bed, args)
    make_plot_by_2(junction_freq, non_junction_freq, args)


if __name__ == '__main__':
    main()
    print('Finished')
