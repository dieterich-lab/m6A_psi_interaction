import pandas as pd
import os
from tqdm import tqdm
import matplotlib
import matplotlib.pyplot as plt
from argparse import ArgumentParser
matplotlib.use('TkAgg')


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


def get_df_bed_mod(in_args):
    print(f'Loading bedmethyl {in_args.bedmethyl}')

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
    sel_bedmethyl_fields = [bedmethyl_fields[i] for i in list(range(6)) + [10]]

    out_df_bed = pd.read_csv(in_args.bedmethyl, sep='\t', usecols=list(range(6)) + [10],
                             names=sel_bedmethyl_fields, dtype={'chrom': str})
    out_df_bed_mod = out_df_bed[out_df_bed['name'] == in_args.mod_code]
    out_df_bed_mod_strand = {
        '+': out_df_bed_mod[out_df_bed_mod['strand'] == '+'],
        '-': out_df_bed_mod[out_df_bed_mod['strand'] == '-']
    }
    return out_df_bed_mod_strand


def get_exon_junction_freq(in_df_gtf, in_df_bed, in_args):
    print(f'Parsing sites in {len(in_df_gtf)} exons...')

    out_junction_freq = []
    out_non_junction_freq = []
    for _, this_row in tqdm(in_df_gtf.iterrows()):
        this_chrom_start = this_row['start'] - 1
        this_chrom_end = this_row['end']
        this_strand = str(this_row['strand'])

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


def make_plot(in_junction_freq, in_non_junction_freq, in_args):
    out_dir = os.path.dirname(in_args.out_file)
    os.makedirs(out_dir, exist_ok=True)

    bin_range = [0, 100]
    num_bins = 10

    plt.figure(figsize=(4, 4))
    # plt.violinplot([junction_freq, non_junction_freq])
    plt.hist(in_junction_freq, density=True, log=True, range=bin_range, bins=num_bins,
             alpha=0.5, color='r', histtype='step', label=f'$\leq${in_args.exon_junction_margin} nt of exon edge')
    plt.hist(in_non_junction_freq, density=True, log=True, range=bin_range, bins=num_bins,
             alpha=0.5, color='b', histtype='step', label=f'>{in_args.exon_junction_margin} nt')
    plt.legend()
    plt.xlim(bin_range)
    plt.xlabel('Site S')
    plt.ylabel('Density')
    plt.savefig(in_args.out_file, bbox_inches='tight')


def main():
    parser = ArgumentParser()
    parser.add_argument('--gtf_exon', type=str, required=True,
                        help='annotation file for exons')
    parser.add_argument('--bedmethyl', type=str, required=True,
                        help='bedmethyl output from modkit pileup')
    parser.add_argument('--mod_code', type=str, required=True,
                        help='modification code in bedmethyl, eg, "a" for m6A, "17802" for pseudouridine')
    parser.add_argument('--exon_junction_margin', type=int, default=20,
                        help='cut-off for classification of sites on exon junction')
    parser.add_argument('--out_file', type=str, default=f'./exon_junction_mod_freq.png',
                        help='path of output image')
    args = parser.parse_args()

    df_gtf = get_df_gtf(args)
    df_bed = get_df_bed_mod(args)
    junction_freq, non_junction_freq = get_exon_junction_freq(df_gtf, df_bed, args)
    make_plot(junction_freq, non_junction_freq, args)


if __name__ == '__main__':
    main()
    print('Finished')
