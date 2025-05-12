import pandas as pd
import os
from tqdm import tqdm
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt


home = os.environ['HOME']
img_out = f'{home}/img_out/RNA004_psi_KD_OE_analysis'
os.makedirs(img_out, exist_ok=True)

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

bedMethyl_fields = [
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
sel_bedMethyl_fields = [bedMethyl_fields[i] for i in list(range(6))+[10]]

dict_mod_display = {
    'a': 'm^6A',
    '17802': '\Psi'
}

gtf_dir = f'{home}/Data/genomes/homo_sapiens/GRCh38_102/exons'
bed_dir = f'{home}/Data/TRR319_RMaP_BaseCalling_RNA004/Adrian/HEK293_psU-KD'

this_chr = 'X'
mod_code = '17802'
# mod_code = 'a'
this_gtf_file = os.path.join(gtf_dir, f'chr{this_chr}.exons.GRCh38.102.gtf')
df_gtf = pd.read_csv(this_gtf_file, sep='\t', names=gtf_fields, dtype={'chrom': str})
df_gtf = df_gtf[df_gtf['source'] == 'ensembl_havana']

this_bed_file = os.path.join(bed_dir, f'chr{this_chr}.CTRL.modkit042.cov10.bedmethyl')
df_bed = pd.read_csv(this_bed_file, sep='\t', usecols=list(range(6))+[10],
                     names=sel_bedMethyl_fields, dtype={'chrom': str})
df_bed_mod = df_bed[df_bed['name'] == mod_code]
df_bed_mod_strand = {
    '+': df_bed_mod[df_bed_mod['strand'] == '+'],
    '-': df_bed_mod[df_bed_mod['strand'] == '-']
}

exon_junction_margin = 20

junction_freq = []
non_junction_freq = []
print(f'Parsing sites in {len(df_gtf)} exons...')
for _, this_row in tqdm(df_gtf.iterrows()):
    # this_chrom = this_row['chrom']
    this_chrom_start = this_row['start'] - 1
    this_chrom_end = this_row['end']
    this_strand = str(this_row['strand'])

    this_strand_bed =  df_bed_mod_strand[this_strand]
    sub_df_bed = this_strand_bed[
        (this_strand_bed['chromStart'] >= this_chrom_start)
        * (this_strand_bed['chromEnd'] <= this_chrom_end)
    ]

    if len(sub_df_bed):
        mask = (sub_df_bed['chromStart'] <= (this_chrom_start+exon_junction_margin))\
               + (sub_df_bed['chromStart'] >= (this_chrom_end-exon_junction_margin))
        junction_freq.extend(sub_df_bed[mask]['frequency'].values)
        non_junction_freq.extend(sub_df_bed[~mask]['frequency'].values)

bin_range = [0, 100]
num_bins = 10

plt.figure(figsize=(4, 4))
# plt.violinplot([junction_freq, non_junction_freq])
plt.hist(junction_freq, density=True, log=True, range=bin_range, bins=num_bins,
         alpha=0.5, color='r', histtype='step', label=f'$\leq${exon_junction_margin} nt of exon edge')
plt.hist(non_junction_freq, density=True, log=True, range=bin_range, bins=num_bins,
         alpha=0.5, color='b', histtype='step', label=f'>{exon_junction_margin} nt')
plt.legend()
plt.xlim(bin_range)
plt.xlabel('Site S')
plt.ylabel('Density')
plt.title(f'chr{this_chr}, ${dict_mod_display[mod_code]}$')
plt.savefig(os.path.join(img_out, f'hist_junction_vs_non_junction_sites_mod{mod_code}_chr{this_chr}.png'), bbox_inches='tight')