import os
import pandas as pd
import matplotlib as mpl
mpl.use('TkAgg')
#######################################################################
cm = 1/2.54  # centimeters in inches
gr = 1.618
dpi = 1200
mpl.rcParams['figure.dpi'] = dpi
mpl.rcParams['savefig.dpi'] = dpi
mpl.rcParams['font.size'] = 8
mpl.rcParams['legend.fontsize'] = 6
mpl.rcParams['xtick.labelsize'] = 8
mpl.rcParams['ytick.labelsize'] = 8
mpl.rcParams['xtick.major.size'] = 4
mpl.rcParams['ytick.major.size'] = 4
mpl.rcParams['lines.linewidth'] = 1
mpl.rcParams['font.family'] = 'Arial'
# FMT = 'svg'
# fig_kwargs = dict(format=FMT, bbox_inches='tight', dpi=dpi, transparent=True)
FMT = 'png'
fig_kwargs = dict(format=FMT, bbox_inches='tight', dpi=dpi)
#######################################################################
import matplotlib.pyplot as plt


def get_longest_isoform(in_df):
    unique_locs = list(set([f"{row['chr']}_{row['coord']}" for _, row in df_metagene[['chr', 'coord']].iterrows()]))

    longest_tx = []
    for this_unique_loc in unique_locs:
        this_chr, this_coord = this_unique_loc.split('_')
        sub_df = in_df[(in_df['chr'] == this_chr) * (in_df['coord'] == int(this_coord))]
        tx_len = sub_df[['utr5_size', 'cds_size', 'utr3_size']].sum(axis=1)
        longest_tx.append(sub_df.iloc[tx_len.argmax()])
    out_df = pd.DataFrame(longest_tx)
    out_df.sort_values(['chr', 'coord'], inplace=True)
    return out_df


bed_file = '/home/adrian/Data/TRR319_RMaP_BaseCalling_RNA004/Adrian/m6a.dist.measures.txt'
img_out = '/home/adrian/img_out/RNA004_psi_KD_OE_analysis'

df_metagene = pd.read_csv(bed_file, sep='\t')
df_metagene_longest_tx = get_longest_isoform(df_metagene)

plt.figure(figsize=(5*cm, 5*cm))
plt.hist(df_metagene_longest_tx['rel_location'], range=[0, 3], bins=30)
plt.xlabel('Gene region')
plt.ylabel('Site count')
plt.xticks([])
plt.xlim([0, 3])
plt.axvline(x=1, c='gray', ls='--')
plt.axvline(x=2, c='gray', ls='--')
plt.savefig(os.path.join(img_out, f'metagene_m6a_up.{FMT}'), **fig_kwargs)