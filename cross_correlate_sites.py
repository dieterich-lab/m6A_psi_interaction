import os
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('TkAGG')
import matplotlib.pyplot as plt


img_out = '/home/adrian/img_out'

ds = 'TRUB1-KD'
cond0 = 'down_psi'
cond1 = 'up_m6A'
sites_file0 = os.path.join(img_out, f'RNA004_psi_KD_OE_analysis/sites_{ds}_{cond0}.tsv')
sites_file1 = os.path.join(img_out, f'RNA004_psi_KD_OE_analysis/sites_{ds}_{cond1}.tsv')

df_sites0 = pd.read_csv(sites_file0, sep='\t', dtype={'chrom': str})
df_sites1 = pd.read_csv(sites_file1, sep='\t', dtype={'chrom': str})

min_dist = []
for _, this_row in df_sites0.iterrows():
    chrom0, chromStart0, strand0 = this_row[['chrom', 'chromStart', 'strand']]
    sub_df_sites1 = df_sites1[
        (df_sites1['chrom'] == chrom0)
        * (df_sites1['strand'] == strand0)
    ]

    distances = (sub_df_sites1['chromStart'] - chromStart0).values
    if strand0 == '-':
        distances = -distances
    if len(distances):
        min_dist.append(distances[np.argmin(np.abs(distances))])

xmax = 10
plt.figure(figsize=(4, 4))
counts, _, _ = plt.hist(min_dist, range=[-xmax, xmax], bins=2*xmax)
plt.xlabel(f'Nearest distance (nt)\nPos({cond1}) - Pos({cond0})')
plt.ylabel('Site count')
plt.title(f'{int(np.sum(counts))} / {len(min_dist)} {cond0} sites')
plt.savefig(os.path.join(img_out, 'RNA004_psi_KD_OE_analysis', f'hist_dist_{ds}_{cond0}_{cond1}.png'), bbox_inches='tight')