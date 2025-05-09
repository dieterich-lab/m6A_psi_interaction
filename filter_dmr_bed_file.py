import os
import pandas as pd
from tqdm import tqdm
from argparse import ArgumentParser


home = os.environ['HOME']
parser = ArgumentParser()
parser.add_argument('in_file')
args = parser.parse_args()


dmr_fields = [
    'chrom',
    'chromStart',
    'chromEnd',
    'name',
    'score',
    'strand',
    'sample_a_counts',
    'sample_a_total',
    'sample_b_counts',
    'sample_b_total',
    'sample_a_percents',
    'sample_b_percents',
    'sample_a_fraction_modified',
    'sample_b_fraction_modified',
    'p_value',
    'effect_size',
    'balanced_p_value',
    'balanced_effect_size',
    'pct_a_samples',
    'pct_b_samples',
    'per_replicate_p_values',
    'per_replicate_effect_sizes'
]


in_file = args.in_file

print(f'Parsing {in_file}...')
all_chunks_filtered = []
for this_chunk in tqdm(pd.read_csv(in_file, sep='\t', names=dmr_fields, dtype={'chrom': str}, iterator=True, chunksize=10000)):
    this_chunk_filtered = this_chunk[(
            (this_chunk['sample_a_fraction_modified'] > 0) | (this_chunk['sample_b_fraction_modified'] > 0)
    )]
    all_chunks_filtered.append(this_chunk_filtered)
df_dmr_filtered = pd.concat(all_chunks_filtered)
df_dmr_filtered.to_csv(in_file.replace('.dmr', '.diff_sites.dmr'), sep='\t', index=False)
print('Done')