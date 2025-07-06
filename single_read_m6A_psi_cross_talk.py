import os
import pandas as pd
import pysam
import numpy as np
from tqdm import tqdm


thresh_valid_reads = 1000
num_top_locs = 5

def get_mean_logit(in_probs):
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

thresh_min_locs = 10
def get_mod_mean_occupancy(in_read, min_locs=thresh_min_locs):
    mod_mean_occupancy = {}
    for this_mod, this_tag in mod_tags.items():
        this_mod_locs = [this_tup[0] for this_tup in in_read.modified_bases.get(this_tag, [])]
        if len(this_mod_locs) < min_locs:
            mod_mean_occupancy[this_mod] = np.nan
            continue
        loc_motifs = [in_read.query_sequence[this_loc - 2:this_loc + 3] for this_loc in this_mod_locs]
        motif_included = np.array([this_motif in mod_motifs[this_mod] for this_motif in loc_motifs])
        this_mod_probs = np.array([this_tup[1] for this_tup in in_read.modified_bases.get(this_tag, [])]) / 255.0
        filtered_mod_probs = this_mod_probs[motif_included]
        if len(filtered_mod_probs) >= min_locs:
            mod_mean_occupancy[this_mod] = np.mean(filtered_mod_probs >= 0.5)
        else:
            mod_mean_occupancy[this_mod] = np.nan
    return mod_mean_occupancy

########################################################################################################################
### R002 ###############################################################################################################
########################################################################################################################
# mod_tags = {
#     'm6A': ('N', 0, 21891),
#     'psi': ('N', 0, 17802)
# }
#
# base_dir = '/home/adrian/Data/TRR319_RMaP_BaseCalling/Adrian/results/psico-mAFiA_v1'
#
# ds = 'HEK293_WT_R002'
# bam_file = os.path.join(base_dir, 'HEK293/WT_P2/chrALL.mAFiA.reads.bam')

# ds = 'M3KO'
# bam_file = os.path.join(base_dir, 'HEK293T_Mettl3_KO/merged/chrALL.mAFiA.reads.bam')

# ds = 'M3KD'
# bam_file = os.path.join(base_dir, 'NanoSPA/HEK_siMETTL3_input_merged/chrALL.mAFiA.reads.bam')

# ds = 'TRUB1KD'
# bam_file = os.path.join(base_dir, 'NanoSPA/HEK_siTRUB1_input_merged/chrALL.mAFiA.reads.bam')

# ds = 'TRUB1OE'
# bam_file = os.path.join(base_dir, 'HEK293_TRUB1_OE/merged/chrALL.mAFiA.reads.bam')

########################################################################################################################
### R004 ###############################################################################################################
########################################################################################################################
mod_tags = {
    'm6A': ('A', 0, 'a'),
    'psi': ('T', 0, 17802)
}

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

base_dir = '/home/adrian/Data/TRR319_RMaP_BaseCalling_RNA004/Isabel/20250224_HEK293_psU_kds_RTA/Dorado_082'

ds = 'HEK293_ctrl_R004'
bam_file = os.path.join(base_dir, 'HEK293_ctrl_RTA/calls_2025-02-26_T06-44-51.bam')

# ds = 'HEK293_TRUB1_kd'
# bam_file = os.path.join(base_dir, 'HEK293_TRUB1_kd_RTA/calls_2025-02-26_T06-43-59.bam')

########################################################################################################################

single_read_mean_occupancy = []
with pysam.AlignmentFile(bam_file, 'rb', check_sq=False) as bam:
    for this_read in tqdm(bam.fetch(until_eof=True)):
        if this_read.modified_bases is not None:
            # single_read_mean_logit.append(get_mean_logit_mod_level(this_read))
            single_read_mean_occupancy.append(get_mod_mean_occupancy(this_read))

vec_m6A, vec_psi = np.vstack([
    (this_read_mean_occupancy['m6A'],  this_read_mean_occupancy['psi'])
    for this_read_mean_occupancy in single_read_mean_occupancy
    if ~np.isnan(this_read_mean_occupancy['m6A']) and ~np.isnan(this_read_mean_occupancy['psi'])
]).T
num_valid_reads = len(vec_m6A)

out_file = os.path.join(base_dir, f'single_read_frac_m6A_psi_{ds}.npz')
np.savez(out_file, vec_m6A=vec_m6A, vec_psi=vec_psi)
