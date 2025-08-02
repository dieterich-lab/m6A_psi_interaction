import os
import pandas as pd

base_dir = "/home/achan/prj/TRR319_RMaP_BaseCalling/Adrian"
res_dir = os.path.join(
    base_dir,
    "results/psico-mAFiA_v1/rRNA/Parental_HCT116_Christiane_Zorbas/NR_003286_RNA18SN5"
)
mafia_path = os.path.join(res_dir, "mAFiA.sites.bed")
annot_path = os.path.join(base_dir, "rRNA/rnam_mod_levels.bed")

df_mafia = pd.read_csv(mafia_path, sep="\t")
df_mafia['name'] = df_mafia['name'].apply(lambda s: s.replace('psi', 'psU'))
df_mafia.drop(columns=['strand', 'score'], inplace=True)

df_annot = pd.read_csv(annot_path, sep="\t",
                       names=['chrom', 'chromStart', 'chromEnd', 'name', 'score']
                       )
df_annot = df_annot[df_annot['chrom'] == 'NR_003286_RNA18SN5']
df_annot = df_annot[df_annot['name'] != 'ND']
df_annot['score'] = df_annot['score'].apply(lambda s: float(s))

df_merged = pd.merge(df_annot, df_mafia,
                     on=['chrom', 'chromStart', 'chromEnd'],
                     how='outer',
                     suffixes=('_annot', '_mafia')
                     )

thresh_mod_ratio = 50.0
thresh_confidence = 50.0

mod = 'm6A'
# mod = 'psU'

out_file = os.path.join(
    res_dir,
    f"{mod}.annot.mAFiA.modRatio{thresh_mod_ratio}.conf{thresh_confidence}.tsv"
)

df_mod = df_merged[
    (
        (df_merged['name_annot'] == mod)
        * (df_merged['score'] >= thresh_mod_ratio)
    )
    | (
        (df_merged['name_mafia'] == mod)
        * (df_merged['modRatio'] >= thresh_mod_ratio)
        * (df_merged['confidence'] >= thresh_confidence)
    )
]
# df_mod_thresh = df_mod[
#     (df_mod['modRatio'] >= thresh_mod_ratio)
#     * (df_mod['confidence'] >= thresh_confidence)
# ]
num_sites_annot = (df_mod['name_annot'] == mod).sum()
num_sites_annot_mafia = ((df_mod['name_annot'] == mod) * (df_mod['name_mafia'] == mod)).sum()
num_sites_detected = (
        (df_mod['name_annot'] == mod)
        * (df_mod['name_mafia'] == mod)
        * (df_mod['score'] >= thresh_mod_ratio)
        * (df_mod['modRatio'] >= thresh_mod_ratio)
        * (df_mod['confidence'] >= thresh_confidence)
).sum()
df_mod.to_csv(out_file, sep='\t', index=False)
with open(out_file, 'a+') as handle_out:
    handle_out.write("\n")
    handle_out.write(f"# {num_sites_annot} annotated {mod} site(s), "
                     f"of which {num_sites_annot_mafia} within mAFiA motifs, "
                     f"{num_sites_detected} detected\n")