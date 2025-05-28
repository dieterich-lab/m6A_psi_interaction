# Analysis scripts for m6A-psi interaction
Required packages (version numbers not critical):
```
biopython==1.80
logomaker==0.8
matplotlib==3.10.3
numpy==1.23.5
pandas==2.2.3
pybedtools==0.10.0
pysam==0.22.0
scipy==1.8.1
tqdm==4.64.1
```

## Volcano plots
Uses output of `modkit dmr pair` to compare sites that are up- or down-regulated in their modification levels.

Basic example:
```
python3 volcano_plot_mod_level_changes.py \
--dmr_base_A /prj/TRR319_RMaP_BaseCalling_RNA004/Adrian/HEK293_M3I/dmr/DMSO_rep1_M3I_6h_rep1.cov10.m6A.diff_sites.dmr \
--dmr_base_U /prj/TRR319_RMaP_BaseCalling_RNA004/Adrian/HEK293_M3I/dmr/DMSO_rep1_M3I_6h_rep1.cov10.psi.diff_sites.dmr \
--out_dir ${out_dir}
```
will create the following volcano plot in the output directory specified.

![](https://github.com/ADHDrian/RNA004_psi_KD_OE_analysis/blob/main/assets/images/volcano_plot.png)

Other options:
- `--thresh_neg_log_pval`: (negative log10) p-value threshold for statistically significant changes, default 2.0
- `--thresh_change`: threshold for absolute change in modification percentage, default 25.0
- `--thresh_count`: minimum coverage cut-off for sites, default 50
- `--filter_by_motifs`: filters sites by preset motifs, DRACH for A and TRUB1 / PUS1 / PUS7 motifs for U; need to supply reference (see next)
- `--reference`: fasta reference file used to align the sites, for example, `/biodb/genomes/homo_sapiens/GRCh38_102/GRCh38_102.fa`
- `--output_logo`: outputs the sequence logo as well as reference locations of up / down sites, needs reference
- `--balanced`: uses sample-size adjusted p-values and changes, only available for DMR output run on replicates
