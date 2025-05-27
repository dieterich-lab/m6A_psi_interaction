# Analysis scripts for m6A-psi interaction
Requires packages `matplotlib` and `logomaker`

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
