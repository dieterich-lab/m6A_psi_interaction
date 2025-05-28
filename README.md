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

## Volcano plot
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

## Exon junction
This script computes the distribution of site modifications that are within / further away from a certain margin of exon edges. The following example uses the zipped data in `assets/data/chr1.exon_jcn_data.zip`:
```
python3 hist_mod_level_near_exon_junction.py \
--gtf_exon
./assets/data/chr1.exons.GRCh38.102.gtf
--bedmethyl
./assets/data/chr1.CTRL.cov10.bedmethyl
--mod_code
a
--out_file
${out_dir}/hist_exon_mod/hist_mod_chr1_a.png
```

m6A (`--mod_code a`) in chr1

![](https://github.com/ADHDrian/RNA004_psi_KD_OE_analysis/blob/main/assets/images/hist_mod_chr1_a.png)

psi (`--mod_code 17802`)

![](https://github.com/ADHDrian/RNA004_psi_KD_OE_analysis/blob/main/assets/images/hist_mod_chr1_17802.png)

The exon junction margin can be adjusted by `--exon_junction_margin`.
Note: For the purpose of speed optimization, the input exon gtf and bedmethyl files are expected to be filtered such sites are within the same chromosome. Inputting cross-chromosome data might produce erroneous results.

# Upstream pipelines with ONT software
The following steps are performed to generate bedmethyl and dmr files for the aforementioned analyses.

## Note on dorado v1.0.0
The new version of dorado implemented new modification models for 2-O-Methyl entities.
To run basecalling:
```
module load dorado/1.0.0

model=${model_dir}/rna004_130bps_sup@v5.2.0

dorado basecaller ${model} ${pod5} \
--output-dir ${outdir} \
--recursive \
--modified-bases m5C_2OmeC inosine_m6A_2OmeA pseU_2OmeU 2OmeG \
--estimate-poly-a
```
This will generate an unmapped modbam file in ${outdir}. To align it *directly* to a reference, without extracting the fastq. For example:
```
index=/biodb/genomes/homo_sapiens/GRCh38_102/GRCh38_102.mmi
junc=/biodb/genomes/homo_sapiens/GRCh38_102/GRCh38_102.bed
in_bam=/prj/TRR319_RMaP_BaseCalling_RNA004/Isabel/20250512_HEK293_M3I/dorado_v1/HEK293_M3I_24h_1/calls_*.bam

dorado aligner ${index} ${in_bam} --mm2-opts "-x splice --junc-bed ${junc} -k 14 --secondary=no" > ${out_bam}
```
According to ONT's documentation, the alignment can be done simultaneously with basecalling. However, I have not tested it myself.

## Note on modkit v0.4.2
To generate site-level modification levels from mapped modbam files:
```
module load modkit/0.4.2

BAM=${BASE}/mapped.bam
OUTPUT=${BASE}/modkit042.bedmethyl

modkit pileup -t 36 --filter-threshold 0.95 ${BAM} ${OUTPUT}
```

To perform differential modification analysis on a specific base (A, T) between two different conditions with replicates:
```
ref="/biodb/genomes/homo_sapiens/GRCh38_102/GRCh38_102.fa"

base=A
#base=T

modkit dmr pair \
  -a ${in_dir}/${cond0}_rep1.cov10.bedmethyl.gz \
  -a ${in_dir}/${cond0}_rep2.cov10.bedmethyl.gz \
  -b ${in_dir}/${cond1}_rep1.cov10.bedmethyl.gz \
  -b ${in_dir}/${cond1}_rep2.cov10.bedmethyl.gz \
  -o ${out_dir}/${cond0}_${cond1}_replicates.cov10.base${base}.dmr \
  --ref ${ref} \
  --base ${base} \
  --threads 36
```
