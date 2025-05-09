#!/bin/bash

#SBATCH --job-name=modkit_dmr
#SBATCH --partition=general
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
#SBATCH --mem=80GB

in_dir="/prj/TRR319_RMaP_BaseCalling_RNA004/Adrian/HEK293_psU-KD"
out_dir="/prj/TRR319_RMaP_BaseCalling_RNA004/Adrian/HEK293_psU-KD"
ref="/biodb/genomes/homo_sapiens/GRCh38_102/GRCh38_102.fa"

cond0=CTRL
cond1=TRUB1-KD

mod=m6A
base=A
#mod=psi
#base=T

modkit dmr pair \
  -a ${in_dir}/${cond0}.modkit042.cov10.bedmethyl.gz \
  -b ${in_dir}/${cond1}.modkit042.cov10.bedmethyl.gz \
  -o ${out_dir}/${cond0}_${cond1}.cov10.${mod}.dmr \
  --ref ${ref} \
  --base ${base} \
  --threads 36 \
  --log-filepath ${out_dir}/${cond0}_${cond1}.cov10.${mod}.log
