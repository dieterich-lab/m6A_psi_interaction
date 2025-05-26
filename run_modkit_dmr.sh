#!/bin/bash

#SBATCH --job-name=modkit_dmr
#SBATCH --partition=general
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
#SBATCH --mem=80GB
#SBATCH --output=/home/achan/slurm/modkit_dmr_%A.out

module load modkit/0.4.2

in_dir="/prj/TRR319_RMaP_BaseCalling_RNA004/Adrian/HEK293_M3I/modkit042"
out_dir="/prj/TRR319_RMaP_BaseCalling_RNA004/Adrian/HEK293_M3I/dmr"
ref="/biodb/genomes/homo_sapiens/GRCh38_102/GRCh38_102.fa"

cond0=DMSO_rep1
cond1=M3I_6h_rep1

#mod=m6A
#base=A
mod=psi
base=T

modkit dmr pair \
  -a ${in_dir}/${cond0}.cov10.bedmethyl.gz \
  -b ${in_dir}/${cond1}.cov10.bedmethyl.gz \
  -o ${out_dir}/${cond0}_${cond1}.cov10.${mod}.dmr \
  --ref ${ref} \
  --assign-code 19227:T \
  --assign-code 19228:C \
  --assign-code 19229:G \
  --assign-code 69426:A \
  --base ${base} \
  --threads 36 \
  --log-filepath ${out_dir}/${cond0}_${cond1}.cov10.${mod}.log
