#!/usr/bin/env bash
#SBATCH --partition=general
#SBATCH --cpus-per-task=40
#SBATCH --mem=100GB
#SBATCH --verbose
#SBATCH --job-name=extract_polyA
#SBATCH --output=/home/achan/slurm/extract_polyA_%A.out

in_dir="/prj/TRR319_RMaP_BaseCalling_RNA004/Isabel/20250224_HEK293_psU_kds_RTA/Dorado_082"
out_dir="/prj/TRR319_RMaP_BaseCalling_RNA004/Adrian/polyA"

in_ds=HEK293_ctrl_RTA
out_ds=HEK293_CTRL

in_bam=${in_dir}/${in_ds}/mapped.bam
out_tsv=${out_dir}/read_polyA_${out_ds}.tsv

samtools view -@36 ${in_bam} | awk '{for (i=12; i<=NF; ++i) { if ($i ~ "^pt:"){ split($i, tc, ":"); td[tc[1]] = tc[3]; } }; print $1"\t"td["pt"] }' > ${out_tsv}