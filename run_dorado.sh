#!/usr/bin/env bash
#SBATCH --partition=gpu
#SBATCH --gres=gpu:hopper:1
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=8GB
#SBATCH --verbose
#SBATCH --job-name=dorado_v1
#SBATCH --output=/home/achan/slurm/dorado_v1_%A.out
#SBATCH --error=/home/achan/slurm/dorado_v1_%A.err

set -e -u

module load dorado/1.0.0

model=/home/achan/dorado_models/rna004_130bps_sup@v5.2.0
pod5=/prj/TRR319_RMaP_BaseCalling_RNA004/Isabel/20250512_HEK293_M3I/HEK293_M3I_24h_1_RTA/20250512_1509_3C_PAY70432_d02d7a17/pod5
outdir=/prj/TRR319_RMaP_BaseCalling_RNA004/Isabel/20250512_HEK293_M3I/dorado_v1/HEK293_M3I_24h_1

dorado basecaller ${model} ${pod5} \
--output-dir ${outdir} \
--recursive \
--modified-bases m5C_2OmeC inosine_m6A_2OmeA pseU_2OmeU 2OmeG \
--estimate-poly-a
