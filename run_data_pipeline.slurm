#!/bin/bash

#SBATCH -J "m6A_data_pipe"
#SBATCH -o "data_pipeline_%j.out"
#SBATCH -e "data_pipeline_%j.err"
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 36
#SBATCH --mem=50G
#SBATCH -p long

# ------------------------------------------------------------------
# Environment Setup
# ------------------------------------------------------------------
source /biosw/__modules__/modules.rc
module load python3/3.10.13_deb12
module load modkit/0.5.1_1d4b492
module load samtools
module load htslib

# ------------------------------------------------------------------
# Configuration & Variables
# ------------------------------------------------------------------
BASE_DIR="/prj/pratikum_ws_2025/project2"
DATA_DIR="${BASE_DIR}/data"
SCRIPTS_DIR="${BASE_DIR}/scripts/m6A_psi_interaction"
RESULTS_DIR="${BASE_DIR}/scripts/results_ajili/final_pipeline_results"

# Reference Files
REF_GENOME="/biodb/genomes/homo_sapiens/GRCh38_102/GRCh38_102.fa"

# Input BAMs
BAM_CTRL_1="${DATA_DIR}/HEK293_ctrl_1_.bam"
BAM_CTRL_2="${DATA_DIR}/HEK293_ctrl_2_.bam"
BAM_TRUB1_1="${DATA_DIR}/HEK293_TRUB1_1_.bam"
BAM_TRUB1_2="${DATA_DIR}/HEK293_TRUB1_2_.bam"

# Create Output Directories
mkdir -p "${RESULTS_DIR}/bedmethyl"
mkdir -p "${RESULTS_DIR}/dmr"
mkdir -p "${RESULTS_DIR}/merged_bam"

# ------------------------------------------------------------------
# Helper Functions
# ------------------------------------------------------------------
log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1"
}

run_pileup() {
    local bam_in=$1
    local name=$2
    local threads=$3
    local out_bed="${RESULTS_DIR}/bedmethyl/${name}.bedmethyl.gz"

    # if [ -f "${out_bed}" ] && [ -f "${out_bed}.tbi" ]; then
    #     log "Skipping pileup for ${name}, output exists."
    # else
    log "Running modkit pileup for ${name} with ${threads} threads..."
    modkit pileup -t 36 --filter-threshold 0.95 "${bam_in}" - | bgzip -c > "${out_bed}"
    tabix -p bed "${out_bed}"
    # fi
}

run_dmr() {
    local base=$1
    local outfile=$2
    local threads=$3
    
    if [ -f "${outfile}" ]; then
        log "Skipping DMR Base ${base}, output exists."
    else
        log "Running modkit dmr for Base ${base} with ${threads} threads..."
        modkit dmr pair \
            -b "${RESULTS_DIR}/bedmethyl/HEK293_TRUB1_1.bedmethyl.gz" \
            -b "${RESULTS_DIR}/bedmethyl/HEK293_TRUB1_2.bedmethyl.gz" \
            -a "${RESULTS_DIR}/bedmethyl/HEK293_ctrl_1.bedmethyl.gz" \
            -a "${RESULTS_DIR}/bedmethyl/HEK293_ctrl_2.bedmethyl.gz" \
            -o "${outfile}" \
            --header \
            --ref "${REF_GENOME}" \
            --base "${base}" \
            --threads "${threads}"
    fi
}

# ------------------------------------------------------------------
# Step 1: Generate Bedmethyl Files (Pileup)
# ------------------------------------------------------------------
log "Starting Phase 1: Primary Data Processing (Pileup)"

# Run 4 pileups in parallel. 
run_pileup "${BAM_CTRL_1}" "HEK293_ctrl_1" 36 &
# run_pileup "${BAM_CTRL_2}" "HEK293_ctrl_2" 9 &
# run_pileup "${BAM_TRUB1_1}" "HEK293_TRUB1_1" 9 &
# run_pileup "${BAM_TRUB1_2}" "HEK293_TRUB1_2" 9 &

wait
log "Phase 1 finished."

# ------------------------------------------------------------------
# Step 2: Merge BAMs for Cross-Talk Analysis
# ------------------------------------------------------------------
log "Starting Phase 2: Merging BAMs"

MERGED_TRUB1_BAM="${RESULTS_DIR}/merged_bam/HEK293_TRUB1_merged.bam"

if [ -f "${MERGED_TRUB1_BAM}" ]; then
    log "Skipping BAM merge, output exists."
else
    log "Merging TRUB1 BAM files..."
    # This is a single heavy job, give it all cores
    samtools merge -@ 36 -o - "${BAM_TRUB1_1}" "${BAM_TRUB1_2}" | samtools sort -@ 36 -o "${MERGED_TRUB1_BAM}" -
    samtools index "${MERGED_TRUB1_BAM}"
fi

# ------------------------------------------------------------------
# Step 3: Differential Modification Analysis (DMR)
# ------------------------------------------------------------------
log "Starting Phase 3: DMR Analysis"

DMR_A="${RESULTS_DIR}/dmr/HEK293_ctrl_vs_TRUB1.baseA.dmr"
DMR_T="${RESULTS_DIR}/dmr/HEK293_ctrl_vs_TRUB1.baseT.dmr"

# Run 2 DMR analyses in parallel.
# Total cores = 36. 36 / 2 = 18 threads per job.
run_dmr "A" "${DMR_A}" 18 &
run_dmr "T" "${DMR_T}" 18 &

wait
log "Phase 3 finished."
log "Data Pipeline finished successfully! You can now run the plotting pipeline."
