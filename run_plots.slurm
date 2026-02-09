#!/bin/bash

#SBATCH -J "m6A_plots"
#SBATCH -o "plots_%j.out"
#SBATCH -e "plots_%j.err"
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 12
#SBATCH --mem=16G
#SBATCH -p long

# ------------------------------------------------------------------
# Environment Setup
# ------------------------------------------------------------------
source /biosw/__modules__/modules.rc
module load python3/3.10.13_deb12
module load samtools
module load htslib

# ------------------------------------------------------------------
# Configuration & Variables
# ------------------------------------------------------------------
BASE_DIR="/prj/pratikum_ws_2025/project2"
SCRIPTS_DIR="${BASE_DIR}/scripts/m6A_psi_interaction"
RESULTS_DIR="${BASE_DIR}/scripts/results_ajili/final_pipeline_results"

# Activate Python Virtual Environment
source "${SCRIPTS_DIR}/.venv/bin/activate"

# Reference Files
REF_GENOME="/biodb/genomes/homo_sapiens/GRCh38_102/GRCh38_102.fa"
GTF_FILE="/biodb/genomes/homo_sapiens/GRCh38_102/GRCh38.102.gtf"

# Input Files (from Data Folder)
DATA_DIR="${BASE_DIR}/data"
DMR_A="${DATA_DIR}/HEK293_ctrl_HEK293_TRUB1_replicates.cov10.baseA_.dmr"
DMR_T="${DATA_DIR}/HEK293_ctrl_HEK293_TRUB1_replicates.cov10.baseT_.dmr"
# Using Replicate 1 as merged BAM is not available in data folder
MERGED_TRUB1_BAM="${DATA_DIR}/HEK293_TRUB1_1_.bam"

# Create Output Directories
mkdir -p "${RESULTS_DIR}/plots"

# ------------------------------------------------------------------
# Helper Functions
# ------------------------------------------------------------------
log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1"
}

# ------------------------------------------------------------------
# Step 4: Downstream Python Analysis
# ------------------------------------------------------------------
log "Starting Phase 4: Python Scripts Analysis"

# Run independent Python scripts in parallel
(
    log "Running volcano_plot_mod_level_changes.py..."
    python3 "${SCRIPTS_DIR}/volcano_plot_mod_level_changes.py" \
        --dmr_base_A "${DMR_A}" \
        --dmr_base_U "${DMR_T}" \
        --reference "${REF_GENOME}" \
        --out_dir "${RESULTS_DIR}/plots" \
        --plot_name "volcano_plot.png" \
        --thresh_count 10 \
        --balanced
) &

(
    log "Running hist_mod_level_near_exon_junction.py (m6A and psi combined)..."
    # Run for all samples
    for bed_file in "${DATA_DIR}"/*.bedmethyl_.gz; do
        sample_name=$(basename "${bed_file}" .bedmethyl_.gz)
        log "Processing ${sample_name}..."
        
        python3 "${SCRIPTS_DIR}/hist_mod_level_near_exon_junction.py" \
            --gtf_exon "${GTF_FILE}" \
            --bedmethyl "${bed_file}" \
            --mod_code  17802 a \
            --name "${sample_name}_combined" \
            --img_out "${RESULTS_DIR}/plots" &
    done
    wait
) &

(
    log "Running single_read_m6A_psi_cross_talk.py..."
    python3 "${SCRIPTS_DIR}/single_read_m6A_psi_cross_talk.py" \
        --bam_file "${MERGED_TRUB1_BAM}" \
        --ds "TRUB1_1" \
        --img_out "${RESULTS_DIR}/plots" \
        --filter_by_motifs \
        --do_both_strands
) &

wait
log "Phase 4 initial scripts finished."

# 4.4 Cross-Correlation of Sites
# This depends on the output of volcano_plot, so it must run after the wait.
BED_UP_A="${RESULTS_DIR}/plots/sites_up_A.bed"
BED_UP_U="${RESULTS_DIR}/plots/sites_up_U.bed"

if [ -f "${BED_UP_A}" ] && [ -f "${BED_UP_U}" ]; then
    log "Running cross_correlate_sites.py..."
    python3 "${SCRIPTS_DIR}/cross_correlate_sites.py" \
        --bed1 "${BED_UP_A}" \
        --bed2 "${BED_UP_U}" \
        --output_file "${RESULTS_DIR}/plots/A_U_up_correlation.png"
else
    log "WARNING: BED files for cross-correlation not found. Skipping this step."
fi

log "Plotting Pipeline finished successfully!"
