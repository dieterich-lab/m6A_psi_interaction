#!/bin/bash

# Batch submission script for Metagene Plots
# Usage: ./run_batch_metaplots.sh

# Paths
# Adjust PROJECT_ROOT as necessary
PROJECT_ROOT=$(readlink -f "$(dirname "$0")/../..")
# Use local copy of GenePred in the scripts dir
GENEPRED="$(dirname "$0")/assets/data/hg38_all_gencode_v36.genePred"
SLURM_SCRIPT="$PROJECT_ROOT/scripts/m6A_psi_interaction/run_metaplot.slurm"
OUT_BASE="$PROJECT_ROOT/results/metaplots"

# Load modules
module load python3/3.10.13_deb12 2>/dev/null || echo "Module load skipped"

echo "Project Root: $PROJECT_ROOT"
echo "Annotation: $GENEPRED"
echo "Output Base: $OUT_BASE"

mkdir -p "$OUT_BASE"

# Check if SLURM script exists
if [ ! -f "$SLURM_SCRIPT" ]; then
    echo "Error: Slurm script not found at $SLURM_SCRIPT"
    exit 1
fi

EXTRACT_SCRIPT="$PROJECT_ROOT/scripts/m6A_psi_interaction/extract_bedmethyl.py"

# --- Part 1: Process Raw Data (bedmethyl) from data/ ---
echo "--- Processing Raw Datasets from data/ ---"
DATA_DIR="$PROJECT_ROOT/data"
BEDMETHYL_FILES=$(ls "$DATA_DIR"/*.bedmethyl_.gz 2>/dev/null)

if [ -n "$BEDMETHYL_FILES" ]; then
    for bf in $BEDMETHYL_FILES; do
        if [ ! -f "$bf" ]; then continue; fi
        
        base=$(basename "$bf" .bedmethyl_.gz)
        echo "Processing $base..."
        
        # Process for both m6A ('a') and PseudoU ('17802')
        
        # define output for extracted beds
        EXTRACT_DIR="$PROJECT_ROOT/results/extracted_beds/$base"
        mkdir -p "$EXTRACT_DIR"
        
        # 1. m6A
        echo "  Extracting m6A..."
        bed_m6a="$EXTRACT_DIR/${base}.m6A.bed"
        if [ ! -f "$bed_m6a" ]; then
            python3 "$EXTRACT_SCRIPT" \
                --input "$bf" \
                --output "$bed_m6a" \
                --mod "a" \
                --min_cov 20 --min_pct 10
        fi
        
        # Submit m6A plot
        out_plot_m6a="$OUT_BASE/${base}_m6A"
        if [ ! -f "$out_plot_m6a/metagene_plot.png" ]; then
            echo "  Submitting m6A plot job..."
            sbatch "$SLURM_SCRIPT" "$bed_m6a" "$GENEPRED" "$out_plot_m6a"
        else
            echo "  [Skip] m6A plot exists."
        fi
        
        # 2. PseudoU (17802)
        echo "  Extracting PseudoU..."
        bed_psi="$EXTRACT_DIR/${base}.psi.bed"
        if [ ! -f "$bed_psi" ]; then
            python3 "$EXTRACT_SCRIPT" \
                --input "$bf" \
                --output "$bed_psi" \
                --mod "17802" \
                --min_cov 20 --min_pct 10
        fi
        
        # Submit PseudoU plot
        out_plot_psi="$OUT_BASE/${base}_psi"
        if [ ! -f "$out_plot_psi/metagene_plot.png" ]; then
            echo "  Submitting PseudoU plot job..."
            sbatch "$SLURM_SCRIPT" "$bed_psi" "$GENEPRED" "$out_plot_psi"
        else
            echo "  [Skip] PseudoU plot exists."
        fi
        
        # 3. Combined Plot Job
        # We need the metagene_data.tsv from both previous steps.
        # Create a separate SLURM script for combination that takes the two metagene_data.tsv files
        
        # The output of processing is $out_plot_m6a/metagene_data.tsv
        DATA_M6A="$out_plot_m6a/metagene_data.tsv"
        DATA_PSI="$out_plot_psi/metagene_data.tsv"
        OUT_COMBINED="$OUT_BASE/${base}_combined/combined_plot.png"
        mkdir -p "$(dirname "$OUT_COMBINED")"
        
        COMBINED_SCRIPT="$PROJECT_ROOT/scripts/m6A_psi_interaction/run_combined_plot.slurm"
        
        if [ ! -f "$OUT_COMBINED" ]; then
            echo "  Submitting Combined plot job (waits for inputs)..."
            sbatch "$COMBINED_SCRIPT" "$DATA_M6A" "$DATA_PSI" "$OUT_COMBINED"
        else
            echo "  [Skip] Combined plot exists."
        fi
        
    done
else
    echo "No .bedmethyl_.gz files found in $DATA_DIR"
fi

# --- Part 2: Process Differential Sites (from results/) ---
echo "--- Processing Differential Sites from results/ ---"
# Looking for BED files in results/
BED_FILES=$(ls "$PROJECT_ROOT/results"/sites_*.bed 2>/dev/null)

if [ -n "$BED_FILES" ]; then
    for bed in $BED_FILES; do
        name=$(basename "$bed" .bed)
        out_dir="$OUT_BASE/$name"
        
        echo "Submitting job for $name..."
        sbatch "$SLURM_SCRIPT" "$bed" "$GENEPRED" "$out_dir"
    done
else
    echo "No sites_*.bed files found in $PROJECT_ROOT/results/"
fi

# 2. Process Bedmethyl files (optional/example)
# Requires high memory/time for full files.
# files: data/*.bedmethyl_.gz

# Example:
# bedmethyl="$PROJECT_ROOT/data/HEK293_TRUB1_1.bedmethyl_.gz"
# if [ -f "$bedmethyl" ]; then
#    echo "Submitting full bedmethyl data..."
#    sbatch --mem=32G --time=4:00:00 "$SLURM_SCRIPT" "$bedmethyl" "$GENEPRED" "$OUT_BASE/HEK293_TRUB1_1"
# fi

echo "All jobs submitted."
