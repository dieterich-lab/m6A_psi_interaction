#!/bin/bash
#SBATCH --job-name=metaplot
#SBATCH --output=metaplot_%j.out
#SBATCH --error=metaplot_%j.err
#SBATCH --mem=8G
#SBATCH --time=2:00:00
#SBATCH --cpus-per-task=1

# Usage: sbatch run_metaplot.slurm <BED_FILE> <GENEPRED_FILE> <OUT_DIR>
# Example: sbatch run_metaplot.slurm data/sample.bed data/hg38.genePred results/sample_metaplot

# Set Paths
# Assuming the script is submitted from project root or scripts directory.
# We try to locate the scripts directory relative to this file or current dir.

if [ -z "$SLURM_SUBMIT_DIR" ]; then
    WORKING_DIR=$(pwd)
else
    WORKING_DIR=$SLURM_SUBMIT_DIR
fi

# Logic to find scripts folder
if [ -f "$WORKING_DIR/scripts/m6A_psi_interaction/generate_metagene_data.py" ]; then
    SCRIPTS_DIR="$WORKING_DIR/scripts/m6A_psi_interaction"
elif [ -f "./generate_metagene_data.py" ]; then
    SCRIPTS_DIR="."
else
    echo "Could not locate scripts directory."
    exit 1
fi

# Activate Python Environment
# Adjust this loading logic to your cluster specifics
module load python3/3.10.13_deb12 2>/dev/null || echo "Module load skipped or failed"

# Check for venv
if [ -d "$WORKING_DIR/.venv" ]; then
    source "$WORKING_DIR/.venv/bin/activate"
elif [ -d "$SCRIPTS_DIR/../.venv" ]; then
    source "$SCRIPTS_DIR/../.venv/bin/activate"
fi

# Arguments
BED_FILE=$1
GENEPRED_FILE=$2
OUT_DIR=$3

if [ -z "$BED_FILE" ] || [ -z "$GENEPRED_FILE" ] || [ -z "$OUT_DIR" ]; then
    echo "Usage: sbatch run_metaplot.slurm <BED_FILE> <GENEPRED_FILE> <OUT_DIR>"
    exit 1
fi

echo "Detailed usage:"
echo "  BED File: $BED_FILE"
echo "  Annotation: $GENEPRED_FILE"
echo "  Output Directory: $OUT_DIR"
echo "  Scripts Directory: $SCRIPTS_DIR"

mkdir -p "$OUT_DIR"

# 1. Generate Metagene Data (Python implementation of MetaPlotR processing)
DATA_FILE="$OUT_DIR/metagene_data.tsv"
echo "[$(date)] Generating metagene data..."

python3 "$SCRIPTS_DIR/generate_metagene_data.py" \
    --bed "$BED_FILE" \
    --genepred "$GENEPRED_FILE" \
    --output "$DATA_FILE"

if [ $? -ne 0 ]; then
    echo "Error in generate_metagene_data.py"
    exit 1
fi

# 2. Plot Metagene
echo "[$(date)] Plotting..."

python3 "$SCRIPTS_DIR/plot_metagene.py" \
    --input "$DATA_FILE" \
    --output_dir "$OUT_DIR" \
    --plot_name "metagene_plot" \
    --use_longest_isoform

if [ $? -ne 0 ]; then
    echo "Error in plot_metagene.py"
    exit 1
fi

echo "[$(date)] Pipeline completed successfully."
