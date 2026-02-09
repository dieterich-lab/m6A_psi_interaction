#!/bin/bash
#SBATCH --job-name=combined_plot
#SBATCH --output=combined_plot_%j.out
#SBATCH --error=combined_plot_%j.err
#SBATCH --mem=8G
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=1

# Usage: sbatch run_combined_plot.slurm <M6A_TSV> <PSI_TSV> <OUT_PNG>

if [ -z "$SLURM_SUBMIT_DIR" ]; then
    WORKING_DIR=$(pwd)
else
    WORKING_DIR=$SLURM_SUBMIT_DIR
fi

# Logic to find scripts folder
if [ -f "$WORKING_DIR/scripts/m6A_psi_interaction/plot_combined_metagene.py" ]; then
    SCRIPTS_DIR="$WORKING_DIR/scripts/m6A_psi_interaction"
else
    # Fallback or assume relative
    SCRIPTS_DIR="."
fi

# Activate Python Environment
module load python3/3.10.13_deb12 2>/dev/null || echo "Module load skipped"

if [ -d "$WORKING_DIR/.venv" ]; then
    source "$WORKING_DIR/.venv/bin/activate"
elif [ -d "$SCRIPTS_DIR/../.venv" ]; then
    source "$SCRIPTS_DIR/../.venv/bin/activate"
fi

M6A_IN=$1
PSI_IN=$2
OUT_FILE=$3

if [ -z "$M6A_IN" ] || [ -z "$PSI_IN" ] || [ -z "$OUT_FILE" ]; then
    echo "Usage: sbatch run_combined_plot.slurm <M6A_TSV> <PSI_TSV> <OUT_PNG>"
    exit 1
fi

echo "Detailed usage:"
echo "  M6A Input: $M6A_IN"
echo "  PSI Input: $PSI_IN"
echo "  Output File: $OUT_FILE"

# Wait for files to exist (simple polling if submitted together)
MAX_WAIT=3600
elapsed=0
while [ ! -f "$M6A_IN" ] || [ ! -f "$PSI_IN" ]; do
    if [ $elapsed -ge $MAX_WAIT ]; then
        echo "Timed out waiting for input files."
        exit 1
    fi
    echo "Waiting for input files ($elapsed / $MAX_WAIT)..."
    sleep 60
    elapsed=$((elapsed+60))
done

echo "[$(date)] plotting combined..."

python3 "$SCRIPTS_DIR/plot_combined_metagene.py" \
    --m6a_input "$M6A_IN" \
    --psi_input "$PSI_IN" \
    --output "$OUT_FILE" \
    --use_longest_isoform

if [ $? -ne 0 ]; then
    echo "Error in plot_combined_metagene.py"
    exit 1
fi

echo "[$(date)] Done."
