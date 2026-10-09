#!/bin/bash
#SBATCH --job-name=metasmash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --partition=cpu_shared
#SBATCH --time=48:00:00
#SBATCH --mem=64G
#SBATCH --array=0-3%4
#SBATCH --output=metasmash_%A_%a.out
#SBATCH --error=metasmash_%A_%a.err

set -eo pipefail

PROJECT="$HOME/MAGs_V2"
OUTPUT_ROOT="${PROJECT}/13.MetaSMASH"

cd "$PROJECT"

# Load Conda
module load conda/24.9.2
source /scratch/conda/miniconda3/20240912/etc/profile.d/conda.sh

# Activate metaSMASH environment
conda activate metasmash

# Check metaSMASH is available
command -v metasmash >/dev/null 2>&1 || {
    echo "ERROR: metaSMASH is not available in the metasmash environment."
    exit 1
}

# -----------------------------
# Samples
# -----------------------------

samples=(
    "GB_Biofilm"
    "GB_Sediment"
    "HS_Biofilm"
    "HS_Sediment"
)

sample="${samples[$SLURM_ARRAY_TASK_ID]}"

# -----------------------------
# Input and output
# -----------------------------

assembly="${PROJECT}/02.MEGAHIT/${sample}/final.contigs.fa"
output_dir="${OUTPUT_ROOT}/${sample}"

echo "========================================"
echo "metaSMASH metagenome BGC analysis"
echo "Sample:     ${sample}"
echo "Assembly:   ${assembly}"
echo "Output:     ${output_dir}"
echo "CPUs:       ${SLURM_CPUS_PER_TASK}"
echo "Started:    $(date)"
echo "Node:       $(hostname)"
echo "========================================"

# Check assembly exists
if [ ! -s "$assembly" ]; then
    echo "ERROR: Assembly missing or empty:"
    echo "$assembly"
    exit 1
fi

# Create main output folder
mkdir -p "$OUTPUT_ROOT"

# Prevent accidental overwrite of previous results
if [ -d "$output_dir" ] &&
   [ -n "$(ls -A "$output_dir" 2>/dev/null)" ]; then

    echo "ERROR: Output directory already contains files:"
    echo "$output_dir"
    echo
    echo "Remove or rename it before rerunning."
    exit 1
fi

# -----------------------------
# Run metaSMASH
# -----------------------------

echo
echo "Running metaSMASH..."
echo

metasmash \
    --streaming on \
    --cpus "$SLURM_CPUS_PER_TASK" \
    --workers 4 \
    --taxon bacteria \
    --genefinding-tool prodigal-m \
    --cb-knownclusters \
    --output-dir "$output_dir" \
    "$assembly"

# -----------------------------
# Check completion
# -----------------------------

if [ ! -d "$output_dir" ]; then
    echo "ERROR: metaSMASH output directory was not produced."
    exit 1
fi

echo
echo "========================================"
echo "metaSMASH completed"
echo "Sample:   ${sample}"
echo "Finished: $(date)"
echo "========================================"

echo
echo "Output files:"
find "$output_dir" \
    -maxdepth 2 \
    -type f \
    | head -n 50

echo
echo "GenBank/BGC region files:"
find "$output_dir" \
    -type f \
    \( -name "*.gbk" -o -name "*.gbff" \) \
    | wc -l