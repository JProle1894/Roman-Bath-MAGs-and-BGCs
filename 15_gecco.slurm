#!/bin/bash
#SBATCH --job-name=gecco
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --partition=cpu_shared
#SBATCH --time=24:00:00
#SBATCH --mem=64G
#SBATCH --array=0-3%4
#SBATCH --output=gecco_%A_%a.out
#SBATCH --error=gecco_%A_%a.err

set -eo pipefail

PROJECT="$HOME/MAGs_V2"
OUTPUT_ROOT="${PROJECT}/14.GECCO"

cd "$PROJECT"

# ------------------------------------------------
# Load Conda and activate GECCO
# ------------------------------------------------

module load conda/24.9.2
source /scratch/conda/miniconda3/20240912/etc/profile.d/conda.sh
conda activate gecco_env

command -v gecco >/dev/null 2>&1 || {
    echo "ERROR: GECCO is not available in gecco_env"
    exit 1
}

# ------------------------------------------------
# Samples
# ------------------------------------------------

samples=(
    "GB_Biofilm"
    "GB_Sediment"
    "HS_Biofilm"
    "HS_Sediment"
)

sample="${samples[$SLURM_ARRAY_TASK_ID]}"

# ------------------------------------------------
# Input / output
# ------------------------------------------------

assembly="${PROJECT}/02.MEGAHIT/${sample}/final.contigs.fa"
output_dir="${OUTPUT_ROOT}/${sample}"

echo "========================================"
echo "GECCO BGC prediction"
echo "Sample:   ${sample}"
echo "Assembly: ${assembly}"
echo "Output:   ${output_dir}"
echo "CPUs:     ${SLURM_CPUS_PER_TASK}"
echo "Node:     $(hostname)"
echo "Started:  $(date)"
echo "========================================"

# ------------------------------------------------
# Check input
# ------------------------------------------------

if [ ! -s "$assembly" ]; then
    echo "ERROR: MEGAHIT assembly missing or empty:"
    echo "$assembly"
    exit 1
fi

mkdir -p "$OUTPUT_ROOT"

# Avoid accidentally overwriting an existing run
if [ -d "$output_dir" ] &&
   [ -n "$(ls -A "$output_dir" 2>/dev/null)" ]; then

    echo "ERROR: GECCO output directory already contains files:"
    echo "$output_dir"
    echo "Remove or rename it before rerunning."
    exit 1
fi

mkdir -p "$output_dir"

# ------------------------------------------------
# Run GECCO
# ------------------------------------------------

echo
echo "Running GECCO..."
echo

gecco run \
    --genome "$assembly" \
    --output-dir "$output_dir" \
    --jobs "$SLURM_CPUS_PER_TASK"

# ------------------------------------------------
# Check output
# ------------------------------------------------

genes_file="${output_dir}/final.contigs.genes.tsv"
clusters_file="${output_dir}/final.contigs.clusters.tsv"

if [ ! -s "$genes_file" ]; then
    echo "ERROR: GECCO genes output was not produced."
    exit 1
fi

echo
echo "========================================"
echo "GECCO completed successfully"
echo "Sample: ${sample}"
echo "Finished: $(date)"
echo "========================================"

echo
echo "Output files:"
find "$output_dir" -maxdepth 1 -type f -printf '%f\n' | sort

echo
if [ -s "$clusters_file" ]; then

    cluster_count=$(tail -n +2 "$clusters_file" | wc -l)

    echo "Predicted BGCs: ${cluster_count}"
    echo "Cluster table: ${clusters_file}"

else

    echo "No clusters.tsv file found."
    echo "This can legitimately mean GECCO predicted no BGCs for this assembly."

fi

echo
echo "GenBank BGC files:"
find "$output_dir" \
    -maxdepth 1 \
    -type f \
    -name "*_cluster_*.gbk" \
    | wc -l