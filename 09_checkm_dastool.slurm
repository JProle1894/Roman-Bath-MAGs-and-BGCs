#!/bin/bash
#SBATCH --job-name=checkm_dastool
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --partition=cpu_shared
#SBATCH --time=24:00:00
#SBATCH --mem=64G
#SBATCH --array=0-3%4
#SBATCH --output=checkm_dastool_%A_%a.out
#SBATCH --error=checkm_dastool_%A_%a.err

set -eo pipefail

PROJECT="$HOME/MAGs_V2"
cd "$PROJECT"

module load conda/24.9.2
source /scratch/conda/miniconda3/20240912/etc/profile.d/conda.sh
conda activate checkm_env

command -v checkm >/dev/null 2>&1 || {
    echo "ERROR: CheckM is not available in checkm_env"
    exit 1
}

samples=(
    "GB_Biofilm"
    "GB_Sediment"
    "HS_Biofilm"
    "HS_Sediment"
)

sample="${samples[$SLURM_ARRAY_TASK_ID]}"

# DAS Tool refined bins
input_dir="08.DASTool/${sample}/${sample}_DASTool_bins"

# CheckM output
output_dir="09.CheckM_DASTool/${sample}"
work_dir="${output_dir}/work"
results_file="${output_dir}/checkm_results.tsv"

echo "========================================"
echo "CheckM assessment of DAS Tool bins"
echo "Sample:  ${sample}"
echo "Input:   ${input_dir}"
echo "Output:  ${output_dir}"
echo "Started: $(date)"
echo "========================================"

# ----------------------------------------
# Check DAS Tool bins exist
# ----------------------------------------

if [ ! -d "$input_dir" ]; then
    echo "ERROR: DAS Tool bins directory not found:"
    echo "$input_dir"
    exit 1
fi

# ----------------------------------------
# Detect FASTA extension
# ----------------------------------------

bin_extension=""

for ext in fa fasta fna; do

    count=$(find "$input_dir" \
        -maxdepth 1 \
        -type f \
        -name "*.${ext}" \
        | wc -l)

    if [ "$count" -gt 0 ]; then
        bin_extension="$ext"
        bin_count="$count"
        break
    fi

done

if [ -z "$bin_extension" ]; then
    echo "ERROR: No FASTA bins found in:"
    echo "$input_dir"
    exit 1
fi

echo "Bins found: ${bin_count}"
echo "Extension: .${bin_extension}"

# ----------------------------------------
# Prevent accidental overwrite
# ----------------------------------------

if [ -d "$output_dir" ] &&
   [ -n "$(ls -A "$output_dir" 2>/dev/null)" ]; then

    echo "ERROR: CheckM output directory already contains files:"
    echo "$output_dir"
    exit 1

fi

mkdir -p "$output_dir"

# ----------------------------------------
# Run CheckM
# ----------------------------------------

checkm lineage_wf \
    -x "$bin_extension" \
    -t "$SLURM_CPUS_PER_TASK" \
    --tab_table \
    -f "$results_file" \
    "$input_dir" \
    "$work_dir"

# ----------------------------------------
# Confirm result
# ----------------------------------------

if [ ! -s "$results_file" ]; then
    echo "ERROR: CheckM results table was not produced."
    exit 1
fi

echo
echo "========================================"
echo "CheckM completed successfully"
echo "Sample:  ${sample}"
echo "Bins:    ${bin_count}"
echo "Results: ${results_file}"
echo "Finished: $(date)"
echo "========================================"