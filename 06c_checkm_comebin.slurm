#!/bin/bash
#SBATCH --job-name=checkm_comebin
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --partition=cpu_shared
#SBATCH --time=24:00:00
#SBATCH --mem=64G
#SBATCH --array=0-3%1
#SBATCH --output=checkm_comebin_%A_%a.out
#SBATCH --error=checkm_comebin_%A_%a.err

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

comebin_root="05.Binning/COMEBin/${sample}"

# Find the COMEBin final bins directory.
candidate_dirs=(
    "${comebin_root}/comebin_res/comebin_res_bins"
    "${comebin_root}/comebin_res_bins"
    "${comebin_root}/bins"
)

input_dir=""

for candidate in "${candidate_dirs[@]}"; do
    if [ -d "$candidate" ]; then
        input_dir="$candidate"
        break
    fi
done

if [ -z "$input_dir" ]; then
    echo "ERROR: Could not find the COMEBin bins directory."
    echo
    echo "Directories checked:"
    printf '  %s\n' "${candidate_dirs[@]}"
    echo
    echo "Existing directories under ${comebin_root}:"
    find "$comebin_root" -maxdepth 3 -type d 2>/dev/null || true
    exit 1
fi

# Determine which FASTA extension COMEBin used.
bin_extension=""

for extension in fa fasta fna; do
    count=$(find "$input_dir" \
        -maxdepth 1 \
        -type f \
        -name "*.${extension}" \
        | wc -l)

    if [ "$count" -gt 0 ]; then
        bin_extension="$extension"
        bin_count="$count"
        break
    fi
done

if [ -z "$bin_extension" ]; then
    echo "ERROR: No .fa, .fasta or .fna bins were found in:"
    echo "$input_dir"
    exit 1
fi

output_dir="06.CheckM/COMEBin/${sample}"
work_dir="${output_dir}/work"
results_file="${output_dir}/checkm_results.tsv"

echo "========================================"
echo "CheckM assessment of COMEBin bins"
echo "Sample:       ${sample}"
echo "Input:        ${input_dir}"
echo "Extension:    .${bin_extension}"
echo "Bins found:   ${bin_count}"
echo "Output:       ${output_dir}"
echo "Started:      $(date)"
echo "========================================"

if [ -d "$output_dir" ] &&
   [ -n "$(ls -A "$output_dir" 2>/dev/null)" ]; then

    echo "ERROR: Output directory already contains files:"
    echo "$output_dir"
    echo
    echo "Remove it before rerunning this sample."
    exit 1
fi

mkdir -p "$output_dir"

checkm lineage_wf \
    -x "$bin_extension" \
    -t "$SLURM_CPUS_PER_TASK" \
    --tab_table \
    -f "$results_file" \
    "$input_dir" \
    "$work_dir"

if [ ! -s "$results_file" ]; then
    echo "ERROR: CheckM results table was not produced."
    exit 1
fi

echo "========================================"
echo "CheckM completed successfully"
echo "Sample:  ${sample}"
echo "Results: ${results_file}"
echo "Finished: $(date)"
echo "========================================"