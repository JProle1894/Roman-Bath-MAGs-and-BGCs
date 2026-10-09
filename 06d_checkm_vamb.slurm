#!/bin/bash
#SBATCH --job-name=checkm_vamb
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --partition=cpu_shared
#SBATCH --time=24:00:00
#SBATCH --mem=64G
#SBATCH --array=0-3%1
#SBATCH --output=checkm_vamb_%A_%a.out
#SBATCH --error=checkm_vamb_%A_%a.err

set -eo pipefail

cd "$SLURM_SUBMIT_DIR"

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

input_dir="05.Binning/VAMB/${sample}/bins"
output_dir="06.CheckM/VAMB/${sample}"
work_dir="${output_dir}/work"
results_file="${output_dir}/checkm_results.tsv"

echo "========================================"
echo "CheckM assessment of VAMB bins"
echo "Sample:     ${sample}"
echo "Input:      ${input_dir}"
echo "Output:     ${output_dir}"
echo "Started:    $(date)"
echo "========================================"

if [ ! -d "$input_dir" ]; then
    echo "ERROR: VAMB bins directory not found:"
    echo "$input_dir"
    exit 1
fi

# VAMB normally writes bins with the .fna extension.
bin_count=$(find "$input_dir" \
    -maxdepth 1 \
    -type f \
    -name "*.fna" \
    | wc -l)

if [ "$bin_count" -eq 0 ]; then
    echo "ERROR: No .fna VAMB bins found in:"
    echo "$input_dir"
    exit 1
fi

echo "VAMB bins found: ${bin_count}"

if [ -d "$output_dir" ] &&
   [ -n "$(ls -A "$output_dir" 2>/dev/null)" ]; then

    echo "ERROR: CheckM output directory already contains files:"
    echo "$output_dir"
    echo "Remove or rename it before rerunning."
    exit 1
fi

mkdir -p "$output_dir"

checkm lineage_wf \
    -x fna \
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