#!/bin/bash
#SBATCH --job-name=checkm_sb2
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --partition=cpu_shared
#SBATCH --time=24:00:00
#SBATCH --mem=64G
#SBATCH --array=0-3%1
#SBATCH --output=checkm_semibin2_%A_%a.out
#SBATCH --error=checkm_semibin2_%A_%a.err

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

base_dir="05.Binning/SemiBin2/${sample}"

if [ -d "${base_dir}/output_recluster_bins" ]; then
    input_dir="${base_dir}/output_recluster_bins"
elif [ -d "${base_dir}/output_bins" ]; then
    input_dir="${base_dir}/output_bins"
else
    echo "ERROR: SemiBin2 bin directory not found in ${base_dir}"
    exit 1
fi

output_dir="06.CheckM/SemiBin2/${sample}"

echo "Sample:     ${sample}"
echo "Input bins: ${input_dir}"
echo "Output:     ${output_dir}"

bin_count=$(find "$input_dir" -maxdepth 1 -type f -name "*.fa" | wc -l)

if [ "$bin_count" -eq 0 ]; then
    echo "ERROR: No .fa bins found in $input_dir"
    exit 1
fi

echo "Bins found: ${bin_count}"

if [ -e "$output_dir" ]; then
    echo "ERROR: Output directory already exists: $output_dir"
    echo "Remove it before rerunning this sample."
    exit 1
fi

mkdir -p "$(dirname "$output_dir")"

checkm lineage_wf \
    -x fa \
    -t "$SLURM_CPUS_PER_TASK" \
    --tab_table \
    -f "${output_dir}_results.tsv" \
    "$input_dir" \
    "$output_dir"

if [ ! -s "${output_dir}_results.tsv" ]; then
    echo "ERROR: CheckM results table was not produced."
    exit 1
fi

echo "CheckM completed successfully for ${sample}"
echo "Results: ${output_dir}_results.tsv"