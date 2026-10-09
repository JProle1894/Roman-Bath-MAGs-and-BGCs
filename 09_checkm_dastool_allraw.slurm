#!/bin/bash
#SBATCH --job-name=checkm_das_allraw
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --partition=cpu_shared
#SBATCH --time=24:00:00
#SBATCH --mem=64G
#SBATCH --array=0-3%4
#SBATCH --output=checkm_das_allraw_%A_%a.out
#SBATCH --error=checkm_das_allraw_%A_%a.err

set -eo pipefail

PROJECT="$HOME/MAGs_V2"
INPUT_ROOT="${PROJECT}/08.DASTool_AllRaw"
OUTPUT_ROOT="${PROJECT}/09.CheckM_DASTool_AllRaw"

cd "$PROJECT"

# ------------------------------------------------
# Load Conda / CheckM
# ------------------------------------------------

module load conda/24.9.2
source /scratch/conda/miniconda3/20240912/etc/profile.d/conda.sh

conda activate checkm_env

command -v checkm >/dev/null 2>&1 || {
    echo "ERROR: CheckM not found in checkm_env"
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

input_dir="${INPUT_ROOT}/${sample}/${sample}_DASTool_bins"

output_dir="${OUTPUT_ROOT}/${sample}"
work_dir="${output_dir}/checkm_work"
results_file="${output_dir}/checkm_results.tsv"

echo "================================================="
echo "CheckM on DAS Tool ALL-RAW output"
echo "Sample: ${sample}"
echo "Input:  ${input_dir}"
echo "Output: ${output_dir}"
echo "CPUs:   ${SLURM_CPUS_PER_TASK}"
echo "Started: $(date)"
echo "================================================="

# ------------------------------------------------
# Check input exists
# ------------------------------------------------

if [ ! -d "$input_dir" ]; then
    echo "ERROR: DAS Tool bins directory not found:"
    echo "$input_dir"
    exit 1
fi

# ------------------------------------------------
# Detect FASTA extension
# ------------------------------------------------

bin_extension=""

for ext in fa fasta fna fas; do

    if find "$input_dir" \
        -maxdepth 1 \
        -type f \
        -name "*.${ext}" \
        -print -quit \
        | grep -q .; then

        bin_extension="$ext"
        break

    fi
done

if [ -z "$bin_extension" ]; then
    echo "ERROR: No bin FASTA files found in:"
    echo "$input_dir"
    exit 1
fi

bin_count=$(find "$input_dir" \
    -maxdepth 1 \
    -type f \
    -name "*.${bin_extension}" \
    | wc -l)

echo
echo "Detected extension: .${bin_extension}"
echo "Number of DAS Tool bins: ${bin_count}"
echo

# ------------------------------------------------
# Prevent accidental overwrite
# ------------------------------------------------

if [ -s "$results_file" ]; then
    echo "ERROR: CheckM result already exists:"
    echo "$results_file"
    echo
    echo "Remove or rename the output directory before rerunning."
    exit 1
fi

mkdir -p "$output_dir"
mkdir -p "$work_dir"

# ------------------------------------------------
# Run CheckM
# ------------------------------------------------

echo "Running CheckM lineage workflow..."
echo

checkm lineage_wf \
    -x "$bin_extension" \
    -t "$SLURM_CPUS_PER_TASK" \
    --tab_table \
    -f "$results_file" \
    "$input_dir" \
    "$work_dir"

# ------------------------------------------------
# Verify output
# ------------------------------------------------

if [ ! -s "$results_file" ]; then
    echo "ERROR: CheckM result file was not generated:"
    echo "$results_file"
    exit 1
fi

echo
echo "================================================="
echo "CheckM completed successfully"
echo "Sample: ${sample}"
echo "Finished: $(date)"
echo "================================================="
echo
echo "Results:"
echo "$results_file"
echo
echo "First few rows:"
head "$results_file"