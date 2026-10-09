#!/bin/bash
#SBATCH --job-name=gunc_das_allraw
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --partition=cpu_shared
#SBATCH --time=24:00:00
#SBATCH --mem=64G
#SBATCH --array=0-3%4
#SBATCH --output=gunc_das_allraw_%A_%a.out
#SBATCH --error=gunc_das_allraw_%A_%a.err

set -eo pipefail

PROJECT="$HOME/MAGs_V2"

INPUT_ROOT="${PROJECT}/08.DASTool_AllRaw"
OUTPUT_ROOT="${PROJECT}/11.GUNC_DASTool_AllRaw"

GUNC_DB_ROOT="$HOME/databases/GUNC"

cd "$PROJECT"

# ============================================================
# Load Conda / GUNC
# ============================================================

module load conda/24.9.2
source /scratch/conda/miniconda3/20240912/etc/profile.d/conda.sh

conda activate gunc_env

command -v gunc >/dev/null 2>&1 || {
    echo "ERROR: GUNC not found in gunc_env"
    exit 1
}

# ============================================================
# Samples
# ============================================================

samples=(
    "GB_Biofilm"
    "GB_Sediment"
    "HS_Biofilm"
    "HS_Sediment"
)

sample="${samples[$SLURM_ARRAY_TASK_ID]}"

# ============================================================
# Input / output
# ============================================================

input_dir="${INPUT_ROOT}/${sample}/${sample}_DASTool_bins"

output_dir="${OUTPUT_ROOT}/${sample}"
temp_dir="${output_dir}/tmp"

echo "================================================="
echo "GUNC on DAS Tool ALL-RAW output"
echo "Sample:  ${sample}"
echo "Input:   ${input_dir}"
echo "Output:  ${output_dir}"
echo "CPUs:    ${SLURM_CPUS_PER_TASK}"
echo "Started: $(date)"
echo "================================================="

# ============================================================
# Check DAS Tool bins
# ============================================================

if [ ! -d "$input_dir" ]; then
    echo "ERROR: DAS Tool bins directory not found:"
    echo "$input_dir"
    exit 1
fi

# ============================================================
# Detect FASTA extension
# ============================================================

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
    echo "ERROR: No FASTA bins found in:"
    echo "$input_dir"
    exit 1
fi

bin_count=$(find "$input_dir" \
    -maxdepth 1 \
    -type f \
    -name "*.${bin_extension}" \
    | wc -l)

echo
echo "FASTA extension: .${bin_extension}"
echo "Number of DAS Tool bins: ${bin_count}"
echo

# ============================================================
# Find GUNC database
# ============================================================

if [ ! -d "$GUNC_DB_ROOT" ]; then
    echo "ERROR: GUNC database directory not found:"
    echo "$GUNC_DB_ROOT"
    exit 1
fi

mapfile -t GUNC_DATABASES < <(
    find "$GUNC_DB_ROOT" \
        -type f \
        -name "*.dmnd" \
        ! -path "*ci_test*" \
        ! -path "*test_data*" \
        | sort
)

if [ "${#GUNC_DATABASES[@]}" -eq 0 ]; then

    echo "ERROR: No GUNC .dmnd database found under:"
    echo "$GUNC_DB_ROOT"
    exit 1

elif [ "${#GUNC_DATABASES[@]}" -gt 1 ]; then

    echo "ERROR: More than one GUNC database was found."
    echo
    printf '  %s\n' "${GUNC_DATABASES[@]}"
    echo
    echo "Specify which database should be used in the script."
    exit 1

fi

GUNC_DB="${GUNC_DATABASES[0]}"

echo "GUNC database:"
echo "$GUNC_DB"
echo

# ============================================================
# Prevent accidental overwrite
# ============================================================

if [ -d "$output_dir" ] &&
   find "$output_dir" \
       -maxdepth 1 \
       -type f \
       -name "GUNC.*.maxCSS_level.tsv" \
       | grep -q .; then

    echo "ERROR: GUNC results already exist for ${sample}:"
    echo "$output_dir"
    echo
    echo "Remove or rename them before rerunning."
    exit 1

fi

mkdir -p "$output_dir"
mkdir -p "$temp_dir"

# ============================================================
# Check GUNC setup
# ============================================================

echo "Checking GUNC environment/database..."
echo

gunc check \
    --db_file "$GUNC_DB" \
    --out_dir "$output_dir"

# ============================================================
# Run GUNC
# ============================================================

echo
echo "Running GUNC..."
echo

gunc run \
    --input_dir "$input_dir" \
    --file_suffix ".${bin_extension}" \
    --db_file "$GUNC_DB" \
    --out_dir "$output_dir" \
    --temp_dir "$temp_dir" \
    --threads "$SLURM_CPUS_PER_TASK"

# ============================================================
# Find final GUNC results
# ============================================================

mapfile -t result_files < <(
    find "$output_dir" \
        -maxdepth 1 \
        -type f \
        -name "GUNC.*.maxCSS_level.tsv"
)

if [ "${#result_files[@]}" -eq 0 ]; then

    echo
    echo "ERROR: No GUNC maxCSS result file was generated."
    exit 1

fi

echo
echo "================================================="
echo "GUNC completed successfully"
echo "================================================="
echo
echo "Sample: ${sample}"
echo "Bins analysed: ${bin_count}"
echo "Finished: $(date)"
echo
echo "Main GUNC result file(s):"

printf '%s\n' "${result_files[@]}"

echo
echo "First few rows:"
head "${result_files[0]}"

echo
echo "================================================="