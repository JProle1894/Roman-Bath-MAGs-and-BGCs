#!/bin/bash
#SBATCH --job-name=gunc_dastool
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --partition=cpu_shared
#SBATCH --time=24:00:00
#SBATCH --mem=64G
#SBATCH --array=0-3%4
#SBATCH --output=gunc_dastool_%A_%a.out
#SBATCH --error=gunc_dastool_%A_%a.err

set -eo pipefail

PROJECT="$HOME/MAGs_V2"
GUNC_DB_DIR="$HOME/databases/GUNC"

cd "$PROJECT"

module load conda/24.9.2
source /scratch/conda/miniconda3/20240912/etc/profile.d/conda.sh
conda activate gunc_env

command -v gunc >/dev/null 2>&1 || {
    echo "ERROR: GUNC is not available in gunc_env"
    exit 1
}

samples=(
    "GB_Biofilm"
    "GB_Sediment"
    "HS_Biofilm"
    "HS_Sediment"
)

sample="${samples[$SLURM_ARRAY_TASK_ID]}"

# ----------------------------------------
# Find GUNC database
# ----------------------------------------

mapfile -t gunc_databases < <(
    find "$GUNC_DB_DIR" \
        -type f \
        -name "*.dmnd" \
        ! -name "ci_test.dmnd" \
        2>/dev/null
)

if [ "${#gunc_databases[@]}" -eq 0 ]; then
    echo "ERROR: No GUNC database found under:"
    echo "$GUNC_DB_DIR"
    exit 1
fi

if [ "${#gunc_databases[@]}" -gt 1 ]; then
    echo "ERROR: More than one GUNC database found:"
    printf '  %s\n' "${gunc_databases[@]}"
    echo "Specify which database should be used."
    exit 1
fi

GUNC_DB="${gunc_databases[0]}"

# ----------------------------------------
# Input/output paths
# ----------------------------------------

input_dir="08.DASTool/${sample}/${sample}_DASTool_bins"
output_dir="11.GUNC_DASTool/${sample}"
temp_dir="${output_dir}/tmp"

echo "========================================"
echo "GUNC assessment of DAS Tool bins"
echo "Sample:   ${sample}"
echo "Input:    ${input_dir}"
echo "Database: ${GUNC_DB}"
echo "Output:   ${output_dir}"
echo "Started:  $(date)"
echo "========================================"

# ----------------------------------------
# Check input
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

    echo "ERROR: GUNC output directory already contains files:"
    echo "$output_dir"
    exit 1

fi

mkdir -p "$output_dir"
mkdir -p "$temp_dir"

# ----------------------------------------
# Validate GUNC setup
# ----------------------------------------

gunc check \
    --db_file "$GUNC_DB" \
    --out_dir "$output_dir"

# ----------------------------------------
# Run GUNC
# ----------------------------------------

gunc run \
    --input_dir "$input_dir" \
    --file_suffix ".${bin_extension}" \
    --db_file "$GUNC_DB" \
    --out_dir "$output_dir" \
    --temp_dir "$temp_dir" \
    --threads "$SLURM_CPUS_PER_TASK"

# ----------------------------------------
# Check result
# ----------------------------------------

result_file=$(find "$output_dir" \
    -maxdepth 1 \
    -type f \
    -name "GUNC.*.maxCSS_level.tsv" \
    -print -quit)

if [ -z "$result_file" ] || [ ! -s "$result_file" ]; then
    echo "ERROR: GUNC maxCSS results table was not produced."
    exit 1
fi

echo
echo "========================================"
echo "GUNC completed successfully"
echo "Sample:  ${sample}"
echo "Bins:    ${bin_count}"
echo "Results: ${result_file}"
echo "Finished: $(date)"
echo "========================================"

echo
echo "GUNC pass/fail summary:"

awk -F'\t' '
NR==1 {
    for(i=1;i<=NF;i++)
        if($i=="pass.GUNC")
            col=i
    next
}
col {
    count[$col]++
}
END {
    for(x in count)
        print x, count[x]
}
' "$result_file"