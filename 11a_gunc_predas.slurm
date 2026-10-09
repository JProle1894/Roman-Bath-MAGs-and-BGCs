#!/bin/bash
#SBATCH --job-name=gunc_predas
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --partition=cpu_shared
#SBATCH --time=24:00:00
#SBATCH --mem=64G
#SBATCH --array=0-15%4
#SBATCH --output=gunc_predas_%A_%a.out
#SBATCH --error=gunc_predas_%A_%a.err

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

# ------------------------------------------------
# Samples and binners
# ------------------------------------------------

samples=(
    "GB_Biofilm"
    "GB_Sediment"
    "HS_Biofilm"
    "HS_Sediment"
)

binners=(
    "MetaBAT2"
    "SemiBin2"
    "COMEBin"
    "VAMB"
)

# 16 total combinations:
# 4 binners x 4 samples

binner_index=$(( SLURM_ARRAY_TASK_ID / 4 ))
sample_index=$(( SLURM_ARRAY_TASK_ID % 4 ))

binner="${binners[$binner_index]}"
sample="${samples[$sample_index]}"

echo "========================================"
echo "Pre-DAS Tool GUNC"
echo "Array task: ${SLURM_ARRAY_TASK_ID}"
echo "Binner:     ${binner}"
echo "Sample:     ${sample}"
echo "Started:    $(date)"
echo "========================================"

# ------------------------------------------------
# Find GUNC database
# ------------------------------------------------

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
    exit 1
fi

GUNC_DB="${gunc_databases[0]}"

echo "GUNC database: ${GUNC_DB}"

# ------------------------------------------------
# Find correct bins directory
# ------------------------------------------------

case "$binner" in

    MetaBAT2)

        input_dir="05.Binning/MetaBAT2/${sample}"
        ;;

    SemiBin2)

        semibin_root="05.Binning/SemiBin2/${sample}"

        if [ -d "${semibin_root}/output_recluster_bins" ]; then
            input_dir="${semibin_root}/output_recluster_bins"

        elif [ -d "${semibin_root}/output_bins" ]; then
            input_dir="${semibin_root}/output_bins"

        else
            echo "ERROR: SemiBin2 bins not found for ${sample}"
            exit 1
        fi
        ;;

    COMEBin)

        comebin_root="05.Binning/COMEBin/${sample}"
        input_dir=""

        candidates=(
            "${comebin_root}/comebin_res/comebin_res_bins"
            "${comebin_root}/comebin_res_bins"
            "${comebin_root}/bins"
        )

        for candidate in "${candidates[@]}"; do

            if [ -d "$candidate" ]; then
                input_dir="$candidate"
                break
            fi

        done

        if [ -z "$input_dir" ]; then
            echo "ERROR: COMEBin bins not found for ${sample}"
            exit 1
        fi
        ;;

    VAMB)

        input_dir="05.Binning/VAMB/${sample}/bins"
        ;;

    *)

        echo "ERROR: Unknown binner: ${binner}"
        exit 1
        ;;

esac

if [ ! -d "$input_dir" ]; then
    echo "ERROR: Input directory not found:"
    echo "$input_dir"
    exit 1
fi

# ------------------------------------------------
# Detect FASTA extension
# ------------------------------------------------

bin_extension=""
bin_count=0

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

echo "Input:     ${input_dir}"
echo "Extension: .${bin_extension}"
echo "Bins:      ${bin_count}"

# ------------------------------------------------
# Output directories
# ------------------------------------------------

output_dir="11.GUNC_PreDAS/${binner}/${sample}"
temp_dir="${output_dir}/tmp"

if [ -d "$output_dir" ] &&
   [ -n "$(ls -A "$output_dir" 2>/dev/null)" ]; then

    echo "ERROR: Output already exists:"
    echo "$output_dir"
    exit 1

fi

mkdir -p "$output_dir"
mkdir -p "$temp_dir"

# ------------------------------------------------
# Check GUNC setup
# ------------------------------------------------

gunc check \
    --db_file "$GUNC_DB" \
    --out_dir "$output_dir"

# ------------------------------------------------
# Run GUNC
# ------------------------------------------------

gunc run \
    --input_dir "$input_dir" \
    --file_suffix ".${bin_extension}" \
    --db_file "$GUNC_DB" \
    --out_dir "$output_dir" \
    --temp_dir "$temp_dir" \
    --threads "$SLURM_CPUS_PER_TASK"

# ------------------------------------------------
# Find results
# ------------------------------------------------

result_file=$(find "$output_dir" \
    -maxdepth 1 \
    -type f \
    -name "GUNC.*.maxCSS_level.tsv" \
    -print -quit)

if [ -z "$result_file" ] || [ ! -s "$result_file" ]; then
    echo "ERROR: GUNC result table was not produced."
    exit 1
fi

result_count=$(tail -n +2 "$result_file" | wc -l)

echo
echo "========================================"
echo "GUNC completed successfully"
echo "Binner:  ${binner}"
echo "Sample:  ${sample}"
echo "Bins in: ${bin_count}"
echo "Bins scored: ${result_count}"
echo "Results: ${result_file}"
echo "Finished: $(date)"
echo "========================================"