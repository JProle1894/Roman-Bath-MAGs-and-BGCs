#!/bin/bash
#SBATCH --job-name=dastool
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --partition=cpu_shared
#SBATCH --time=24:00:00
#SBATCH --mem=64G
#SBATCH --array=0-3%4
#SBATCH --output=dastool_%A_%a.out
#SBATCH --error=dastool_%A_%a.err

set -eo pipefail

PROJECT="$HOME/MAGs_V2"
cd "$PROJECT"

module load conda/24.9.2
source /scratch/conda/miniconda3/20240912/etc/profile.d/conda.sh
conda activate dastool_env

command -v DAS_Tool >/dev/null 2>&1 || {
    echo "ERROR: DAS_Tool not found in dastool_env"
    exit 1
}

command -v diamond >/dev/null 2>&1 || {
    echo "ERROR: DIAMOND not found in dastool_env"
    exit 1
}

samples=(
    "GB_Biofilm"
    "GB_Sediment"
    "HS_Biofilm"
    "HS_Sediment"
)

sample="${samples[$SLURM_ARRAY_TASK_ID]}"

assembly="02.MEGAHIT/${sample}/final.contigs.fa"

output_dir="08.DASTool/${sample}"
input_dir="${output_dir}/contigs2bin"

echo "========================================"
echo "DAS Tool refinement"
echo "Sample:   ${sample}"
echo "Assembly: ${assembly}"
echo "Started:  $(date)"
echo "========================================"

# ----------------------------------------
# Check assembly
# ----------------------------------------

if [ ! -s "$assembly" ]; then
    echo "ERROR: Assembly missing or empty:"
    echo "$assembly"
    exit 1
fi

# ----------------------------------------
# MetaBAT2 bins
# ----------------------------------------

metabat_dir="05.Binning/MetaBAT2/${sample}"

if [ ! -d "$metabat_dir" ]; then
    echo "ERROR: MetaBAT2 directory missing:"
    echo "$metabat_dir"
    exit 1
fi

# ----------------------------------------
# SemiBin2 bins
# ----------------------------------------

semibin_root="05.Binning/SemiBin2/${sample}"

if [ -d "${semibin_root}/output_recluster_bins" ]; then

    semibin_dir="${semibin_root}/output_recluster_bins"

elif [ -d "${semibin_root}/output_bins" ]; then

    semibin_dir="${semibin_root}/output_bins"

else

    echo "ERROR: Could not find SemiBin2 bins for ${sample}"
    exit 1

fi

# ----------------------------------------
# COMEBin bins
# ----------------------------------------

comebin_root="05.Binning/COMEBin/${sample}"

comebin_dir=""

comebin_candidates=(
    "${comebin_root}/comebin_res/comebin_res_bins"
    "${comebin_root}/comebin_res_bins"
    "${comebin_root}/bins"
)

for candidate in "${comebin_candidates[@]}"; do

    if [ -d "$candidate" ]; then
        comebin_dir="$candidate"
        break
    fi

done

if [ -z "$comebin_dir" ]; then

    echo "ERROR: Could not find COMEBin bins for ${sample}"

    echo "Directories under:"
    echo "$comebin_root"

    find "$comebin_root" -maxdepth 3 -type d 2>/dev/null || true

    exit 1
fi

# ----------------------------------------
# VAMB bins
# ----------------------------------------

vamb_dir="05.Binning/VAMB/${sample}/bins"

if [ ! -d "$vamb_dir" ]; then
    echo "ERROR: VAMB bins directory missing:"
    echo "$vamb_dir"
    exit 1
fi

# ----------------------------------------
# Detect FASTA extension
# ----------------------------------------

detect_extension () {

    local directory="$1"

    for ext in fa fasta fna fas; do

        if find "$directory" \
            -maxdepth 1 \
            -type f \
            -name "*.${ext}" \
            -print -quit \
            | grep -q .; then

            echo "$ext"
            return 0

        fi

    done

    return 1
}

metabat_ext=$(detect_extension "$metabat_dir") || {
    echo "ERROR: No MetaBAT2 FASTA bins found"
    exit 1
}

semibin_ext=$(detect_extension "$semibin_dir") || {
    echo "ERROR: No SemiBin2 FASTA bins found"
    exit 1
}

comebin_ext=$(detect_extension "$comebin_dir") || {
    echo "ERROR: No COMEBin FASTA bins found"
    exit 1
}

vamb_ext=$(detect_extension "$vamb_dir") || {
    echo "ERROR: No VAMB FASTA bins found"
    exit 1
}

echo
echo "Bin directories detected:"
echo "MetaBAT2: ${metabat_dir} (.${metabat_ext})"
echo "SemiBin2: ${semibin_dir} (.${semibin_ext})"
echo "COMEBin:  ${comebin_dir} (.${comebin_ext})"
echo "VAMB:     ${vamb_dir} (.${vamb_ext})"

# ----------------------------------------
# Stop if previous DAS Tool output exists
# ----------------------------------------

if [ -d "$output_dir" ] &&
   [ -n "$(ls -A "$output_dir" 2>/dev/null)" ]; then

    echo
    echo "ERROR: Output directory already contains files:"
    echo "$output_dir"
    echo
    echo "Remove the old/failed output before rerunning."
    exit 1
fi

mkdir -p "$input_dir"

# ----------------------------------------
# Convert FASTA bins to contig2bin tables
# ----------------------------------------

fasta_to_contigs2bin () {

    local bin_dir="$1"
    local extension="$2"
    local output="$3"

    > "$output"

    for fasta in "${bin_dir}"/*."${extension}"; do

        [ -f "$fasta" ] || continue

        bin_id=$(basename "$fasta" ".${extension}")

        awk -v bin="$bin_id" '
            /^>/ {
                contig=$1
                sub(/^>/, "", contig)
                print contig "\t" bin
            }
        ' "$fasta" >> "$output"

    done
}

echo
echo "Creating contig-to-bin tables..."

fasta_to_contigs2bin \
    "$metabat_dir" \
    "$metabat_ext" \
    "${input_dir}/MetaBAT2_contigs2bin.tsv"

fasta_to_contigs2bin \
    "$semibin_dir" \
    "$semibin_ext" \
    "${input_dir}/SemiBin2_contigs2bin.tsv"

fasta_to_contigs2bin \
    "$comebin_dir" \
    "$comebin_ext" \
    "${input_dir}/COMEBin_contigs2bin.tsv"

fasta_to_contigs2bin \
    "$vamb_dir" \
    "$vamb_ext" \
    "${input_dir}/VAMB_contigs2bin.tsv"

# ----------------------------------------
# Check tables
# ----------------------------------------

for table in \
    "${input_dir}/MetaBAT2_contigs2bin.tsv" \
    "${input_dir}/SemiBin2_contigs2bin.tsv" \
    "${input_dir}/COMEBin_contigs2bin.tsv" \
    "${input_dir}/VAMB_contigs2bin.tsv"
do

    if [ ! -s "$table" ]; then

        echo "ERROR: Empty contig2bin table:"
        echo "$table"
        exit 1

    fi

done

echo
echo "Contig-to-bin tables created:"
wc -l "${input_dir}"/*_contigs2bin.tsv

# ----------------------------------------
# Run DAS Tool
# ----------------------------------------

prefix="${output_dir}/${sample}"

echo
echo "Running DAS Tool..."
echo

DAS_Tool \
    -i "${input_dir}/MetaBAT2_contigs2bin.tsv,${input_dir}/SemiBin2_contigs2bin.tsv,${input_dir}/COMEBin_contigs2bin.tsv,${input_dir}/VAMB_contigs2bin.tsv" \
    -l "MetaBAT2,SemiBin2,COMEBin,VAMB" \
    -c "$assembly" \
    -o "$prefix" \
    --search_engine diamond \
    --write_bin_evals \
    --write_bins \
    -t "$SLURM_CPUS_PER_TASK"

# ----------------------------------------
# Check final output
# ----------------------------------------

summary="${prefix}_DASTool_summary.tsv"

if [ ! -s "$summary" ]; then

    echo
    echo "ERROR: DAS Tool summary was not produced:"
    echo "$summary"
    exit 1

fi

bins_dir="${prefix}_DASTool_bins"

echo
echo "========================================"
echo "DAS Tool completed successfully"
echo "Sample: ${sample}"
echo "Finished: $(date)"
echo
echo "Summary:"
echo "$summary"
echo
echo "Refined bins:"
echo "$bins_dir"

if [ -d "$bins_dir" ]; then

    bin_count=$(find "$bins_dir" \
        -maxdepth 1 \
        -type f \
        \( -name "*.fa" -o -name "*.fasta" -o -name "*.fna" \) \
        | wc -l)

    echo "Number of DAS Tool bins: ${bin_count}"

else

    echo "WARNING: DAS Tool bins directory was not found."

fi

echo "========================================"