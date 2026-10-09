#!/bin/bash
#SBATCH --job-name=dastool_allraw
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --partition=cpu_shared
#SBATCH --time=24:00:00
#SBATCH --mem=64G
#SBATCH --array=0-3%4
#SBATCH --output=dastool_allraw_%A_%a.out
#SBATCH --error=dastool_allraw_%A_%a.err

set -eo pipefail

PROJECT="$HOME/MAGs_V2"
OUTPUT_ROOT="${PROJECT}/08.DASTool_AllRaw"

cd "$PROJECT"

# ============================================================
# Load DAS Tool environment
# ============================================================

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

assembly="${PROJECT}/02.MEGAHIT/${sample}/final.contigs.fa"

output_dir="${OUTPUT_ROOT}/${sample}"
contig2bin_dir="${output_dir}/contigs2bin"

echo "===================================================="
echo "DAS Tool - ALL RAW BINS"
echo "Sample:   ${sample}"
echo "Assembly: ${assembly}"
echo "Started:  $(date)"
echo "===================================================="

# ============================================================
# Check assembly
# ============================================================

if [ ! -s "$assembly" ]; then
    echo "ERROR: MEGAHIT assembly missing:"
    echo "$assembly"
    exit 1
fi

# ============================================================
# MetaBAT2 bins
# ============================================================

metabat_dir="${PROJECT}/05.Binning/MetaBAT2/${sample}"

if [ ! -d "$metabat_dir" ]; then
    echo "ERROR: MetaBAT2 directory not found:"
    echo "$metabat_dir"
    exit 1
fi

# ============================================================
# SemiBin2 bins
# ============================================================

semibin_root="${PROJECT}/05.Binning/SemiBin2/${sample}"

if [ -d "${semibin_root}/output_recluster_bins" ]; then

    semibin_dir="${semibin_root}/output_recluster_bins"

elif [ -d "${semibin_root}/output_bins" ]; then

    semibin_dir="${semibin_root}/output_bins"

else

    echo "ERROR: SemiBin2 bins not found for ${sample}"
    exit 1

fi

# ============================================================
# COMEBin bins
# ============================================================

comebin_root="${PROJECT}/05.Binning/COMEBin/${sample}"

comebin_candidates=(
    "${comebin_root}/comebin_res/comebin_res_bins"
    "${comebin_root}/comebin_res_bins"
    "${comebin_root}/bins"
)

comebin_dir=""

for candidate in "${comebin_candidates[@]}"; do

    if [ -d "$candidate" ]; then
        comebin_dir="$candidate"
        break
    fi

done

if [ -z "$comebin_dir" ]; then

    echo "ERROR: COMEBin bins not found for ${sample}"
    echo
    echo "Directories searched:"
    printf '  %s\n' "${comebin_candidates[@]}"

    exit 1
fi

# ============================================================
# VAMB bins
# ============================================================

vamb_dir="${PROJECT}/05.Binning/VAMB/${sample}/bins"

if [ ! -d "$vamb_dir" ]; then
    echo "ERROR: VAMB bins not found:"
    echo "$vamb_dir"
    exit 1
fi

# ============================================================
# Detect FASTA extension
# ============================================================

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

# ============================================================
# Count raw bins
# ============================================================

metabat_count=$(find "$metabat_dir" -maxdepth 1 \
    -type f -name "*.${metabat_ext}" | wc -l)

semibin_count=$(find "$semibin_dir" -maxdepth 1 \
    -type f -name "*.${semibin_ext}" | wc -l)

comebin_count=$(find "$comebin_dir" -maxdepth 1 \
    -type f -name "*.${comebin_ext}" | wc -l)

vamb_count=$(find "$vamb_dir" -maxdepth 1 \
    -type f -name "*.${vamb_ext}" | wc -l)

echo
echo "RAW BIN SETS:"
echo "----------------------------------------------------"
echo "MetaBAT2: ${metabat_count} bins"
echo "SemiBin2: ${semibin_count} bins"
echo "COMEBin:  ${comebin_count} bins"
echo "VAMB:     ${vamb_count} bins"
echo "----------------------------------------------------"

total_raw=$((metabat_count + semibin_count + comebin_count + vamb_count))

echo "Total raw bin predictions: ${total_raw}"
echo

# ============================================================
# Prevent accidental overwrite
# ============================================================

if [ -d "$output_dir" ] &&
   [ -n "$(ls -A "$output_dir" 2>/dev/null)" ]; then

    echo "ERROR: DAS Tool output already exists:"
    echo "$output_dir"
    echo
    echo "Remove it before rerunning this sample."
    exit 1

fi

mkdir -p "$contig2bin_dir"

# ============================================================
# Convert FASTA bins into DAS Tool contig-to-bin tables
#
# Output format:
# contig_ID    bin_ID
# ============================================================

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

echo "Creating contig-to-bin tables..."

fasta_to_contigs2bin \
    "$metabat_dir" \
    "$metabat_ext" \
    "${contig2bin_dir}/MetaBAT2_contigs2bin.tsv"

fasta_to_contigs2bin \
    "$semibin_dir" \
    "$semibin_ext" \
    "${contig2bin_dir}/SemiBin2_contigs2bin.tsv"

fasta_to_contigs2bin \
    "$comebin_dir" \
    "$comebin_ext" \
    "${contig2bin_dir}/COMEBin_contigs2bin.tsv"

fasta_to_contigs2bin \
    "$vamb_dir" \
    "$vamb_ext" \
    "${contig2bin_dir}/VAMB_contigs2bin.tsv"

# ============================================================
# Confirm tables contain data
# ============================================================

for table in \
    "${contig2bin_dir}/MetaBAT2_contigs2bin.tsv" \
    "${contig2bin_dir}/SemiBin2_contigs2bin.tsv" \
    "${contig2bin_dir}/COMEBin_contigs2bin.tsv" \
    "${contig2bin_dir}/VAMB_contigs2bin.tsv"
do

    if [ ! -s "$table" ]; then

        echo "ERROR: Empty contig-to-bin table:"
        echo "$table"
        exit 1

    fi

done

echo
echo "Contigs represented in each binning prediction:"

wc -l "${contig2bin_dir}"/*_contigs2bin.tsv

# ============================================================
# Run DAS Tool
# ============================================================

prefix="${output_dir}/${sample}"

echo
echo "===================================================="
echo "Running DAS Tool"
echo "===================================================="
echo

DAS_Tool \
    -i "${contig2bin_dir}/MetaBAT2_contigs2bin.tsv,${contig2bin_dir}/SemiBin2_contigs2bin.tsv,${contig2bin_dir}/COMEBin_contigs2bin.tsv,${contig2bin_dir}/VAMB_contigs2bin.tsv" \
    -l "MetaBAT2,SemiBin2,COMEBin,VAMB" \
    -c "$assembly" \
    -o "$prefix" \
    --search_engine diamond \
    --write_bin_evals \
    --write_bins \
    -t "$SLURM_CPUS_PER_TASK"

# ============================================================
# Check DAS Tool output
# ============================================================

summary="${prefix}_DASTool_summary.tsv"
bins_dir="${prefix}_DASTool_bins"

if [ ! -s "$summary" ]; then

    echo
    echo "ERROR: DAS Tool summary was not generated:"
    echo "$summary"
    exit 1

fi

if [ ! -d "$bins_dir" ]; then

    echo
    echo "ERROR: DAS Tool bins directory was not generated:"
    echo "$bins_dir"
    exit 1

fi

final_bin_count=$(find "$bins_dir" \
    -maxdepth 1 \
    -type f \
    \( -name "*.fa" -o -name "*.fasta" -o -name "*.fna" \) \
    | wc -l)

echo
echo "===================================================="
echo "DAS Tool completed successfully"
echo "===================================================="
echo
echo "Sample: ${sample}"
echo
echo "Raw candidate bins entering DAS Tool:"
echo "  MetaBAT2 = ${metabat_count}"
echo "  SemiBin2 = ${semibin_count}"
echo "  COMEBin  = ${comebin_count}"
echo "  VAMB     = ${vamb_count}"
echo
echo "Total raw predictions = ${total_raw}"
echo
echo "DAS Tool selected bins = ${final_bin_count}"
echo
echo "Summary:"
echo "$summary"
echo
echo "Final bins:"
echo "$bins_dir"
echo
echo "Finished: $(date)"
echo "===================================================="