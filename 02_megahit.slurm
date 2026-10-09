#!/bin/bash
#SBATCH --job-name=megahit
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --partition=cpu_shared
#SBATCH --time=24:00:00
#SBATCH --mem=64G
#SBATCH --array=0-3
#SBATCH --output=megahit_%A_%a.out
#SBATCH --error=megahit_%A_%a.err

set -eo pipefail

cd "$SLURM_SUBMIT_DIR"

module load conda/24.9.2
source /scratch/conda/miniconda3/20240912/etc/profile.d/conda.sh
conda activate mags

samples=(
    "GB_Biofilm"
    "GB_Sediment"
    "HS_Biofilm"
    "HS_Sediment"
)

sample="${samples[$SLURM_ARRAY_TASK_ID]}"

READS_DIR="00.CleanData"
OUTPUT_ROOT="02.MEGAHIT"

input_dir="${READS_DIR}/${sample}"
output_dir="${OUTPUT_ROOT}/${sample}"

echo "Sample: $sample"
echo "Input:  $input_dir"
echo "Output: $output_dir"

if [ ! -d "$input_dir" ]; then
    echo "ERROR: Input directory does not exist: $input_dir"
    exit 1
fi

# Match files ending specifically in _1.clean.fq.gz or _1.fastq.gz.
mapfile -t r1_files < <(
    find "$input_dir" -maxdepth 1 -type f \
        \( -name "*_1.clean.fq.gz" -o -name "*_1.fastq.gz" \) \
        | sort
)

# Match files ending specifically in _2.clean.fq.gz or _2.fastq.gz.
mapfile -t r2_files < <(
    find "$input_dir" -maxdepth 1 -type f \
        \( -name "*_2.clean.fq.gz" -o -name "*_2.fastq.gz" \) \
        | sort
)

echo "R1 files found: ${#r1_files[@]}"
printf "  %s\n" "${r1_files[@]}"

echo "R2 files found: ${#r2_files[@]}"
printf "  %s\n" "${r2_files[@]}"

if [ "${#r1_files[@]}" -eq 0 ] || \
   [ "${#r2_files[@]}" -eq 0 ]; then
    echo "ERROR: No paired reads found for $sample"
    exit 1
fi

if [ "${#r1_files[@]}" -ne "${#r2_files[@]}" ]; then
    echo "ERROR: R1 and R2 file counts do not match."
    exit 1
fi

# Convert multiple lanes into comma-separated MEGAHIT input lists.
r1_csv=$(IFS=,; echo "${r1_files[*]}")
r2_csv=$(IFS=,; echo "${r2_files[*]}")

mkdir -p "$OUTPUT_ROOT"

if [ -e "$output_dir" ]; then
    echo "ERROR: Output directory already exists:"
    echo "$output_dir"
    echo "Rename or remove the existing directory before rerunning."
    exit 1
fi

megahit \
    -1 "$r1_csv" \
    -2 "$r2_csv" \
    -o "$output_dir" \
    --min-contig-len 1000 \
    -t "$SLURM_CPUS_PER_TASK"

if [ ! -s "${output_dir}/final.contigs.fa" ]; then
    echo "ERROR: final.contigs.fa was not produced."
    exit 1
fi

echo "MEGAHIT completed successfully for $sample"
ls -lh "${output_dir}/final.contigs.fa"