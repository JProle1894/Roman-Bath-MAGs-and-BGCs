#!/bin/bash
#SBATCH --job-name=fastqc
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --partition=cpu_shared
#SBATCH --time=02:00:00
#SBATCH --output=fastqc_%j.out
#SBATCH --error=fastqc_%j.err

set -eo pipefail

cd "$SLURM_SUBMIT_DIR"

module load conda/24.9.2
source /scratch/conda/miniconda3/20240912/etc/profile.d/conda.sh
conda activate mags

READS_DIR="00.CleanData"
OUTPUT_DIR="01.FastQC"

mkdir -p "$OUTPUT_DIR"

mapfile -t reads < <(
    find "$READS_DIR" -type f \
        \( -name "*.fq.gz" -o -name "*.fastq.gz" \) \
        | sort
)

if [ "${#reads[@]}" -eq 0 ]; then
    echo "ERROR: No FASTQ files found in $READS_DIR"
    exit 1
fi

echo "FastQC input files: ${#reads[@]}"

fastqc \
    --threads "$SLURM_CPUS_PER_TASK" \
    --outdir "$OUTPUT_DIR" \
    "${reads[@]}"

html_count=$(find "$OUTPUT_DIR" -type f \
    -name "*_fastqc.html" | wc -l)

zip_count=$(find "$OUTPUT_DIR" -type f \
    -name "*_fastqc.zip" | wc -l)

echo "Input files:  ${#reads[@]}"
echo "HTML reports: $html_count"
echo "ZIP reports:  $zip_count"

if [ "$html_count" -ne "${#reads[@]}" ] || \
   [ "$zip_count" -ne "${#reads[@]}" ]; then
    echo "ERROR: The number of FastQC reports does not match the inputs."
    exit 1
fi

echo "FastQC completed successfully."