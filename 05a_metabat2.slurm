#!/bin/bash
#SBATCH --job-name=metabat2
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --partition=cpu_shared
#SBATCH --time=08:00:00
#SBATCH --mem=32G
#SBATCH --array=0-3%2
#SBATCH --output=logs/metabat2_%A_%a.out
#SBATCH --error=logs/metabat2_%A_%a.err

set -eo pipefail

ROOT="$SLURM_SUBMIT_DIR"
cd "$ROOT"

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

assembly="${ROOT}/02.MEGAHIT/${sample}/final.contigs.fa"
depth="${ROOT}/04.ContigDepth/${sample}/${sample}_depth.txt"
output_dir="${ROOT}/05.Binning/MetaBAT2/${sample}"

echo "Sample:   $sample"
echo "Assembly: $assembly"
echo "Depth:    $depth"
echo "Output:   $output_dir"

command -v metabat2 >/dev/null 2>&1 || {
    echo "ERROR: metabat2 is not available in the active environment."
    exit 1
}

if [ ! -s "$assembly" ]; then
    echo "ERROR: Assembly not found or empty: $assembly"
    exit 1
fi

if [ ! -s "$depth" ]; then
    echo "ERROR: Depth file not found or empty: $depth"
    exit 1
fi

if [ -d "$output_dir" ] && [ -n "$(ls -A "$output_dir")" ]; then
    echo "ERROR: Output directory already contains files: $output_dir"
    exit 1
fi

mkdir -p "$output_dir"

metabat2 \
    -i "$assembly" \
    -a "$depth" \
    -o "${output_dir}/${sample}_bin" \
    -m 1500 \
    -t "$SLURM_CPUS_PER_TASK" \
    --seed 12345

bin_count=$(find "$output_dir" -maxdepth 1 -type f -name "*.fa" | wc -l)

echo "MetaBAT2 completed for $sample"
echo "Bins produced: $bin_count"