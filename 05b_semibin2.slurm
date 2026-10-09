#!/bin/bash
#SBATCH --job-name=semibin2
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --partition=cpu_shared
#SBATCH --time=24:00:00
#SBATCH --mem=64G
#SBATCH --array=0-3%1
#SBATCH --output=logs/semibin2_%A_%a.out
#SBATCH --error=logs/semibin2_%A_%a.err

set -eo pipefail

ROOT="$SLURM_SUBMIT_DIR"
cd "$ROOT"

module load conda/24.9.2
source /scratch/conda/miniconda3/20240912/etc/profile.d/conda.sh
conda activate semibin2_env

samples=(
    "GB_Biofilm"
    "GB_Sediment"
    "HS_Biofilm"
    "HS_Sediment"
)

sample="${samples[$SLURM_ARRAY_TASK_ID]}"

assembly="${ROOT}/02.MEGAHIT/${sample}/final.contigs.fa"
bam="${ROOT}/03.Bowtie2/${sample}/${sample}.sorted.bam"
output_dir="${ROOT}/05.Binning/SemiBin2/${sample}"

echo "Sample:   $sample"
echo "Assembly: $assembly"
echo "BAM:      $bam"
echo "Output:   $output_dir"

command -v SemiBin2 >/dev/null 2>&1 || {
    echo "ERROR: SemiBin2 is not available in the active environment."
    exit 1
}

if [ ! -s "$assembly" ]; then
    echo "ERROR: Assembly not found or empty: $assembly"
    exit 1
fi

if [ ! -s "$bam" ]; then
    echo "ERROR: BAM file not found or empty: $bam"
    exit 1
fi

if [ -d "$output_dir" ] && [ -n "$(ls -A "$output_dir")" ]; then
    echo "ERROR: Output directory already contains files: $output_dir"
    exit 1
fi

SemiBin2 single_easy_bin \
    --input-fasta "$assembly" \
    --input-bam "$bam" \
    --output "$output_dir" \
    --environment global \
    --sequencing-type=short_read \
    --engine=cpu \
    --compression=none \
    --random-seed 12345 \
    --threads "$SLURM_CPUS_PER_TASK"

if [ -d "${output_dir}/output_bins" ]; then
    bins_dir="${output_dir}/output_bins"
elif [ -d "${output_dir}/output_recluster_bins" ]; then
    bins_dir="${output_dir}/output_recluster_bins"
else
    echo "ERROR: SemiBin2 output-bin directory was not found."
    exit 1
fi

bin_count=$(find "$bins_dir" -maxdepth 1 -type f | wc -l)

echo "SemiBin2 completed for $sample"
echo "Bins produced: $bin_count"
echo "Bin directory: $bins_dir"