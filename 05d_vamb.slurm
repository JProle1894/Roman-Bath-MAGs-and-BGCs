#!/bin/bash
#SBATCH --job-name=vamb
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --partition=cpu_shared
#SBATCH --time=48:00:00
#SBATCH --mem=64G
#SBATCH --array=0-3%1
#SBATCH --output=vamb_%A_%a.out
#SBATCH --error=vamb_%A_%a.err

set -eo pipefail

cd "$SLURM_SUBMIT_DIR"

module load conda/24.9.2
source /scratch/conda/miniconda3/20240912/etc/profile.d/conda.sh
conda activate vamb_env

command -v vamb >/dev/null 2>&1 || {
    echo "ERROR: VAMB is not available in vamb_env"
    exit 1
}

samples=(
    "GB_Biofilm"
    "GB_Sediment"
    "HS_Biofilm"
    "HS_Sediment"
)

sample="${samples[$SLURM_ARRAY_TASK_ID]}"

ASSEMBLY_ROOT="02.MEGAHIT"
BAM_ROOT="03.Bowtie2"
OUTPUT_ROOT="05.Binning/VAMB"

assembly="${ASSEMBLY_ROOT}/${sample}/final.contigs.fa"
bam_dir="${BAM_ROOT}/${sample}"
bam="${bam_dir}/${sample}.sorted.bam"
output_dir="${OUTPUT_ROOT}/${sample}"

echo "========================================"
echo "Independent VAMB binning"
echo "Sample:   ${sample}"
echo "Assembly: ${assembly}"
echo "BAM:      ${bam}"
echo "Output:   ${output_dir}"
echo "========================================"

if [ ! -s "$assembly" ]; then
    echo "ERROR: Assembly missing or empty: $assembly"
    exit 1
fi

if [ ! -s "$bam" ]; then
    echo "ERROR: Sorted BAM missing or empty: $bam"
    exit 1
fi

bam_count=$(find "$bam_dir" -maxdepth 1 \
    -type f -name "*.bam" | wc -l)

if [ "$bam_count" -ne 1 ]; then
    echo "ERROR: Expected one BAM file in $bam_dir"
    echo "BAM files found: $bam_count"
    exit 1
fi

mkdir -p "$OUTPUT_ROOT"

if [ -e "$output_dir" ]; then
    echo "ERROR: Output directory already exists: $output_dir"
    echo "Remove or rename it before rerunning."
    exit 1
fi

vamb bin default \
    --outdir "$output_dir" \
    --fasta "$assembly" \
    --bamdir "$bam_dir" \
    --minfasta 250000 \
    -p "$SLURM_CPUS_PER_TASK" \
    --seed 12345

bins_dir="${output_dir}/bins"

if [ ! -d "$bins_dir" ]; then
    echo "ERROR: VAMB bins directory was not produced."
    exit 1
fi

bin_count=$(find "$bins_dir" -maxdepth 1 \
    -type f -name "*.fna" | wc -l)

echo "========================================"
echo "VAMB completed successfully"
echo "Sample: ${sample}"
echo "Bins produced: ${bin_count}"
echo "Bin directory: ${bins_dir}"
echo "========================================"