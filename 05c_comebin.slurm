#!/bin/bash
#SBATCH --job-name=comebin
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --partition=cpu_shared
#SBATCH --time=48:00:00
#SBATCH --mem=64G
#SBATCH --array=0-3%1
#SBATCH --output=logs/comebin_%A_%a.out
#SBATCH --error=logs/comebin_%A_%a.err

set -eo pipefail

ROOT="$SLURM_SUBMIT_DIR"
cd "$ROOT"

module load conda/24.9.2
source /scratch/conda/miniconda3/20240912/etc/profile.d/conda.sh
conda activate comebin_env

samples=(
    "GB_Biofilm"
    "GB_Sediment"
    "HS_Biofilm"
    "HS_Sediment"
)

sample="${samples[$SLURM_ARRAY_TASK_ID]}"

assembly="${ROOT}/02.MEGAHIT/${sample}/final.contigs.fa"
bam_dir="${ROOT}/03.Bowtie2/${sample}"
bam="${bam_dir}/${sample}.sorted.bam"
output_dir="${ROOT}/05.Binning/COMEBin/${sample}"

echo "Sample:        $sample"
echo "Assembly:      $assembly"
echo "BAM directory: $bam_dir"
echo "Output:        $output_dir"

command -v run_comebin.sh >/dev/null 2>&1 || {
    echo "ERROR: run_comebin.sh is not available in the active environment."
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

bam_count=$(find "$bam_dir" -maxdepth 1 -type f -name "*.bam" | wc -l)

if [ "$bam_count" -ne 1 ]; then
    echo "ERROR: Expected one BAM file in $bam_dir, but found $bam_count."
    exit 1
fi

if [ -d "$output_dir" ] && [ -n "$(ls -A "$output_dir")" ]; then
    echo "ERROR: Output directory already contains files: $output_dir"
    exit 1
fi

mkdir -p "$output_dir"

run_comebin.sh \
    -a "$assembly" \
    -p "$bam_dir" \
    -o "$output_dir" \
    -t "$SLURM_CPUS_PER_TASK"

bins_dir="${output_dir}/comebin_res/comebin_res_bins"

if [ ! -d "$bins_dir" ]; then
    echo "ERROR: COMEBin output-bin directory was not produced."
    exit 1
fi

bin_count=$(find "$bins_dir" -maxdepth 1 -type f | wc -l)

echo "COMEBin completed for $sample"
echo "Bins produced: $bin_count"
echo "Bin directory: $bins_dir"