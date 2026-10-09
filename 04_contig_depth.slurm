#!/bin/bash
#SBATCH --job-name=contig_depth
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --partition=cpu_shared
#SBATCH --time=04:00:00
#SBATCH --mem=16G
#SBATCH --array=0-3
#SBATCH --output=depth_%A_%a.out
#SBATCH --error=depth_%A_%a.err

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

BAM_ROOT="03.Bowtie2"
OUTPUT_ROOT="04.ContigDepth"

bam="${BAM_ROOT}/${sample}/${sample}.sorted.bam"
output_dir="${OUTPUT_ROOT}/${sample}"
depth_file="${output_dir}/${sample}_depth.txt"

echo "Sample:     ${sample}"
echo "Input BAM:  ${bam}"
echo "Depth file: ${depth_file}"

if [ ! -s "$bam" ]; then
    echo "ERROR: BAM file not found or empty: $bam"
    exit 1
fi

mkdir -p "$output_dir"

jgi_summarize_bam_contig_depths \
    --outputDepth "$depth_file" \
    "$bam"

if [ ! -s "$depth_file" ]; then
    echo "ERROR: Depth file was not produced."
    exit 1
fi

echo "Depth calculation completed for ${sample}"
head "$depth_file"