#!/bin/bash
#SBATCH --job-name=bowtie2
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --partition=cpu_shared
#SBATCH --time=24:00:00
#SBATCH --mem=64G
#SBATCH --array=0-3
#SBATCH --output=bowtie2_%A_%a.out
#SBATCH --error=bowtie2_%A_%a.err

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

READS_ROOT="00.CleanData"
ASSEMBLY_ROOT="02.MEGAHIT"
OUTPUT_ROOT="03.Bowtie2"

reads_dir="${READS_ROOT}/${sample}"
assembly="${ASSEMBLY_ROOT}/${sample}/final.contigs.fa"
output_dir="${OUTPUT_ROOT}/${sample}"

echo "========================================"
echo "Sample:   ${sample}"
echo "Reads:    ${reads_dir}"
echo "Assembly: ${assembly}"
echo "Output:   ${output_dir}"
echo "========================================"

# Check that the required programs are available
command -v bowtie2
command -v bowtie2-build
command -v samtools

# Check the input folders and assembly
if [ ! -d "$reads_dir" ]; then
    echo "ERROR: Reads directory not found: $reads_dir"
    exit 1
fi

if [ ! -s "$assembly" ]; then
    echo "ERROR: Assembly not found or empty: $assembly"
    exit 1
fi

# Locate forward reads
mapfile -t r1_files < <(
    find "$reads_dir" -maxdepth 1 -type f \
        \( -name "*_1.fq.gz" \
        -o -name "*_1.fastq.gz" \
        -o -name "*_1.clean.fq.gz" \
        -o -name "*_1.clean.fastq.gz" \) \
        | sort
)

# Locate reverse reads
mapfile -t r2_files < <(
    find "$reads_dir" -maxdepth 1 -type f \
        \( -name "*_2.fq.gz" \
        -o -name "*_2.fastq.gz" \
        -o -name "*_2.clean.fq.gz" \
        -o -name "*_2.clean.fastq.gz" \) \
        | sort
)

echo "R1 files found: ${#r1_files[@]}"
printf "  %s\n" "${r1_files[@]}"

echo "R2 files found: ${#r2_files[@]}"
printf "  %s\n" "${r2_files[@]}"

if [ "${#r1_files[@]}" -eq 0 ] || \
   [ "${#r2_files[@]}" -eq 0 ]; then
    echo "ERROR: Paired reads were not found for ${sample}"
    exit 1
fi

if [ "${#r1_files[@]}" -ne "${#r2_files[@]}" ]; then
    echo "ERROR: R1 and R2 file counts do not match."
    exit 1
fi

# Convert multiple read files into comma-separated lists
r1_csv=$(IFS=,; echo "${r1_files[*]}")
r2_csv=$(IFS=,; echo "${r2_files[*]}")

# Prevent accidental overwriting
if [ -e "$output_dir" ]; then
    echo "ERROR: Output directory already exists: $output_dir"
    echo "Remove or rename it before rerunning this sample."
    exit 1
fi

mkdir -p "${output_dir}/index"

index_prefix="${output_dir}/index/${sample}"
bam="${output_dir}/${sample}.sorted.bam"
summary="${output_dir}/${sample}_bowtie2_summary.txt"
flagstat="${output_dir}/${sample}_flagstat.txt"

echo "Building Bowtie2 index..."

bowtie2-build \
    "$assembly" \
    "$index_prefix"

echo "Mapping reads..."

# Use 12 threads for Bowtie2 and 4 for BAM sorting
bowtie2 \
    --very-sensitive \
    -x "$index_prefix" \
    -1 "$r1_csv" \
    -2 "$r2_csv" \
    -p 12 \
    2> "$summary" \
    | samtools sort \
        -@ 4 \
        -o "$bam" \
        -

echo "Indexing BAM file..."

samtools index "$bam"

echo "Generating mapping statistics..."

samtools flagstat \
    -@ 4 \
    "$bam" \
    > "$flagstat"

if [ ! -s "$bam" ]; then
    echo "ERROR: Sorted BAM file was not produced."
    exit 1
fi

if [ ! -s "${bam}.bai" ]; then
    echo "ERROR: BAM index was not produced."
    exit 1
fi

echo "========================================"
echo "Bowtie2 completed successfully: ${sample}"
ls -lh "$bam" "${bam}.bai"
echo
cat "$summary"
echo "========================================"