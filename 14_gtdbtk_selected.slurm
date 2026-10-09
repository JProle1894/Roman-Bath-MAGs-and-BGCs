#!/bin/bash
#SBATCH --job-name=GTDBTk_selected
#SBATCH --output=GTDBTk_selected_%j.out
#SBATCH --error=GTDBTk_selected_%j.err
#SBATCH --time=48:00:00
#SBATCH --mem=128G
#SBATCH --cpus-per-task=16

set -euo pipefail

cd "$HOME/MAGs_V2"

module load conda/24.9.2

# ------------------------------------------------------------
# SETTINGS
# ------------------------------------------------------------

BASE="$HOME/MAGs_V2"

WORKDIR="$BASE/14.GTDBTk_Selected"
BATCHFILE="$WORKDIR/selected_71_MAGs.batch.tsv"
OUTDIR="$WORKDIR/results"
SCRATCH="$WORKDIR/scratch"

mkdir -p "$WORKDIR"
mkdir -p "$SCRATCH"

# Do not accidentally overwrite an existing GTDB-Tk run
if [ -d "$OUTDIR" ]; then
    echo "ERROR: GTDB-Tk output directory already exists:"
    echo "$OUTDIR"
    echo "Remove or rename it before rerunning."
    exit 1
fi

# Start fresh batch file
> "$BATCHFILE"

echo "========================================"
echo "Preparing selected MAGs for GTDB-Tk"
echo "========================================"
echo


# ------------------------------------------------------------
# FUNCTION TO FIND EACH DAS TOOL FASTA
# ------------------------------------------------------------

add_bin () {

    sample="$1"
    bin="$2"

    sample_dir="$BASE/08.DASTool/$sample"

    if [ ! -d "$sample_dir" ]; then
        echo "ERROR: Sample directory not found:"
        echo "$sample_dir"
        exit 1
    fi

    mapfile -t matches < <(
        find "$sample_dir" -type f \
        \( \
            -name "${bin}.fa" \
            -o -name "${bin}.fasta" \
            -o -name "${bin}.fna" \
            -o -name "${bin}.fa.gz" \
            -o -name "${bin}.fasta.gz" \
            -o -name "${bin}.fna.gz" \
        \) \
        -path "*DASTool_bins*"
    )

    if [ "${#matches[@]}" -eq 0 ]; then

        echo
        echo "ERROR: Could not find:"
        echo "Sample: $sample"
        echo "Bin:    $bin"
        echo
        exit 1

    fi

    if [ "${#matches[@]}" -gt 1 ]; then

        echo
        echo "ERROR: More than one FASTA found for:"
        echo "Sample: $sample"
        echo "Bin:    $bin"
        printf '%s\n' "${matches[@]}"
        echo
        exit 1

    fi

    fasta=$(readlink -f "${matches[0]}")

    # Make unique GTDB-Tk genome ID
    genome_id="${sample}_${bin}"

    # Replace dots with underscores
    genome_id="${genome_id//./_}"

    printf "%s\t%s\n" \
        "$fasta" \
        "$genome_id" \
        >> "$BATCHFILE"

    echo "FOUND: $sample  $bin"
}


# ------------------------------------------------------------
# 71 SELECTED MAGs
# ------------------------------------------------------------

add_bin HS_Sediment SemiBin_97
add_bin GB_Biofilm SemiBin_94
add_bin GB_Sediment SemiBin_82
add_bin GB_Biofilm SemiBin_81
add_bin GB_Biofilm SemiBin_79
add_bin GB_Biofilm SemiBin_77
add_bin GB_Biofilm SemiBin_76
add_bin HS_Sediment SemiBin_67
add_bin GB_Biofilm SemiBin_49_sub
add_bin HS_Sediment SemiBin_44
add_bin HS_Biofilm SemiBin_35
add_bin HS_Biofilm SemiBin_25
add_bin GB_Biofilm SemiBin_25
add_bin GB_Biofilm SemiBin_142
add_bin GB_Sediment SemiBin_124
add_bin GB_Biofilm SemiBin_117
add_bin GB_Biofilm SemiBin_103_sub

add_bin HS_Biofilm HS_Biofilm_bin.67_sub
add_bin HS_Biofilm HS_Biofilm_bin.61
add_bin HS_Biofilm HS_Biofilm_bin.40

add_bin GB_Sediment GB_Sediment_bin.7
add_bin GB_Sediment GB_Sediment_bin.69
add_bin GB_Sediment GB_Sediment_bin.47
add_bin GB_Sediment GB_Sediment_bin.15

add_bin GB_Biofilm GB_Biofilm_bin.70
add_bin GB_Biofilm GB_Biofilm_bin.58
add_bin GB_Biofilm GB_Biofilm_bin.52
add_bin GB_Biofilm GB_Biofilm_bin.20
add_bin GB_Biofilm GB_Biofilm_bin.17

add_bin GB_Biofilm 26_sub
add_bin GB_Sediment 114568_sub
add_bin GB_Sediment 107697_sub
add_bin GB_Biofilm 100174_sub

add_bin GB_Sediment 114762
add_bin GB_Sediment 114649

add_bin GB_Biofilm 101173
add_bin GB_Biofilm 100879
add_bin GB_Biofilm 100488
add_bin GB_Biofilm 100131
add_bin GB_Biofilm 99857
add_bin GB_Biofilm 99524
add_bin GB_Biofilm 99507
add_bin GB_Biofilm 99438

add_bin HS_Biofilm 92833
add_bin HS_Biofilm 92700
add_bin HS_Biofilm 91961
add_bin HS_Biofilm 91221
add_bin HS_Biofilm 87913

add_bin HS_Sediment 84750
add_bin HS_Sediment 84346
add_bin HS_Sediment 84102
add_bin HS_Sediment 80937

add_bin GB_Sediment 71549
add_bin GB_Sediment 71404

add_bin HS_Biofilm 69898
add_bin HS_Sediment 62655
add_bin GB_Biofilm 43561
add_bin HS_Sediment 5397

add_bin GB_Sediment 34
add_bin GB_Sediment 28
add_bin GB_Biofilm 24
add_bin HS_Sediment 10
add_bin HS_Biofilm 9
add_bin HS_Biofilm 8
add_bin GB_Sediment 7
add_bin GB_Sediment 6
add_bin GB_Sediment 3
add_bin GB_Sediment 2
add_bin HS_Biofilm 2
add_bin GB_Biofilm 1
add_bin HS_Sediment 1


# ------------------------------------------------------------
# CHECK THAT ALL 71 WERE FOUND
# ------------------------------------------------------------

N=$(wc -l < "$BATCHFILE")

echo
echo "========================================"
echo "Batch file created"
echo "========================================"
echo "MAGs found: $N"
echo "Expected:   71"
echo

if [ "$N" -ne 71 ]; then
    echo "ERROR: Expected 71 MAGs but found $N."
    exit 1
fi

echo "All 71 selected MAGs successfully located."
echo

echo "Batch file:"
echo "$BATCHFILE"
echo


# ------------------------------------------------------------
# CHECK GTDB-TK
# ------------------------------------------------------------

echo "========================================"
echo "Checking GTDB-Tk"
echo "========================================"

conda run --no-capture-output \
    -n gtdbtk_env \
    gtdbtk --version

conda run --no-capture-output \
    -n gtdbtk_env \
    gtdbtk check_install

echo
echo "GTDB-Tk installation check passed."
echo


# ------------------------------------------------------------
# RUN GTDB-TK
# ------------------------------------------------------------

echo "========================================"
echo "Starting GTDB-Tk classify_wf"
echo "========================================"
echo

conda run --no-capture-output \
    -n gtdbtk_env \
    gtdbtk classify_wf \
    --batchfile "$BATCHFILE" \
    --out_dir "$OUTDIR" \
    --prefix RomanBaths \
    --cpus "$SLURM_CPUS_PER_TASK" \
    --pplacer_cpus 4 \
    --scratch_dir "$SCRATCH"

echo
echo "========================================"
echo "GTDB-Tk COMPLETE"
echo "========================================"
echo
echo "Results directory:"
echo "$OUTDIR"
echo

find "$OUTDIR" -maxdepth 1 -type f -name "*summary.tsv" -print