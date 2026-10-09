#!/bin/bash
#SBATCH --job-name=GTDBTk_22
#SBATCH --output=16.GTDBTk_MIMAG_Selected/gtdbtk_%j.out
#SBATCH --error=16.GTDBTk_MIMAG_Selected/gtdbtk_%j.err
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=48:00:00

set -euo pipefail


# ============================================================
# PROJECT PATHS
# ============================================================

PROJECT="/users/jprole/MAGs_V2"

BIN_ROOT="${PROJECT}/08.DASTool_AllRaw"

OUTDIR="${PROJECT}/16.GTDBTk_MIMAG_Selected/results"

WORKDIR="${PROJECT}/16.GTDBTk_MIMAG_Selected"

BATCHFILE="${WORKDIR}/selected_22_MAGs.batch.tsv"

MAPFILE="${WORKDIR}/selected_22_MAGs_ID_mapping.tsv"

SCRATCH="${WORKDIR}/scratch"


mkdir -p "${WORKDIR}"
mkdir -p "${SCRATCH}"


# ============================================================
# LOAD CONDA
# ============================================================

source /scratch/conda/miniconda3/20240912/etc/profile.d/conda.sh

conda activate gtdbtk_env


# ============================================================
# CHECK GTDB-TK
# ============================================================

echo "============================================================"
echo "GTDB-Tk selected MAG classification"
echo "============================================================"
echo

echo "GTDB-Tk version:"
gtdbtk --version

echo


# GTDB-Tk requires its reference database to be defined by
# GTDBTK_DATA_PATH.
#
# If your gtdbtk_env was already configured previously,
# this should already exist.

if [[ -z "${GTDBTK_DATA_PATH:-}" ]]; then

    echo "ERROR:"
    echo "GTDBTK_DATA_PATH is not set."
    echo
    echo "Your GTDB-Tk reference database needs to be configured"
    echo "before this job can run."

    exit 1

fi


echo "GTDB-Tk database:"
echo "${GTDBTK_DATA_PATH}"

echo


# ============================================================
# PREVENT ACCIDENTAL OVERWRITE
# ============================================================

if [[ -d "${OUTDIR}" ]]; then

    echo "ERROR: GTDB-Tk output directory already exists:"
    echo "${OUTDIR}"
    echo
    echo "Remove or rename it before rerunning."

    exit 1

fi


# ============================================================
# SELECTED 22 MAGS
#
# Format:
#
# Sample<TAB>Genome
# ============================================================

readarray -t MAGS <<'EOF'
GB_Biofilm	99507
GB_Biofilm	100879
GB_Biofilm	100174_sub
GB_Biofilm	87587_sub
GB_Biofilm	SemiBin_142
GB_Biofilm	SemiBin_49_sub
GB_Biofilm	SemiBin_76
GB_Biofilm	SemiBin_78_sub
GB_Biofilm	SemiBin_79
GB_Sediment	28
GB_Sediment	114649
GB_Sediment	114568_sub
GB_Sediment	GB_Sediment_bin.15
GB_Sediment	GB_Sediment_bin.69
HS_Biofilm	2
HS_Biofilm	9
HS_Biofilm	92833
HS_Biofilm	HS_Biofilm_bin.12
HS_Biofilm	HS_Biofilm_bin.61
HS_Sediment	1
HS_Sediment	37083
HS_Sediment	SemiBin_68
EOF


# ============================================================
# CREATE GTDB-TK BATCH FILE
# ============================================================

echo "Creating batch file..."

> "${BATCHFILE}"

printf "GTDB_ID\tSample\tOriginal_Genome\tFASTA\n" > "${MAPFILE}"


for ENTRY in "${MAGS[@]}"
do

    IFS=$'\t' read -r SAMPLE GENOME <<< "${ENTRY}"

    SAMPLE_DIR="${BIN_ROOT}/${SAMPLE}/${SAMPLE}_DASTool_bins"


    echo
    echo "------------------------------------------------------------"
    echo "Sample: ${SAMPLE}"
    echo "Genome: ${GENOME}"


    # ========================================================
    # CHECK SAMPLE DIRECTORY EXISTS
    # ========================================================

    if [[ ! -d "${SAMPLE_DIR}" ]]; then

        echo "ERROR: Sample directory not found:"
        echo "${SAMPLE_DIR}"

        exit 1

    fi


    # ========================================================
    # FIND FASTA
    # ========================================================

    FASTA=""


    # First try exact filenames

    for EXT in fa fna fasta
    do

        TEST="${SAMPLE_DIR}/${GENOME}.${EXT}"

        if [[ -f "${TEST}" ]]; then
            FASTA="${TEST}"
            break
        fi

    done


    # Also allow gzipped FASTAs

    if [[ -z "${FASTA}" ]]; then

        for EXT in fa.gz fna.gz fasta.gz
        do

            TEST="${SAMPLE_DIR}/${GENOME}.${EXT}"

            if [[ -f "${TEST}" ]]; then
                FASTA="${TEST}"
                break
            fi

        done

    fi


    # ========================================================
    # STOP IF FASTA CANNOT BE FOUND
    # ========================================================

    if [[ -z "${FASTA}" ]]; then

        echo
        echo "ERROR: FASTA not found for:"
        echo "Sample: ${SAMPLE}"
        echo "Genome: ${GENOME}"
        echo
        echo "Expected somewhere in:"
        echo "${SAMPLE_DIR}"

        exit 1

    fi


    echo "FASTA: ${FASTA}"


    # ========================================================
    # CREATE SAFE UNIQUE GTDB ID
    #
    # Example:
    #
    # GB_Biofilm + GB_Biofilm_bin.15
    #
    # becomes:
    #
    # GB_Biofilm_GB_Biofilm_bin_15
    #
    # This prevents punctuation or repeated bin numbers from
    # causing ambiguity.
    # ========================================================

    SAFE_GENOME=$(echo "${GENOME}" | sed 's/[^A-Za-z0-9_]/_/g')

    GTDB_ID="${SAMPLE}_${SAFE_GENOME}"


    # GTDB-Tk batch file:
    #
    # FASTA_PATH<TAB>GENOME_ID

    printf "%s\t%s\n" \
        "${FASTA}" \
        "${GTDB_ID}" \
        >> "${BATCHFILE}"


    # Keep mapping back to original names

    printf "%s\t%s\t%s\t%s\n" \
        "${GTDB_ID}" \
        "${SAMPLE}" \
        "${GENOME}" \
        "${FASTA}" \
        >> "${MAPFILE}"

done


# ============================================================
# CHECK NUMBER OF GENOMES
# ============================================================

N_GENOMES=$(wc -l < "${BATCHFILE}")


echo
echo "============================================================"
echo "Batch file created"
echo "============================================================"
echo
echo "Number of genomes: ${N_GENOMES}"
echo "Batch file:"
echo "${BATCHFILE}"
echo
echo "Genome ID mapping:"
echo "${MAPFILE}"
echo


if [[ "${N_GENOMES}" -ne 22 ]]; then

    echo "ERROR:"
    echo "Expected 22 genomes but batch file contains ${N_GENOMES}"

    exit 1

fi


echo "Batch file contents:"
cat "${BATCHFILE}"

echo


# ============================================================
# RUN GTDB-TK
# ============================================================

echo "============================================================"
echo "Starting GTDB-Tk classify_wf"
echo "============================================================"

gtdbtk classify_wf \
    --batchfile "${BATCHFILE}" \
    --out_dir "${OUTDIR}" \
    --cpus "${SLURM_CPUS_PER_TASK}" \
    --pplacer_cpus "${SLURM_CPUS_PER_TASK}" \
    --scratch_dir "${SCRATCH}"


# ============================================================
# FINISHED
# ============================================================

echo
echo "============================================================"
echo "GTDB-Tk FINISHED"
echo "============================================================"
echo
echo "Results directory:"
echo "${OUTDIR}"
echo
echo "Genome mapping file:"
echo "${MAPFILE}"
echo
echo "Look for classification summaries:"
echo

find "${OUTDIR}" \
    -type f \
    \( -name "*bac120.summary.tsv" \
       -o -name "*ar53.summary.tsv" \) \
    -print

echo
echo "============================================================"