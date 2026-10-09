#!/bin/bash
#SBATCH --job-name=NCBI_ANI_HQ
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --partition=cpu_shared
#SBATCH --time=24:00:00
#SBATCH --mem=64G
#SBATCH --output=NCBI_ANI_HQ_%j.out
#SBATCH --error=NCBI_ANI_HQ_%j.err

set -eo pipefail

# ============================================================
# PATHS
# ============================================================

PROJECT="$HOME/MAGs_V2"

HQ_DIR="${PROJECT}/12.HQ_MAGs"

# Your existing GTDB-Tk results
GTDB_ROOT="${PROJECT}/14.GTDBTk_Selected/results"

OUTDIR="${PROJECT}/15.NCBI_ANI"

mkdir -p "$OUTDIR"

QUERY_LIST="${OUTDIR}/HQ_MAGs.list"
HQ_IDS="${OUTDIR}/HQ_MAG_IDs.txt"

ACCESSIONS="${OUTDIR}/NCBI_reference_accessions.txt"

ZIPFILE="${OUTDIR}/NCBI_references.zip"
NCBI_DIR="${OUTDIR}/NCBI_references"

REF_LIST="${OUTDIR}/NCBI_references.list"

ANI_OUT="${OUTDIR}/HQ_MAGs_vs_NCBI.fastANI.tsv"

BEST_OUT="${OUTDIR}/HQ_MAGs_best_NCBI_matches.tsv"

METADATA="${OUTDIR}/NCBI_reference_metadata.tsv"

# How many "other related references" to retain from GTDB-Tk
# per MAG in addition to closest genome + placement reference
N_RELATED=10


# ============================================================
# CONDA
# ============================================================

cd "$PROJECT"

module load conda/24.9.2
source /scratch/conda/miniconda3/20240912/etc/profile.d/conda.sh

conda activate ncbi_ani_env

command -v datasets >/dev/null 2>&1 || {
    echo "ERROR: datasets not found"
    exit 1
}

command -v dataformat >/dev/null 2>&1 || {
    echo "ERROR: dataformat not found"
    exit 1
}

command -v fastANI >/dev/null 2>&1 || {
    echo "ERROR: fastANI not found"
    exit 1
}


# ============================================================
# CHECK HQ MAG DIRECTORY
# ============================================================

if [ ! -d "$HQ_DIR" ]; then

    echo "ERROR: HQ MAG directory missing:"
    echo "$HQ_DIR"

    exit 1
fi


# ============================================================
# MAKE QUERY LIST
# ============================================================

find "$HQ_DIR" \
    -type f \
    \( \
        -name "*.fa" \
        -o -name "*.fna" \
        -o -name "*.fasta" \
        -o -name "*.fas" \
    \) \
    | sort \
    > "$QUERY_LIST"

QUERY_COUNT=$(wc -l < "$QUERY_LIST")

echo "============================================"
echo "NCBI ANI analysis"
echo "============================================"
echo
echo "HQ MAGs found: $QUERY_COUNT"
echo

if [ "$QUERY_COUNT" -ne 22 ]; then

    echo "WARNING:"
    echo "Expected 22 HQ MAGs but found $QUERY_COUNT."
    echo
fi


# ============================================================
# GET BIN IDS FROM FASTA FILENAMES
# ============================================================

awk -F/ '
{
    file=$NF
    sub(/\.(fa|fna|fasta|fas)$/, "", file)
    print file
}
' "$QUERY_LIST" \
| sort -u \
> "$HQ_IDS"

echo "HQ MAG IDs:"
cat "$HQ_IDS"

echo


# ============================================================
# FIND GTDB-TK SUMMARY FILES
# ============================================================

mapfile -t GTDB_SUMMARIES < <(
    find "$GTDB_ROOT" \
        -type f \
        \( \
            -name "*bac120.summary.tsv" \
            -o -name "*ar53.summary.tsv" \
        \) \
        | sort
)

if [ "${#GTDB_SUMMARIES[@]}" -eq 0 ]; then

    echo "ERROR: No GTDB-Tk summary files found under:"
    echo "$GTDB_ROOT"
    echo
    echo "Search manually with:"
    echo "find $PROJECT -type f -name '*summary.tsv'"

    exit 1
fi

echo "GTDB-Tk summary files:"
printf '  %s\n' "${GTDB_SUMMARIES[@]}"

echo


# ============================================================
# EXTRACT NCBI ACCESSIONS FROM GTDB-TK
#
# Uses:
#   closest_genome_reference
#   closest_placement_reference
#   first 10 other_related_references
# ============================================================

python3 - \
    "$HQ_IDS" \
    "$ACCESSIONS" \
    "$N_RELATED" \
    "${GTDB_SUMMARIES[@]}" <<'PY'

import csv
import re
import sys

hq_file = sys.argv[1]
out_file = sys.argv[2]
n_related = int(sys.argv[3])
summary_files = sys.argv[4:]

with open(hq_file) as fh:
    hq_ids = {line.strip() for line in fh if line.strip()}

accessions = set()
found_genomes = set()

acc_regex = re.compile(r'GC[AF]_\d+\.\d+')

for summary in summary_files:

    print(f"Reading: {summary}")

    with open(summary, newline="") as fh:

        reader = csv.DictReader(fh, delimiter="\t")

        for row in reader:

            genome = row.get("user_genome", "").strip()

            if genome not in hq_ids:
                continue

            found_genomes.add(genome)

            # Closest ANI genome
            closest = row.get(
                "closest_genome_reference", ""
            )

            if closest:
                matches = acc_regex.findall(closest)

                for acc in matches:
                    accessions.add(acc)

            # Closest phylogenetic placement
            placement = row.get(
                "closest_placement_reference", ""
            )

            if placement:
                matches = acc_regex.findall(placement)

                for acc in matches:
                    accessions.add(acc)

            # Other related genomes
            related = row.get(
                "other_related_references", ""
            )

            matches = acc_regex.findall(related)

            for acc in matches[:n_related]:
                accessions.add(acc)


print()
print(f"HQ genomes found in GTDB-Tk: {len(found_genomes)} / {len(hq_ids)}")

missing = hq_ids - found_genomes

if missing:

    print()
    print("WARNING: HQ MAGs absent from GTDB-Tk results:")

    for genome in sorted(missing):
        print("  " + genome)

print()
print(f"Unique NCBI reference accessions: {len(accessions)}")

with open(out_file, "w") as fh:

    for acc in sorted(accessions):
        fh.write(acc + "\n")

PY


# ============================================================
# CHECK ACCESSION LIST
# ============================================================

REF_ACCESSION_COUNT=$(wc -l < "$ACCESSIONS")

echo
echo "============================================"
echo "NCBI references to download: $REF_ACCESSION_COUNT"
echo "============================================"
echo

cat "$ACCESSIONS"

if [ "$REF_ACCESSION_COUNT" -eq 0 ]; then

    echo "ERROR: No NCBI accessions extracted."
    exit 1
fi


# ============================================================
# DOWNLOAD NCBI GENOMES
# ============================================================

echo
echo "Downloading NCBI genomes..."
echo

rm -rf "$NCBI_DIR"
rm -f "$ZIPFILE"

datasets download genome accession \
    --inputfile "$ACCESSIONS" \
    --include genome \
    --filename "$ZIPFILE" \
    --no-progressbar

if [ ! -s "$ZIPFILE" ]; then

    echo "ERROR: NCBI download failed."
    exit 1
fi

mkdir -p "$NCBI_DIR"

unzip -q "$ZIPFILE" \
    -d "$NCBI_DIR"


# ============================================================
# NCBI METADATA
# ============================================================

REPORT="${NCBI_DIR}/ncbi_dataset/data/assembly_data_report.jsonl"

if [ -s "$REPORT" ]; then

    dataformat tsv genome \
        --inputfile "$REPORT" \
        --fields \
accession,organism-name,type_material-display_text,source_database \
        > "$METADATA"

else

    echo "WARNING: NCBI metadata report not found."
fi


# ============================================================
# MAKE NCBI REFERENCE LIST
# ============================================================

find "$NCBI_DIR" \
    -type f \
    -name "*_genomic.fna" \
    | sort \
    > "$REF_LIST"

REF_COUNT=$(wc -l < "$REF_LIST")

echo
echo "Downloaded NCBI genome FASTAs: $REF_COUNT"
echo

if [ "$REF_COUNT" -eq 0 ]; then

    echo "ERROR: No NCBI genome FASTAs found."
    exit 1
fi


# ============================================================
# FASTANI
# ============================================================

echo
echo "============================================"
echo "Running FastANI"
echo "Queries:    $QUERY_COUNT"
echo "References: $REF_COUNT"
echo "============================================"
echo

fastANI \
    --ql "$QUERY_LIST" \
    --rl "$REF_LIST" \
    -t "$SLURM_CPUS_PER_TASK" \
    -o "$ANI_OUT"


# ============================================================
# BEST NCBI MATCH FOR EACH MAG
# ============================================================

python3 - \
    "$ANI_OUT" \
    "$METADATA" \
    "$BEST_OUT" <<'PY'

import csv
import os
import sys

ani_file = sys.argv[1]
metadata_file = sys.argv[2]
output = sys.argv[3]

# ------------------------------------------------------------
# Read NCBI metadata
# ------------------------------------------------------------

metadata = {}

if os.path.exists(metadata_file):

    with open(metadata_file, newline="") as fh:

        reader = csv.DictReader(fh, delimiter="\t")

        for row in reader:

            accession = (
                row.get("Assembly Accession")
                or row.get("accession")
                or ""
            )

            if accession:
                metadata[accession] = row


# ------------------------------------------------------------
# Parse FastANI results
# ------------------------------------------------------------

best = {}

if os.path.exists(ani_file):

    with open(ani_file) as fh:

        for line in fh:

            parts = line.rstrip().split("\t")

            if len(parts) < 5:
                continue

            query = parts[0]
            ref = parts[1]
            ani = float(parts[2])
            mapped = int(parts[3])
            total = int(parts[4])

            af = mapped / total if total else 0

            # NCBI directory immediately above FASTA
            accession = os.path.basename(
                os.path.dirname(ref)
            )

            if (
                query not in best
                or ani > best[query]["ani"]
            ):

                best[query] = {
                    "accession": accession,
                    "ani": ani,
                    "mapped": mapped,
                    "total": total,
                    "af": af
                }


# ------------------------------------------------------------
# Output
# ------------------------------------------------------------

with open(output, "w", newline="") as fh:

    writer = csv.writer(
        fh,
        delimiter="\t"
    )

    writer.writerow([
        "Sample",
        "MAG",
        "Best_NCBI_accession",
        "NCBI_organism",
        "ANI",
        "Alignment_fraction",
        "Mapped_fragments",
        "Total_fragments"
    ])

    for query in sorted(best):

        result = best[query]

        mag = os.path.basename(query)

        for ext in [
            ".fasta",
            ".fna",
            ".fas",
            ".fa"
        ]:
            if mag.endswith(ext):
                mag = mag[:-len(ext)]

        sample = os.path.basename(
            os.path.dirname(query)
        )

        accession = result["accession"]

        meta = metadata.get(accession, {})

        organism = (
            meta.get("Organism Name")
            or meta.get("organism-name")
            or ""
        )

        writer.writerow([
            sample,
            mag,
            accession,
            organism,
            round(result["ani"], 4),
            round(result["af"], 4),
            result["mapped"],
            result["total"]
        ])

PY


# ============================================================
# FINISH
# ============================================================

echo
echo "============================================"
echo "NCBI ANI comparison COMPLETE"
echo "============================================"
echo
echo "All comparisons:"
echo "$ANI_OUT"
echo
echo "Best NCBI match for each MAG:"
echo "$BEST_OUT"
echo
echo "NCBI metadata:"
echo "$METADATA"
echo
echo "Finished: $(date)"
echo

echo "Preview:"
column -t -s $'\t' "$BEST_OUT" | head -30