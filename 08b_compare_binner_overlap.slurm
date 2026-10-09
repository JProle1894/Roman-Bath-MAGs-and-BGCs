#!/bin/bash
#SBATCH --job-name=bin_overlap
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=cpu_shared
#SBATCH --time=02:00:00
#SBATCH --mem=16G
#SBATCH --array=0-3%4
#SBATCH --output=bin_overlap_%A_%a.out
#SBATCH --error=bin_overlap_%A_%a.err

set -eo pipefail

PROJECT="$HOME/MAGs_V2"

samples=(
    "GB_Biofilm"
    "GB_Sediment"
    "HS_Biofilm"
    "HS_Sediment"
)

sample="${samples[$SLURM_ARRAY_TASK_ID]}"

echo "============================================"
echo "Binner overlap analysis"
echo "Sample: ${sample}"
echo "Started: $(date)"
echo "============================================"

# ------------------------------------------------------------
# Find DAS Tool results
# ------------------------------------------------------------

DASTOOL_DIR=""

for root in \
    "${PROJECT}/08.DASTool_AllRaw" \
    "${PROJECT}/08.DASTool"
do
    if [ -d "${root}/${sample}" ]; then
        DASTOOL_DIR="${root}/${sample}"
        break
    fi
done

if [ -z "$DASTOOL_DIR" ]; then
    echo "ERROR: Could not find DAS Tool directory for ${sample}"
    exit 1
fi

echo "Using DAS Tool directory:"
echo "$DASTOOL_DIR"

# ------------------------------------------------------------
# Input files
# ------------------------------------------------------------

assembly="${PROJECT}/02.MEGAHIT/${sample}/final.contigs.fa"

raw_dir="${DASTOOL_DIR}/contigs2bin"

summary="${DASTOOL_DIR}/${sample}_DASTool_summary.tsv"

# DAS Tool normally calls this contigs2bin,
# but check both possible names
final_map=""

for candidate in \
    "${DASTOOL_DIR}/${sample}_DASTool_contigs2bin.tsv" \
    "${DASTOOL_DIR}/${sample}_DASTool_contig2bin.tsv"
do
    if [ -s "$candidate" ]; then
        final_map="$candidate"
        break
    fi
done

# ------------------------------------------------------------
# Check files exist
# ------------------------------------------------------------

if [ ! -s "$assembly" ]; then
    echo "ERROR: Assembly not found:"
    echo "$assembly"
    exit 1
fi

if [ ! -d "$raw_dir" ]; then
    echo "ERROR: contigs2bin directory not found:"
    echo "$raw_dir"
    exit 1
fi

if [ ! -s "$summary" ]; then
    echo "ERROR: DAS Tool summary not found:"
    echo "$summary"
    exit 1
fi

if [ -z "$final_map" ]; then
    echo "ERROR: Final DAS Tool contigs2bin file not found"
    exit 1
fi

for binner in MetaBAT2 SemiBin2 COMEBin VAMB
do
    file="${raw_dir}/${binner}_contigs2bin.tsv"

    if [ ! -s "$file" ]; then
        echo "ERROR: Missing:"
        echo "$file"
        exit 1
    fi
done

# ------------------------------------------------------------
# Output
# ------------------------------------------------------------

OUTROOT="${PROJECT}/08b.BinnerOverlap"
OUTDIR="${OUTROOT}/${sample}"

mkdir -p "$OUTDIR"

# Threshold for calling two bins the same reconstructed genome
THRESHOLD=0.80

# ------------------------------------------------------------
# Run comparison
# ------------------------------------------------------------

python3 - \
    "$sample" \
    "$assembly" \
    "$raw_dir" \
    "$final_map" \
    "$summary" \
    "$OUTDIR" \
    "$THRESHOLD" <<'PY'

import sys
import os
import csv
from collections import defaultdict

sample = sys.argv[1]
assembly = sys.argv[2]
raw_dir = sys.argv[3]
final_map_file = sys.argv[4]
summary_file = sys.argv[5]
outdir = sys.argv[6]
threshold = float(sys.argv[7])

binners = [
    "MetaBAT2",
    "SemiBin2",
    "COMEBin",
    "VAMB"
]

# ============================================================
# Get contig lengths from MEGAHIT assembly
# ============================================================

print("Reading contig lengths...")

contig_lengths = {}

current = None
length = 0

with open(assembly) as handle:

    for line in handle:

        if line.startswith(">"):

            if current is not None:
                contig_lengths[current] = length

            current = line[1:].strip().split()[0]
            length = 0

        else:
            length += len(line.strip())

    if current is not None:
        contig_lengths[current] = length

print(f"Assembly contigs loaded: {len(contig_lengths):,}")

# ============================================================
# Read final DAS Tool contig assignments
# ============================================================

final_of_contig = {}
final_bp = defaultdict(int)

with open(final_map_file) as handle:

    for line in handle:

        if not line.strip():
            continue

        fields = line.rstrip("\n").split("\t")

        if len(fields) < 2:
            continue

        contig = fields[0]
        final_bin = fields[1]

        final_of_contig[contig] = final_bin

        if contig in contig_lengths:
            final_bp[final_bin] += contig_lengths[contig]

print(f"Final DAS bins: {len(final_bp)}")

# ============================================================
# Read DAS Tool summary to determine which binner supplied
# each selected bin
# ============================================================

source_binner = {}

with open(summary_file, newline="") as handle:

    reader = csv.DictReader(handle, delimiter="\t")

    for row in reader:

        bin_id = row.get("bin")
        bin_set = row.get("bin_set")

        if bin_id:
            source_binner[bin_id] = bin_set

# ============================================================
# Compare each raw bin against each final DAS bin
# ============================================================

all_results = {}

for binner in binners:

    filename = os.path.join(
        raw_dir,
        f"{binner}_contigs2bin.tsv"
    )

    print(f"Reading {binner}...")

    candidate_bp = defaultdict(int)

    # shared_bp[(final_bin, candidate_bin)]
    shared_bp = defaultdict(int)

    with open(filename) as handle:

        for line in handle:

            if not line.strip():
                continue

            fields = line.rstrip("\n").split("\t")

            if len(fields) < 2:
                continue

            contig = fields[0]
            candidate = fields[1]

            length = contig_lengths.get(contig, 0)

            candidate_bp[candidate] += length

            if contig in final_of_contig:

                final_bin = final_of_contig[contig]

                shared_bp[(final_bin, candidate)] += length

    # --------------------------------------------------------
    # Find best matching raw bin for every final DAS bin
    # --------------------------------------------------------

    for final_bin in final_bp:

        best = None

        for (fbin, candidate), shared in shared_bp.items():

            if fbin != final_bin:
                continue

            final_size = final_bp[final_bin]
            candidate_size = candidate_bp[candidate]

            if final_size == 0 or candidate_size == 0:
                continue

            # Fraction of final DAS bin represented in candidate
            final_fraction = shared / final_size

            # Fraction of candidate represented in DAS bin
            candidate_fraction = shared / candidate_size

            # Reciprocal overlap
            reciprocal = min(
                final_fraction,
                candidate_fraction
            )

            # Base-pair weighted Jaccard similarity
            union = final_size + candidate_size - shared

            if union > 0:
                jaccard = shared / union
            else:
                jaccard = 0

            result = {
                "candidate": candidate,
                "shared_bp": shared,
                "final_fraction": final_fraction,
                "candidate_fraction": candidate_fraction,
                "reciprocal": reciprocal,
                "jaccard": jaccard,
            }

            if best is None:
                best = result

            else:

                # Prefer highest reciprocal overlap,
                # then highest Jaccard
                old_score = (
                    best["reciprocal"],
                    best["jaccard"],
                    best["shared_bp"]
                )

                new_score = (
                    result["reciprocal"],
                    result["jaccard"],
                    result["shared_bp"]
                )

                if new_score > old_score:
                    best = result

        all_results[(final_bin, binner)] = best

# ============================================================
# LONG FORMAT OUTPUT
# ============================================================

long_file = os.path.join(
    outdir,
    f"{sample}_binner_overlap_long.tsv"
)

with open(long_file, "w", newline="") as handle:

    fields = [
        "Sample",
        "DAS_bin",
        "DAS_source_binner",
        "Binner",
        "Best_matching_bin",
        "DAS_bin_bp",
        "Shared_bp",
        "Fraction_DAS_bin_shared",
        "Fraction_candidate_shared",
        "Reciprocal_overlap",
        "Jaccard_bp",
        "Same_genome_80pct"
    ]

    writer = csv.DictWriter(
        handle,
        fieldnames=fields,
        delimiter="\t"
    )

    writer.writeheader()

    for final_bin in sorted(final_bp):

        for binner in binners:

            best = all_results.get(
                (final_bin, binner)
            )

            if best is None:

                row = {
                    "Sample": sample,
                    "DAS_bin": final_bin,
                    "DAS_source_binner":
                        source_binner.get(final_bin, ""),
                    "Binner": binner,
                    "Best_matching_bin": "",
                    "DAS_bin_bp": final_bp[final_bin],
                    "Shared_bp": 0,
                    "Fraction_DAS_bin_shared": 0,
                    "Fraction_candidate_shared": 0,
                    "Reciprocal_overlap": 0,
                    "Jaccard_bp": 0,
                    "Same_genome_80pct": "NO"
                }

            else:

                same = (
                    "YES"
                    if best["reciprocal"] >= threshold
                    else "NO"
                )

                row = {
                    "Sample": sample,
                    "DAS_bin": final_bin,
                    "DAS_source_binner":
                        source_binner.get(final_bin, ""),
                    "Binner": binner,
                    "Best_matching_bin":
                        best["candidate"],
                    "DAS_bin_bp":
                        final_bp[final_bin],
                    "Shared_bp":
                        best["shared_bp"],
                    "Fraction_DAS_bin_shared":
                        round(best["final_fraction"], 4),
                    "Fraction_candidate_shared":
                        round(best["candidate_fraction"], 4),
                    "Reciprocal_overlap":
                        round(best["reciprocal"], 4),
                    "Jaccard_bp":
                        round(best["jaccard"], 4),
                    "Same_genome_80pct": same
                }

            writer.writerow(row)

# ============================================================
# WIDE / EASY-TO-READ OUTPUT
# ============================================================

wide_file = os.path.join(
    outdir,
    f"{sample}_binner_overlap_summary.tsv"
)

fields = [
    "Sample",
    "DAS_bin",
    "Selected_from",
    "DAS_bin_bp"
]

for binner in binners:

    fields.extend([
        f"{binner}_bin",
        f"{binner}_reciprocal_overlap",
        f"{binner}_same_genome"
    ])

fields.extend([
    "Recovered_by_n_binners",
    "Recovered_by_binners"
])

with open(wide_file, "w", newline="") as handle:

    writer = csv.DictWriter(
        handle,
        fieldnames=fields,
        delimiter="\t"
    )

    writer.writeheader()

    for final_bin in sorted(final_bp):

        row = {
            "Sample": sample,
            "DAS_bin": final_bin,
            "Selected_from":
                source_binner.get(final_bin, ""),
            "DAS_bin_bp":
                final_bp[final_bin]
        }

        recovered = []

        for binner in binners:

            best = all_results.get(
                (final_bin, binner)
            )

            if best is None:

                row[f"{binner}_bin"] = ""
                row[f"{binner}_reciprocal_overlap"] = 0
                row[f"{binner}_same_genome"] = "NO"

                continue

            reciprocal = best["reciprocal"]

            same = reciprocal >= threshold

            row[f"{binner}_bin"] = best["candidate"]

            row[
                f"{binner}_reciprocal_overlap"
            ] = round(reciprocal, 4)

            row[
                f"{binner}_same_genome"
            ] = "YES" if same else "NO"

            if same:
                recovered.append(binner)

        row["Recovered_by_n_binners"] = len(recovered)

        row["Recovered_by_binners"] = ",".join(recovered)

        writer.writerow(row)

print()
print("======================================")
print("Analysis complete")
print("======================================")
print()
print("Long table:")
print(long_file)
print()
print("Summary table:")
print(wide_file)

PY

echo
echo "Finished: $(date)"