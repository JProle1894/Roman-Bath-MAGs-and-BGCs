#!/bin/bash
#SBATCH --job-name=merge_MAG_tables
#SBATCH --output=merge_MAG_tables_%j.out
#SBATCH --error=merge_MAG_tables_%j.err
#SBATCH --time=00:10:00
#SBATCH --mem=2G
#SBATCH --cpus-per-task=1

set -euo pipefail

cd "$HOME/MAGs_V2"

echo "Starting MAG quality merge"
echo "Working directory: $(pwd)"

# ------------------------------------------------------------
# INPUT FILES
# ------------------------------------------------------------

checkm="09.CheckM_DASTool/all_dastool_checkm.tsv"
gunc="11.GUNC_DASTool/all_dastool_gunc.tsv"

summary_files=(08.DASTool/*/*_DASTool_summary.tsv)

# Check files exist
if [ ! -f "$checkm" ]; then
    echo "ERROR: CheckM file not found:"
    echo "$checkm"
    exit 1
fi

if [ ! -f "$gunc" ]; then
    echo "ERROR: GUNC file not found:"
    echo "$gunc"
    exit 1
fi

if [ "${#summary_files[@]}" -eq 0 ]; then
    echo "ERROR: No DAS Tool summary files found."
    exit 1
fi

echo "CheckM file: $checkm"
echo "GUNC file: $gunc"
echo "DAS Tool summaries found: ${#summary_files[@]}"

# ------------------------------------------------------------
# BUILD DAS TOOL BIN PROVENANCE TABLE
# ------------------------------------------------------------

provenance="dastool_bin_provenance.tsv"

echo -e "Sample\tBin_ID\tBinning_method" > "$provenance"

for file in "${summary_files[@]}"
do
    sample=$(basename "$(dirname "$file")")

    awk -v sample="$sample" '
        BEGIN {
            FS=OFS="\t"
        }

        NR > 1 && $1 != "" {
            print sample, $1, $2
        }
    ' "$file" >> "$provenance"
done

echo "Created provenance table: $provenance"

# ------------------------------------------------------------
# MERGE CHECKM + GUNC + BINNING METHOD
# ------------------------------------------------------------

output="all_MAG_quality_with_binner.tsv"

awk \
    -v provenance="$provenance" \
    -v guncfile="$gunc" \
    '
    BEGIN {
        FS=OFS="\t"
    }

    # --------------------------------------------------------
    # READ DAS TOOL PROVENANCE
    # --------------------------------------------------------

    FILENAME == provenance {

        if (FNR == 1)
            next

        key = $1 SUBSEP $2
        method[key] = $3

        next
    }

    # --------------------------------------------------------
    # READ GUNC TABLE
    # --------------------------------------------------------

    FILENAME == guncfile {

        if (FNR == 1) {

            gunc_nf = NF

            for (i = 3; i <= NF; i++)
                gunc_header[i] = $i

            next
        }

        key = $1 SUBSEP $2

        value = ""

        for (i = 3; i <= NF; i++) {

            if (i == 3)
                value = $i
            else
                value = value OFS $i
        }

        gunc_data[key] = value

        next
    }

    # --------------------------------------------------------
    # READ CHECKM AND WRITE FINAL TABLE
    # --------------------------------------------------------

    FNR == 1 {

        printf "Genome_ID\tSample\tBin_ID\tBinning_method"

        for (i = 3; i <= NF; i++)
            printf "\tCheckM_%s", $i

        for (i = 3; i <= gunc_nf; i++)
            printf "\tGUNC_%s", gunc_header[i]

        printf "\n"

        next
    }

    {

        sample = $1
        bin_id = $2

        key = sample SUBSEP bin_id

        if (key in method)
            binner = method[key]
        else
            binner = "UNKNOWN"

        genome_id = sample "_" binner "_" bin_id

        printf "%s\t%s\t%s\t%s",
            genome_id,
            sample,
            bin_id,
            binner

        # Add CheckM fields
        for (i = 3; i <= NF; i++)
            printf "\t%s", $i

        # Add GUNC fields
        if (key in gunc_data)
            printf "\t%s", gunc_data[key]
        else {

            for (i = 3; i <= gunc_nf; i++)
                printf "\t"
        }

        printf "\n"
    }

    ' \
    "$provenance" \
    "$gunc" \
    "$checkm" \
    > "$output"

# ------------------------------------------------------------
# CHECK FINAL OUTPUT
# ------------------------------------------------------------

rows=$(( $(wc -l < "$output") - 1 ))

unknown=$(awk -F'\t' '
    NR > 1 && $4 == "UNKNOWN" {
        n++
    }

    END {
        print n+0
    }
' "$output")

echo
echo "----------------------------------------"
echo "MERGE FINISHED"
echo "----------------------------------------"
echo "Output file: $output"
echo "MAGs written: $rows"
echo "Bins with UNKNOWN binner: $unknown"
echo

echo "Contribution from each binning method:"

awk -F'\t' '
    NR > 1 {
        count[$4]++
    }

    END {
        for (x in count)
            print x, count[x]
    }
' "$output" | sort

echo
echo "First five lines of output:"
head -5 "$output"