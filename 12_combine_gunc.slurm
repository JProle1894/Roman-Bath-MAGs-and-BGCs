#!/bin/bash
#SBATCH --job-name=combine_gunc_das
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=cpu_shared
#SBATCH --time=00:10:00
#SBATCH --mem=1G
#SBATCH --output=combine_gunc_das_%j.out
#SBATCH --error=combine_gunc_das_%j.err

set -eo pipefail

PROJECT="$HOME/MAGs_V2"
cd "$PROJECT"

samples=(
    "GB_Biofilm"
    "GB_Sediment"
    "HS_Biofilm"
    "HS_Sediment"
)

output="11.GUNC_DASTool/all_dastool_gunc.tsv"

rm -f "$output"

first_file=true

echo "========================================"
echo "Combining post-DAS Tool GUNC results"
echo "Started: $(date)"
echo "========================================"

for sample in "${samples[@]}"; do

    sample_dir="11.GUNC_DASTool/${sample}"

    file=$(find "$sample_dir" \
        -maxdepth 1 \
        -type f \
        -name "GUNC.*.maxCSS_level.tsv" \
        -print -quit 2>/dev/null)

    if [ -z "$file" ] || [ ! -s "$file" ]; then
        echo "WARNING: No GUNC results found for ${sample}"
        continue
    fi

    echo "Adding: $file"

    if $first_file; then

        awk -F'\t' -v OFS='\t' -v sample="$sample" '
        NR==1 {
            # Replace the GUNC column name "genome"
            # with the clearer name "Bin_ID"
            for(i=1;i<=NF;i++) {
                if($i=="genome")
                    $i="Bin_ID"
            }

            print "Sample",$0
            next
        }

        NF>1 {
            print sample,$0
        }
        ' "$file" > "$output"

        first_file=false

    else

        awk -F'\t' -v OFS='\t' -v sample="$sample" '
        NR>1 && NF>1 {
            print sample,$0
        }
        ' "$file" >> "$output"

    fi

done

if [ ! -s "$output" ]; then
    echo "ERROR: No GUNC results were combined."
    exit 1
fi

total_bins=$(tail -n +2 "$output" | wc -l)

echo
echo "========================================"
echo "Finished"
echo "Output: $output"
echo "Total DAS Tool bins: $total_bins"
echo "========================================"

echo
echo "Bins per sample:"

tail -n +2 "$output" \
    | cut -f1 \
    | sort \
    | uniq -c

echo
echo "First five rows:"
head -n 6 "$output"