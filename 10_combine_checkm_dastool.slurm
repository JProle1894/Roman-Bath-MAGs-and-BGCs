#!/bin/bash
#SBATCH --job-name=combine_das_checkm
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=cpu_shared
#SBATCH --time=00:10:00
#SBATCH --mem=1G
#SBATCH --output=combine_das_checkm_%j.out
#SBATCH --error=combine_das_checkm_%j.err

set -eo pipefail

PROJECT="$HOME/MAGs_V2"
cd "$PROJECT"

samples=(
    "GB_Biofilm"
    "GB_Sediment"
    "HS_Biofilm"
    "HS_Sediment"
)

output="09.CheckM_DASTool/all_dastool_checkm.tsv"

rm -f "$output"

first_file=true

echo "========================================"
echo "Combining DAS Tool CheckM results"
echo "Started: $(date)"
echo "========================================"

for sample in "${samples[@]}"; do

    file="09.CheckM_DASTool/${sample}/checkm_results.tsv"

    if [ ! -s "$file" ]; then
        echo "WARNING: No CheckM results found for ${sample}"
        echo "Expected: ${file}"
        continue
    fi

    echo "Adding: ${file}"

    if $first_file; then

        awk -v s="$sample" '
            BEGIN {OFS="\t"}
            NR==1 {
                print "Sample",$0
                next
            }
            NF>1 {
                print s,$0
            }
        ' "$file" > "$output"

        first_file=false

    else

        awk -v s="$sample" '
            BEGIN {OFS="\t"}
            NR>1 && NF>1 {
                print s,$0
            }
        ' "$file" >> "$output"

    fi

done

if [ ! -s "$output" ]; then
    echo "ERROR: No CheckM results were combined."
    exit 1
fi

total_bins=$(tail -n +2 "$output" | wc -l)

echo
echo "========================================"
echo "Combination complete"
echo "Output: ${output}"
echo "Total DAS Tool bins: ${total_bins}"
echo "Finished: $(date)"
echo "========================================"

echo
echo "Bins per sample:"

tail -n +2 "$output" \
    | cut -f1 \
    | sort \
    | uniq -c