#!/bin/bash
#SBATCH --job-name=combine_checkm
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=cpu_shared
#SBATCH --time=00:10:00
#SBATCH --mem=1G
#SBATCH --output=combine_checkm_%j.out
#SBATCH --error=combine_checkm_%j.err

set -eo pipefail

PROJECT="$HOME/MAGs_V2"
cd "$PROJECT"

binners=(
    "MetaBAT2"
    "SemiBin2"
    "COMEBin"
    "VAMB"
)

samples=(
    "GB_Biofilm"
    "GB_Sediment"
    "HS_Biofilm"
    "HS_Sediment"
)

output="06.CheckM/all_binners_checkm.tsv"

# Remove an older combined file if present
rm -f "$output"

first_file=true

echo "========================================"
echo "Combining CheckM results"
echo "Started: $(date)"
echo "========================================"

for binner in "${binners[@]}"; do

    for sample in "${samples[@]}"; do

        # Check both output layouts used by our CheckM scripts
        candidate1="06.CheckM/${binner}/${sample}/checkm_results.tsv"
        candidate2="06.CheckM/${binner}/${sample}_results.tsv"

        if [ -s "$candidate1" ]; then
            file="$candidate1"

        elif [ -s "$candidate2" ]; then
            file="$candidate2"

        else
            echo "WARNING: No CheckM results found for ${binner} ${sample}"
            continue
        fi

        echo "Adding: $file"

        if $first_file; then

            awk -v b="$binner" -v s="$sample" '
                BEGIN {OFS="\t"}
                NR==1 {
                    print "Binner","Sample",$0
                    next
                }
                NF>1 {
                    print b,s,$0
                }
            ' "$file" > "$output"

            first_file=false

        else

            awk -v b="$binner" -v s="$sample" '
                BEGIN {OFS="\t"}
                NR>1 && NF>1 {
                    print b,s,$0
                }
            ' "$file" >> "$output"

        fi

    done

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
echo "Total bins: ${total_bins}"
echo "Finished: $(date)"
echo "========================================"

echo
echo "Bins per binner:"
tail -n +2 "$output" \
    | cut -f1 \
    | sort \
    | uniq -c

echo
echo "Bins per binner and sample:"
tail -n +2 "$output" \
    | cut -f1,2 \
    | sort \
    | uniq -c