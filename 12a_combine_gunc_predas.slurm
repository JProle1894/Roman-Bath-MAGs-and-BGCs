#!/bin/bash
#SBATCH --job-name=combine_gunc_predas
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=cpu_shared
#SBATCH --time=00:10:00
#SBATCH --mem=1G
#SBATCH --output=combine_gunc_predas_%j.out
#SBATCH --error=combine_gunc_predas_%j.err

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

output="11.GUNC_PreDAS/all_binners_gunc_predas.tsv"

rm -f "$output"

first_file=true
missing_count=0
files_added=0

echo "========================================"
echo "Combining pre-DAS GUNC results"
echo "Started: $(date)"
echo "========================================"

for binner in "${binners[@]}"; do

    for sample in "${samples[@]}"; do

        sample_dir="11.GUNC_PreDAS/${binner}/${sample}"

        # Find the GUNC maxCSS result table
        file=$(find "$sample_dir" \
            -maxdepth 1 \
            -type f \
            -name "GUNC.*.maxCSS_level.tsv" \
            -print -quit 2>/dev/null)

        if [ -z "$file" ] || [ ! -s "$file" ]; then

            echo "WARNING: No GUNC results found for:"
            echo "  Binner: ${binner}"
            echo "  Sample: ${sample}"

            missing_count=$((missing_count + 1))
            continue

        fi

        echo "Adding: ${file}"

        files_added=$((files_added + 1))

        if $first_file; then

            awk -v b="$binner" -v s="$sample" '
                BEGIN {
                    OFS="\t"
                }

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
                BEGIN {
                    OFS="\t"
                }

                NR>1 && NF>1 {
                    print b,s,$0
                }
            ' "$file" >> "$output"

        fi

    done

done

# ------------------------------------------------
# Check combined file
# ------------------------------------------------

if [ ! -s "$output" ]; then

    echo
    echo "ERROR: No GUNC results were combined."
    exit 1

fi

total_bins=$(tail -n +2 "$output" | wc -l)

echo
echo "========================================"
echo "Combination complete"
echo "Output: ${output}"
echo "GUNC tables added: ${files_added} / 16"
echo "Missing tables: ${missing_count}"
echo "Total bins/genomes: ${total_bins}"
echo "Finished: $(date)"
echo "========================================"

# ------------------------------------------------
# Counts by binner
# ------------------------------------------------

echo
echo "Genomes per binner:"

tail -n +2 "$output" \
    | cut -f1 \
    | sort \
    | uniq -c

# ------------------------------------------------
# Counts by binner + sample
# ------------------------------------------------

echo
echo "Genomes per binner and sample:"

tail -n +2 "$output" \
    | cut -f1,2 \
    | sort \
    | uniq -c

# ------------------------------------------------
# GUNC pass/fail summary
# ------------------------------------------------

echo
echo "Overall GUNC pass/fail:"

awk -F'\t' '
NR==1 {
    for(i=1;i<=NF;i++) {
        if($i=="pass.GUNC")
            pass_col=i
    }
    next
}

pass_col {
    count[$pass_col]++
}

END {
    for(result in count)
        print result, count[result]
}
' "$output"

# ------------------------------------------------
# Pass/fail by binner
# ------------------------------------------------

echo
echo "GUNC pass/fail by binner:"

awk -F'\t' '
BEGIN {
    OFS="\t"
}

NR==1 {
    for(i=1;i<=NF;i++) {
        if($i=="Binner")
            binner_col=i
        if($i=="pass.GUNC")
            pass_col=i
    }
    next
}

{
    key=$binner_col SUBSEP $pass_col
    count[key]++
}

END {
    for(key in count) {
        split(key,a,SUBSEP)
        print a[1],a[2],count[key]
    }
}
' "$output" | sort