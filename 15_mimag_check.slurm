#!/bin/bash
#SBATCH --job-name=MIMAG_check
#SBATCH --output=MIMAG_check_%j.out
#SBATCH --error=MIMAG_check_%j.err
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=24:00:00

# ============================================================
# MIMAG QUALITY CHECK
#
# Uses:
#   - Existing CheckM results for completeness/contamination
#   - Barrnap for 5S, 16S and 23S rRNA genes
#   - tRNAscan-SE for tRNA amino-acid types
#
# MIMAG high-quality requirements:
#   Completeness >90%
#   Contamination <5%
#   5S rRNA present
#   16S rRNA present
#   23S rRNA present
#   tRNAs representing >=18 of 20 amino acids
# ============================================================


set -o pipefail


# ============================================================
# PATHS
# ============================================================

PROJECT="/users/jprole/MAGs_V2"

BIN_ROOT="${PROJECT}/08.DASTool_AllRaw"

CHECKM_ROOT="${PROJECT}/09.CheckM_DASTool_AllRaw"

OUTDIR="${PROJECT}/15.MIMAG_RNA"

mkdir -p "${OUTDIR}"
mkdir -p "${OUTDIR}/barrnap"
mkdir -p "${OUTDIR}/trnascan"
mkdir -p "${OUTDIR}/logs"


# ============================================================
# CONDA
# ============================================================

source ~/miniconda3/etc/profile.d/conda.sh

conda activate mimag_check


# ============================================================
# CHECK PROGRAMS
# ============================================================

if ! command -v barrnap >/dev/null 2>&1; then
    echo "ERROR: barrnap not found in environment mimag_check"
    exit 1
fi

if ! command -v tRNAscan-SE >/dev/null 2>&1; then
    echo "ERROR: tRNAscan-SE not found in environment mimag_check"
    exit 1
fi


echo "============================================================"
echo "MIMAG QUALITY CHECK"
echo "============================================================"
echo
echo "Barrnap:"
barrnap --version
echo
echo "tRNAscan-SE:"
tRNAscan-SE --version 2>&1 | head -1
echo


# ============================================================
# GENOME LIST
# ============================================================

readarray -t MAGS <<'EOF'
GB_Biofilm	1
GB_Biofilm	24
GB_Biofilm	25
GB_Biofilm	43561
GB_Biofilm	84913
GB_Biofilm	99438
GB_Biofilm	99507
GB_Biofilm	99524
GB_Biofilm	99857
GB_Biofilm	100131
GB_Biofilm	100488
GB_Biofilm	100879
GB_Biofilm	101173
GB_Biofilm	100174_sub
GB_Biofilm	26_sub
GB_Biofilm	61871_sub
GB_Biofilm	87587_sub
GB_Biofilm	GB_Biofilm_bin.17
GB_Biofilm	GB_Biofilm_bin.20
GB_Biofilm	GB_Biofilm_bin.25
GB_Biofilm	GB_Biofilm_bin.41
GB_Biofilm	GB_Biofilm_bin.52
GB_Biofilm	GB_Biofilm_bin.58
GB_Biofilm	GB_Biofilm_bin.70
GB_Biofilm	SemiBin_103_sub
GB_Biofilm	SemiBin_117
GB_Biofilm	SemiBin_142
GB_Biofilm	SemiBin_25
GB_Biofilm	SemiBin_49_sub
GB_Biofilm	SemiBin_76
GB_Biofilm	SemiBin_77
GB_Biofilm	SemiBin_78_sub
GB_Biofilm	SemiBin_79
GB_Biofilm	SemiBin_81
GB_Biofilm	SemiBin_94
GB_Sediment	2
GB_Sediment	3
GB_Sediment	6
GB_Sediment	7
GB_Sediment	28
GB_Sediment	34
GB_Sediment	71404
GB_Sediment	71549
GB_Sediment	76311
GB_Sediment	114649
GB_Sediment	114762
GB_Sediment	107697_sub
GB_Sediment	114568_sub
GB_Sediment	GB_Sediment_bin.15
GB_Sediment	GB_Sediment_bin.47
GB_Sediment	GB_Sediment_bin.53
GB_Sediment	GB_Sediment_bin.69
GB_Sediment	GB_Sediment_bin.7
GB_Sediment	SemiBin_124
GB_Sediment	SemiBin_37_sub
GB_Sediment	SemiBin_82
HS_Biofilm	2
HS_Biofilm	8
HS_Biofilm	9
HS_Biofilm	69898
HS_Biofilm	87913
HS_Biofilm	91221
HS_Biofilm	91961
HS_Biofilm	92700
HS_Biofilm	92833
HS_Biofilm	HS_Biofilm_bin.12
HS_Biofilm	HS_Biofilm_bin.40
HS_Biofilm	HS_Biofilm_bin.61
HS_Biofilm	HS_Biofilm_bin.67_sub
HS_Biofilm	SemiBin_25
HS_Biofilm	SemiBin_35
HS_Sediment	1
HS_Sediment	10
HS_Sediment	13
HS_Sediment	5397
HS_Sediment	37083
HS_Sediment	62655
HS_Sediment	75559
HS_Sediment	80937
HS_Sediment	81275
HS_Sediment	84102
HS_Sediment	84346
HS_Sediment	84750
HS_Sediment	80022_sub
HS_Sediment	HS_Sediment_bin.62_sub
HS_Sediment	SemiBin_22
HS_Sediment	SemiBin_38
HS_Sediment	SemiBin_44
HS_Sediment	SemiBin_67
HS_Sediment	SemiBin_68
HS_Sediment	SemiBin_97
EOF


# ============================================================
# OUTPUT SUMMARY
# ============================================================

SUMMARY="${OUTDIR}/MIMAG_quality_summary.tsv"

printf "Sample\tGenome\tCompleteness\tContamination\t5S\t5S_copies\t16S\t16S_copies\t23S\t23S_copies\ttRNA_total\ttRNA_AA_types\ttRNA_AA_list\tRNA_PASS\tMIMAG_HQ\tStatus\n" \
    > "${SUMMARY}"


# ============================================================
# FUNCTION: GET CHECKM RESULTS
# ============================================================

get_checkm_values() {

    local CHECKM_FILE="$1"
    local TARGET="$2"

    if [[ ! -f "${CHECKM_FILE}" ]]; then
        printf "NA\tNA\n"
        return
    fi

    awk -F '\t' -v target="${TARGET}" '

    function trim(x) {
        gsub(/^[[:space:]]+/, "", x)
        gsub(/[[:space:]]+$/, "", x)
        return x
    }

    NR == 1 {

        bin_col = 1
        comp_col = 0
        contam_col = 0

        for (i = 1; i <= NF; i++) {

            header = trim($i)

            if (header == "Bin Id" || header == "Bin ID") {
                bin_col = i
            }

            if (header == "Completeness") {
                comp_col = i
            }

            if (header == "Contamination") {
                contam_col = i
            }
        }

        next
    }


    comp_col > 0 && contam_col > 0 {

        id = trim($bin_col)

        stem = id

        sub(/\.(fa|fna|fasta)$/, "", stem)

        if (id == target || stem == target) {

            print trim($comp_col) "\t" trim($contam_col)

            found = 1
            exit
        }
    }


    END {

        if (!found) {
            print "NA\tNA"
        }
    }

    ' "${CHECKM_FILE}"
}


# ============================================================
# LOOP THROUGH ALL 91 MAGS
# ============================================================

TOTAL=${#MAGS[@]}

echo "Number of MAGs to analyse: ${TOTAL}"
echo


for ((i=0; i<TOTAL; i++)); do

    IFS=$'\t' read -r SAMPLE GENOME <<< "${MAGS[$i]}"

    ID="${SAMPLE}__${GENOME}"

    echo
    echo "============================================================"
    echo "MAG $((i+1)) / ${TOTAL}"
    echo "Sample: ${SAMPLE}"
    echo "Genome: ${GENOME}"
    echo "============================================================"


    # ========================================================
    # SAMPLE FASTA DIRECTORY
    # ========================================================

    SAMPLE_DIR="${BIN_ROOT}/${SAMPLE}/${SAMPLE}_DASTool_bins"

    if [[ ! -d "${SAMPLE_DIR}" ]]; then

        echo "ERROR: directory not found:"
        echo "${SAMPLE_DIR}"

        printf "%s\t%s\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tDIRECTORY_NOT_FOUND\n" \
            "${SAMPLE}" "${GENOME}" >> "${SUMMARY}"

        continue
    fi


    # ========================================================
    # FIND FASTA
    # ========================================================

    FASTA=""

    for EXT in fa fna fasta; do

        if [[ -f "${SAMPLE_DIR}/${GENOME}.${EXT}" ]]; then

            FASTA="${SAMPLE_DIR}/${GENOME}.${EXT}"
            break

        fi

    done


    # Check gzipped versions as well
    if [[ -z "${FASTA}" ]]; then

        for EXT in fa.gz fna.gz fasta.gz; do

            if [[ -f "${SAMPLE_DIR}/${GENOME}.${EXT}" ]]; then

                FASTA="${SAMPLE_DIR}/${GENOME}.${EXT}"
                break

            fi

        done

    fi


    if [[ -z "${FASTA}" ]]; then

        echo "ERROR: FASTA not found for ${GENOME}"
        echo "Expected in:"
        echo "${SAMPLE_DIR}"

        printf "%s\t%s\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tFASTA_NOT_FOUND\n" \
            "${SAMPLE}" "${GENOME}" >> "${SUMMARY}"

        continue
    fi


    echo "FASTA found:"
    echo "${FASTA}"


    # ========================================================
    # HANDLE GZIPPED FASTA
    # ========================================================

    INPUT_FASTA="${FASTA}"
    TEMP_FASTA=""

    if [[ "${FASTA}" == *.gz ]]; then

        TEMP_FASTA="${OUTDIR}/${ID}.temporary.fna"

        gunzip -c "${FASTA}" > "${TEMP_FASTA}"

        INPUT_FASTA="${TEMP_FASTA}"

    fi


    # ========================================================
    # GET EXISTING CHECKM RESULTS
    # ========================================================

    CHECKM_FILE="${CHECKM_ROOT}/${SAMPLE}/checkm_results.tsv"

    CHECKM_VALUES=$(get_checkm_values "${CHECKM_FILE}" "${GENOME}")

    COMPLETENESS=$(echo "${CHECKM_VALUES}" | cut -f1)
    CONTAMINATION=$(echo "${CHECKM_VALUES}" | cut -f2)

    echo
    echo "Existing CheckM results:"
    echo "Completeness:  ${COMPLETENESS}"
    echo "Contamination: ${CONTAMINATION}"


    # ========================================================
    # BARRNAP
    #
    # Search with BOTH bacterial and archaeal rRNA models.
    # Presence is accepted if either model detects the gene.
    #
    # Copy number is taken as the larger count rather than
    # adding the two models together, avoiding duplicate counts.
    # ========================================================

    BAC_GFF="${OUTDIR}/barrnap/${ID}.bacteria.gff"
    ARC_GFF="${OUTDIR}/barrnap/${ID}.archaea.gff"


    echo
    echo "Running Barrnap bacterial model..."

    if ! barrnap \
        --kingdom bac \
        --threads "${SLURM_CPUS_PER_TASK}" \
        "${INPUT_FASTA}" \
        > "${BAC_GFF}" \
        2> "${OUTDIR}/logs/${ID}.barrnap_bac.log"
    then

        echo "ERROR: bacterial Barrnap search failed"

        printf "%s\t%s\t%s\t%s\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tBARRNAP_ERROR\n" \
            "${SAMPLE}" \
            "${GENOME}" \
            "${COMPLETENESS}" \
            "${CONTAMINATION}" \
            >> "${SUMMARY}"

        [[ -n "${TEMP_FASTA}" ]] && rm -f "${TEMP_FASTA}"

        continue
    fi


    echo "Running Barrnap archaeal model..."

    if ! barrnap \
        --kingdom arc \
        --threads "${SLURM_CPUS_PER_TASK}" \
        "${INPUT_FASTA}" \
        > "${ARC_GFF}" \
        2> "${OUTDIR}/logs/${ID}.barrnap_arc.log"
    then

        echo "ERROR: archaeal Barrnap search failed"

        printf "%s\t%s\t%s\t%s\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tBARRNAP_ERROR\n" \
            "${SAMPLE}" \
            "${GENOME}" \
            "${COMPLETENESS}" \
            "${CONTAMINATION}" \
            >> "${SUMMARY}"

        [[ -n "${TEMP_FASTA}" ]] && rm -f "${TEMP_FASTA}"

        continue
    fi


    # ========================================================
    # COUNT RRNA GENES
    # ========================================================

    BAC_5S=$(grep -c "Name=5S_rRNA" "${BAC_GFF}" || true)
    BAC_16S=$(grep -c "Name=16S_rRNA" "${BAC_GFF}" || true)
    BAC_23S=$(grep -c "Name=23S_rRNA" "${BAC_GFF}" || true)

    ARC_5S=$(grep -c "Name=5S_rRNA" "${ARC_GFF}" || true)
    ARC_16S=$(grep -c "Name=16S_rRNA" "${ARC_GFF}" || true)
    ARC_23S=$(grep -c "Name=23S_rRNA" "${ARC_GFF}" || true)


    # Take maximum count from either model
    if (( BAC_5S > ARC_5S )); then
        N_5S=${BAC_5S}
    else
        N_5S=${ARC_5S}
    fi

    if (( BAC_16S > ARC_16S )); then
        N_16S=${BAC_16S}
    else
        N_16S=${ARC_16S}
    fi

    if (( BAC_23S > ARC_23S )); then
        N_23S=${BAC_23S}
    else
        N_23S=${ARC_23S}
    fi


    if (( N_5S > 0 )); then
        RRNA_5S="YES"
    else
        RRNA_5S="NO"
    fi

    if (( N_16S > 0 )); then
        RRNA_16S="YES"
    else
        RRNA_16S="NO"
    fi

    if (( N_23S > 0 )); then
        RRNA_23S="YES"
    else
        RRNA_23S="NO"
    fi


    echo
    echo "rRNA results:"
    echo "5S:  ${RRNA_5S} (${N_5S} copies)"
    echo "16S: ${RRNA_16S} (${N_16S} copies)"
    echo "23S: ${RRNA_23S} (${N_23S} copies)"


    # ========================================================
    # TRNASCAN-SE
    #
    # -G = general model containing tRNAs from all three
    # domains. This avoids assuming Bacteria vs Archaea.
    # ========================================================

    TRNA_OUT="${OUTDIR}/trnascan/${ID}.trnascan.tsv"
    TRNA_STATS="${OUTDIR}/trnascan/${ID}.stats"


    echo
    echo "Running tRNAscan-SE..."

    if ! tRNAscan-SE \
        -G \
        --thread "${SLURM_CPUS_PER_TASK}" \
        -o "${TRNA_OUT}" \
        -m "${TRNA_STATS}" \
        "${INPUT_FASTA}" \
        > "${OUTDIR}/logs/${ID}.trnascan.log" \
        2>&1
    then

        echo "ERROR: tRNAscan-SE failed"

        printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\tNA\tNA\tNA\tNA\tNA\tTRNASCAN_ERROR\n" \
            "${SAMPLE}" \
            "${GENOME}" \
            "${COMPLETENESS}" \
            "${CONTAMINATION}" \
            "${RRNA_5S}" \
            "${N_5S}" \
            "${RRNA_16S}" \
            "${N_16S}" \
            "${RRNA_23S}" \
            "${N_23S}" \
            >> "${SUMMARY}"

        [[ -n "${TEMP_FASTA}" ]] && rm -f "${TEMP_FASTA}"

        continue
    fi


    # ========================================================
    # COUNT TOTAL NON-PSEUDOGENE TRNAS
    # ========================================================

    TRNA_TOTAL=$(awk '

        $2 ~ /^[0-9]+$/ &&
        NF >= 5 &&
        tolower($0) !~ /pseudo/ {

            count++
        }

        END {
            print count + 0
        }

    ' "${TRNA_OUT}")


    # ========================================================
    # COUNT UNIQUE CANONICAL AMINO-ACID TYPES
    #
    # MIMAG requires tRNAs corresponding to >=18 of the
    # 20 canonical amino acids.
    # ========================================================

    TRNA_AA_TYPES=$(awk '

        BEGIN {

            split(
                "Ala Arg Asn Asp Cys Gln Glu Gly His Ile Leu Lys Met Phe Pro Ser Thr Trp Tyr Val",
                aa,
                " "
            )

            for (i = 1; i <= 20; i++) {
                valid[aa[i]] = 1
            }
        }


        $2 ~ /^[0-9]+$/ &&
        NF >= 5 &&
        tolower($0) !~ /pseudo/ {

            type = $5

            if (type in valid) {
                found[type] = 1
            }
        }


        END {

            count = 0

            for (type in found) {
                count++
            }

            print count
        }

    ' "${TRNA_OUT}")


    # ========================================================
    # MAKE LIST OF AA TYPES FOUND
    # ========================================================

    TRNA_AA_LIST=$(awk '

        BEGIN {

            split(
                "Ala Arg Asn Asp Cys Gln Glu Gly His Ile Leu Lys Met Phe Pro Ser Thr Trp Tyr Val",
                aa,
                " "
            )

            for (i = 1; i <= 20; i++) {
                valid[aa[i]] = 1
            }
        }


        $2 ~ /^[0-9]+$/ &&
        NF >= 5 &&
        tolower($0) !~ /pseudo/ {

            type = $5

            if (type in valid) {
                found[type] = 1
            }
        }


        END {

            first = 1

            for (i = 1; i <= 20; i++) {

                type = aa[i]

                if (type in found) {

                    if (!first) {
                        printf ","
                    }

                    printf "%s", type

                    first = 0
                }
            }

            printf "\n"
        }

    ' "${TRNA_OUT}")


    echo
    echo "tRNA results:"
    echo "Total predicted tRNAs: ${TRNA_TOTAL}"
    echo "Amino-acid types:      ${TRNA_AA_TYPES}/20"
    echo "Types found:           ${TRNA_AA_LIST}"


    # ========================================================
    # RNA PASS / FAIL
    # ========================================================

    if [[ "${RRNA_5S}" == "YES" ]] && \
       [[ "${RRNA_16S}" == "YES" ]] && \
       [[ "${RRNA_23S}" == "YES" ]] && \
       (( TRNA_AA_TYPES >= 18 ))
    then

        RNA_PASS="PASS"

    else

        RNA_PASS="FAIL"

    fi


    # ========================================================
    # FINAL MIMAG HIGH-QUALITY ASSESSMENT
    # ========================================================

    if [[ "${COMPLETENESS}" == "NA" ]] || \
       [[ "${CONTAMINATION}" == "NA" ]]
    then

        MIMAG_HQ="CHECKM_NOT_FOUND"

    else

        MIMAG_HQ=$(awk \
            -v comp="${COMPLETENESS}" \
            -v contam="${CONTAMINATION}" \
            -v rna="${RNA_PASS}" \
            'BEGIN {

                if (comp > 90 &&
                    contam < 5 &&
                    rna == "PASS") {

                    print "PASS"

                } else {

                    print "FAIL"
                }
            }'
        )

    fi


    # ========================================================
    # WRITE RESULT TO SUMMARY
    # ========================================================

    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\tOK\n" \
        "${SAMPLE}" \
        "${GENOME}" \
        "${COMPLETENESS}" \
        "${CONTAMINATION}" \
        "${RRNA_5S}" \
        "${N_5S}" \
        "${RRNA_16S}" \
        "${N_16S}" \
        "${RRNA_23S}" \
        "${N_23S}" \
        "${TRNA_TOTAL}" \
        "${TRNA_AA_TYPES}" \
        "${TRNA_AA_LIST}" \
        "${RNA_PASS}" \
        "${MIMAG_HQ}" \
        >> "${SUMMARY}"


    echo
    echo "------------------------------------------------------------"
    echo "RNA criterion: ${RNA_PASS}"
    echo "MIMAG HQ:      ${MIMAG_HQ}"
    echo "------------------------------------------------------------"


    # ========================================================
    # REMOVE TEMPORARY UNZIPPED FASTA
    # ========================================================

    if [[ -n "${TEMP_FASTA}" ]]; then
        rm -f "${TEMP_FASTA}"
    fi

done


# ============================================================
# FINISHED
# ============================================================

echo
echo "============================================================"
echo "ALL MAGS FINISHED"
echo "============================================================"
echo
echo "Final table:"
echo "${SUMMARY}"
echo
echo "Number passing MIMAG high-quality criteria:"
awk -F '\t' 'NR > 1 && $15 == "PASS" {n++} END {print n+0}' "${SUMMARY}"
echo
echo "Number failing MIMAG high-quality criteria:"
awk -F '\t' 'NR > 1 && $15 == "FAIL" {n++} END {print n+0}' "${SUMMARY}"
echo
echo "============================================================"