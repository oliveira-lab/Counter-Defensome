#!/bin/bash

usage() {
    cat <<EOF
Usage:
  checkv_v4_global_prefix_nc85.sh -i INPUT_DIR -e CONDA_ENV_CHECKV -t CPUS -d CHECKV_DB -p PREFIX

Required arguments:
  -i  Directory containing the input FASTA files.
  -e  Conda environment containing CheckV.
  -t  Number of CPUs to use.
  -d  CheckV database directory.
  -p  One global prefix applied to every selected viral contig.
EOF
    exit 1
}

while getopts "i:e:t:d:p:" opt; do
    case $opt in
        i) INPUT_DIR="$OPTARG" ;;
        e) CONDA_ENV_CHECKV="$OPTARG" ;;
        t) CPUS_PER_TASK="$OPTARG" ;;
        d) CHECKVDB="$OPTARG" ;;
        p) PREFIX_SAMPLE="$OPTARG" ;;
        *) usage ;;
    esac
done

if [[ -z "$INPUT_DIR" || -z "$CONDA_ENV_CHECKV" || -z "$CPUS_PER_TASK" || -z "$CHECKVDB" || -z "$PREFIX_SAMPLE" ]]; then
    usage
fi

PREFIX_SAMPLE="${PREFIX_SAMPLE%_}"

OUTPUT_BASE="${INPUT_DIR}/checkv_result"
mkdir -p "$OUTPUT_BASE"

find "$INPUT_DIR" -maxdepth 1 -type f \( -name "*.fasta" -o -name "*.fna" -o -name "*.fa" \) | while read -r FASTA_FILE; do
    SAMPLE_ID=$(basename "${FASTA_FILE%.*}")
    TARGET_DIR="${OUTPUT_BASE}/${SAMPLE_ID}_checkv"

    module unload conda 2>/dev/null || true
    ml conda
    conda run --name "$CONDA_ENV_CHECKV" checkv end_to_end \
        "$FASTA_FILE" "$TARGET_DIR" -t "$CPUS_PER_TASK" -d "$CHECKVDB" --remove_tmp

    INPUT_TSV="$TARGET_DIR/quality_summary.tsv"
    INPUT_FNA="$TARGET_DIR/viruses.fna"
    OUTPUT_DIR="$TARGET_DIR/checkv_filtered_results"

    if [[ -f "$INPUT_TSV" ]]; then
        mkdir -p "$OUTPUT_DIR"

        awk -F'\t' 'NR==1 || $8 == "High-quality" || $8 == "Complete"' \
            "$INPUT_TSV" > "$OUTPUT_DIR/selected_hits.tsv"

        mapfile -t IDs < <(
            awk -F'\t' 'NR>1 && ($8 == "High-quality" || $8 == "Complete") {print $1}' \
                "$INPUT_TSV"
        )

        if [[ ${#IDs[@]} -gt 0 ]]; then
            for id in "${IDs[@]}"; do
                NEW_ID="${PREFIX_SAMPLE}_${SAMPLE_ID}_${id}"
                TMP_FILE="$OUTPUT_DIR/${NEW_ID}.fna.tmp"
                OUTPUT_FILE="$OUTPUT_DIR/${NEW_ID}.fna"

                awk -v id="$id" -v newid="$NEW_ID" '
                    BEGIN {RS=">"; FS="\n"}
                    {
                        split($1, h, /[[:space:]]+/)
                        if (h[1] == id) {
                            print ">" newid
                            for (i=2; i<=NF; i++) {
                                if ($i != "") print $i
                            }
                        }
                    }
                ' "$INPUT_FNA" > "$TMP_FILE"

                if [[ -s "$TMP_FILE" ]]; then
                    mv "$TMP_FILE" "$OUTPUT_FILE"
                else
                    rm -f "$TMP_FILE"
                    echo "[WARNING] ID not found in viruses.fna: $id" >&2
                fi
            done
        fi
    fi
done

module unload conda 2>/dev/null || true
ml drep

DREP_OUTPUT_DIR="${OUTPUT_BASE}/drep_results"

mapfile -t FILES < <(find "${OUTPUT_BASE}" -type f \( -name "*.fna" -o -name "*.fasta" -o -name "*.fa" \))

if [ ${#FILES[@]} -gt 0 ]; then
    dRep dereplicate "${DREP_OUTPUT_DIR}" \
        -g "${FILES[@]}" \
        --processors "$CPUS_PER_TASK" -sa 0.95 -nc 0.85 -l 18000 --ignoreGenomeQuality
else
    echo "No genome file found in ${OUTPUT_BASE} to dereplicate."
fi

echo "Workflow finished $(date)"