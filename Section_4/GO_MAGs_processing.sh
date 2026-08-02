#!/bin/bash

usage() {
    cat <<'USAGE'
Usage:
  checkm_v2-2_lineage_header.sh -i INPUT_DIR -f EXTENSION -t THREADS

Required arguments:
  -i  Directory containing the input genome FASTA files.
  -f  Input file extension without the leading dot (eg fna or fa).
  -t  Number of CPUs/threads to use.

USAGE
    exit 1
}

while getopts "i:f:t:h" opt; do
    case $opt in
        i) BASE_DIR="$OPTARG" ;;
        f) EXTENSION="$OPTARG" ;;
        t) THREADS="$OPTARG" ;;
        h) usage ;;
        *) usage ;;
    esac
done

if [[ -z "$BASE_DIR" || -z "$EXTENSION" || -z "$THREADS" ]]; then
    usage
fi

OUT_DIR="${BASE_DIR}/checkm_results"
SCRATCH_DIR="/tmp/checkm_tmp_$RANDOM"
CHECKM_FILE="${OUT_DIR}/checkm_summary_results.txt"
FINAL_DIR="${OUT_DIR}/selected_mags_10_90"

mkdir -p "$SCRATCH_DIR" "$OUT_DIR" "$FINAL_DIR"
export MPLCONFIGDIR="${SCRATCH_DIR}/matplotlib_cache"

ml checkm
checkm lineage_wf "$BASE_DIR" "$OUT_DIR" -x "$EXTENSION" -t "$THREADS" --tmpdir "$SCRATCH_DIR"
checkm qa "$OUT_DIR/lineage.ms" "$OUT_DIR" --out_format 2 --tab_table --file "$CHECKM_FILE"

awk -F '\t' '
NR == 1 {
    for (i = 1; i <= NF; i++) {
        gsub(/\r$/, "", $i)
        if ($i == "Bin Id") bin_col = i
        if ($i == "Completeness") completeness_col = i
        if ($i == "Contamination") contamination_col = i
    }
    next
}
$(completeness_col) > 90 && $(contamination_col) < 10 {
    print $(bin_col)
}
' "$CHECKM_FILE" | tr -d '\r' | xargs -L1 | while read -r bin_id; do
    for ext in fna fasta fa; do
        if [ -f "${BASE_DIR}/${bin_id}.${ext}" ]; then
            cp "${BASE_DIR}/${bin_id}.${ext}" "$FINAL_DIR/"
            break
        fi
    done
done

module unload checkm
ml drep

DREP_WORK="${OUT_DIR}/drep"
mapfile -t GENOME_FILES < <(find "$FINAL_DIR" -maxdepth 1 -name "*.${EXTENSION}" -o -name "*.fna" -o -name "*.fa" -o -name "*.fasta")

if [ "${#GENOME_FILES[@]}" -gt 0 ]; then
    dRep dereplicate "$DREP_WORK" -g "${GENOME_FILES[@]}" -p "$THREADS" --ignoreGenomeQuality -sa 0.95
fi

echo "Workflow finished $(date)"
