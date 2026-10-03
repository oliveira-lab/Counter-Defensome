#!/bin/bash


#########################################################################
# File Name: GO_prophages.sh
# Author(s): Lucas da Silva
# Institution: Genoscope, Evry, France
# Mail: ldasilva@genoscope.cns.fr
# Date: 25/09/2026
#########################################################################


set -uo pipefail

[[ "${1:-}" == "--vs2-db" ]] || { echo "ERROR: --vs2-db PATH is required (VirSorter2 database directory)"; exit 1; }
VS2_DB="${2:?--vs2-db requires a path}"; shift 2

INPUT_DIR="${1:?Usage: GO_prophages.sh --vs2-db PATH <INPUT_DIR> [RUN_LABEL] [OUTPUT_BASE]}"
INPUT_DIR="${INPUT_DIR%/}"
RUN_LABEL="${2:-$(basename "$INPUT_DIR")}"
OUTPUT_BASE="${3:-$(dirname "$INPUT_DIR")/virsorter2_results_${RUN_LABEL}}"
[[ -d "$INPUT_DIR" ]] || { echo "ERROR: INPUT_DIR not found: $INPUT_DIR"; exit 1; }

mkdir -p "$OUTPUT_BASE"
COMBINED_FASTA="$OUTPUT_BASE/all_bacterial_mags_prefixed.fna"
VS2_OUT="$OUTPUT_BASE/vs2_out"
HC_DIR="$OUTPUT_BASE/high_confidence"
mkdir -p "$HC_DIR"

#concatenate MAGs and prefix headers with source MAG_ID
> "$COMBINED_FASTA"
N_MAGS=0
for FASTA_FILE in "$INPUT_DIR"/*.{fasta,fa,fna}; do
    [ -e "$FASTA_FILE" ] || continue
    MAG_ID=$(basename "$FASTA_FILE"); MAG_ID="${MAG_ID%.*}"
    awk -v mag="$MAG_ID" '/^>/{sub(/^>/, ">" mag "__")} 1' "$FASTA_FILE" >> "$COMBINED_FASTA"
    N_MAGS=$((N_MAGS + 1))
done
[ "$N_MAGS" -gt 0 ] || { echo "ERROR: no .fasta/.fa/.fna files in $INPUT_DIR"; exit 1; }
echo "Concatenated $N_MAGS MAGs -> $COMBINED_FASTA"

#run VirSorter2
VS2_JOBS=$(( SLURM_CPUS_PER_TASK / 2 )); [ "$VS2_JOBS" -lt 1 ] && VS2_JOBS=1
ml virsorter
virsorter run \
    -w "$VS2_OUT" -i "$COMBINED_FASTA" -d "$VS2_DB" -j "$VS2_JOBS" -l "$RUN_LABEL" \
    --include-groups dsDNAphage,ssDNA --min-score 0.5 --min-length 0 --rm-tmpdir all

FINAL_SCORE_TSV="$VS2_OUT/${RUN_LABEL}-final-viral-score.tsv"
FINAL_VIRAL_FASTA="$VS2_OUT/${RUN_LABEL}-final-viral-combined.fa"
FINAL_BOUNDARY_TSV="$VS2_OUT/${RUN_LABEL}-final-viral-boundary.tsv"
[[ -f "$FINAL_SCORE_TSV" ]] || FINAL_SCORE_TSV="$VS2_OUT/final-viral-score.tsv"
[[ -f "$FINAL_VIRAL_FASTA" ]] || FINAL_VIRAL_FASTA="$VS2_OUT/final-viral-combined.fa"
[[ -f "$FINAL_BOUNDARY_TSV" ]] || FINAL_BOUNDARY_TSV="$VS2_OUT/final-viral-boundary.tsv"
[[ -f "$FINAL_SCORE_TSV" ]] || { echo "ERROR: VirSorter2 produced no final-viral-score.tsv"; exit 1; }

#high-confidence filter: max_score>=0.9 
HC_SCORE_TSV="$HC_DIR/high_confidence_score.tsv"
HC_IDS="$HC_DIR/high_confidence_ids.txt"
HC_FASTA="$HC_DIR/high_confidence.fna"

awk -F'\t' '
    NR==1 { for (i=1;i<=NF;i++) col[$i]=i; print; next }
    { s=$(col["max_score"])+0; if (s>=0.9) print }
' "$FINAL_SCORE_TSV" > "$HC_SCORE_TSV"
tail -n +2 "$HC_SCORE_TSV" | cut -f1 > "$HC_IDS"

python3 - "$HC_IDS" "$FINAL_VIRAL_FASTA" "$HC_FASTA" <<'PYEOF'
import sys
ids_file, fasta_in, fasta_out = sys.argv[1:4]
wanted = set(l.strip() for l in open(ids_file) if l.strip())
def read_fasta(p):
    h, seq = None, []
    for line in open(p):
        line = line.rstrip("\n")
        if line.startswith(">"):
            if h is not None: yield h, "".join(seq)
            h, seq = line[1:], []
        else:
            seq.append(line)
    if h is not None: yield h, "".join(seq)
n = 0
with open(fasta_out, "w") as out:
    for h, seq in read_fasta(fasta_in):
        if h.split()[0] in wanted:
            out.write(f">{h}\n{seq}\n"); n += 1
print(f"Extracted {n} high-confidence sequences -> {fasta_out}")
PYEOF

if [[ -f "$FINAL_BOUNDARY_TSV" ]]; then
    awk -F'\t' -v idfile="$HC_IDS" '
        BEGIN { while ((getline id < idfile) > 0) wanted[id]=1 }
        NR==1 {
            for (i=1;i<=NF;i++) col[$i]=i
            joincol = ("seqname_new" in col) ? col["seqname_new"] : col["seqname"]
            print; next
        }
        { if ($(joincol) in wanted) print }
    ' "$FINAL_BOUNDARY_TSV" > "$HC_DIR/high_confidence_boundary.tsv"
fi
