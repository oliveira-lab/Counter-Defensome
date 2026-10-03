#!/bin/bash

#########################################################################
# File Name: GO_tag_prophage_genes.sh
# Author(s): Lucas da Silva
# Institution: Genoscope, Evry, France
# Mail: ldasilva@genoscope.cns.fr
# Date: 25/09/2026
#########################################################################


set -uo pipefail

USAGE="Usage: bash GO_tag_prophage_genes.sh [--full-only|--partial-only] <HIGH_CONFIDENCE_SCORE_TSV> <HIGH_CONFIDENCE_BOUNDARY_TSV> <BACT_FFN_DIR> <OUT_TSV> [NEW_BIOGROUP] [MIN_LENGTH_BP]"

FULL_ONLY=0
PARTIAL_ONLY=0
POSITIONAL=()
for arg in "$@"; do
    case "$arg" in
        --full-only) FULL_ONLY=1 ;;
        --partial-only) PARTIAL_ONLY=1 ;;
        -h|--help) echo "$USAGE"; exit 0 ;;
        *) POSITIONAL+=("$arg") ;;
    esac
done
[[ "$FULL_ONLY" -eq 1 && "$PARTIAL_ONLY" -eq 1 ]] && { echo "ERROR: --full-only and --partial-only are mutually exclusive." >&2; exit 1; }
if [[ "$FULL_ONLY" -eq 1 ]]; then MODE="full"
elif [[ "$PARTIAL_ONLY" -eq 1 ]]; then MODE="partial"
else MODE="both"
fi
set -- "${POSITIONAL[@]}"

SCORE_TSV="${1:?$USAGE}"
BOUNDARY_TSV="${2:?$USAGE}"
BACT_FFN_DIR="${3:?$USAGE}"
OUT_TSV="${4:?$USAGE}"
NEW_BIOGROUP="${5:-prophage}"
MIN_LENGTH_BP="${6:-1000}"

[[ -s "$SCORE_TSV" ]] || { echo "ERROR: score tsv not found or empty: $SCORE_TSV" >&2; exit 1; }
[[ -d "$BACT_FFN_DIR" ]] || { echo "ERROR: BACT_FFN_DIR not found: $BACT_FFN_DIR" >&2; exit 1; }
if [[ "$MODE" != "full" && ( "$BOUNDARY_TSV" == "-" || -z "$BOUNDARY_TSV" || ! -s "$BOUNDARY_TSV" ) ]]; then
    echo "ERROR: HIGH_CONFIDENCE_BOUNDARY_TSV required for mode '$MODE' (run rebuild_high_confidence_boundary.sh, or use --full-only with '-')." >&2
    exit 1
fi
mkdir -p "$(dirname "$OUT_TSV")"

python3 - "$SCORE_TSV" "$BOUNDARY_TSV" "$BACT_FFN_DIR" "$OUT_TSV" "$NEW_BIOGROUP" "$MIN_LENGTH_BP" "$MODE" <<'PYEOF'
import sys, glob
from pathlib import Path
from collections import defaultdict

score_tsv, boundary_tsv, bact_ffn_dir, out_tsv, new_biogroup, min_length_bp_str, mode = sys.argv[1:8]
bact_ffn_dir = Path(bact_ffn_dir)
min_length_bp = float(min_length_bp_str)

TRIM_BP_START = "trim_bp_start"
TRIM_BP_END = "trim_bp_end"


def suffix_kind(suffix):
    s = suffix.lower()
    if s == "full":
        return "full"
    if "partial" in s:
        return "partial"
    return None


def bp_window(row):
    try:
        lo = int(float(row.get(TRIM_BP_START, "")))
        hi = int(float(row.get(TRIM_BP_END, "")))
        if lo > 0 and hi > 0:
            return min(lo, hi), max(lo, hi)
    except (TypeError, ValueError):
        pass
    return None


def load_boundary_map(path):
    if not path or path == "-" or not Path(path).is_file():
        return {}
    m = {}
    with open(path) as f:
        header = f.readline().rstrip("\n").split("\t")
        col = {name: i for i, name in enumerate(header)}
        joincol_name = "seqname_new" if "seqname_new" in col else ("seqname" if "seqname" in col else None)
        if joincol_name is None:
            return {}
        joincol = col[joincol_name]
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) <= joincol:
                continue
            m[parts[joincol]] = {name: (parts[i] if i < len(parts) else "") for name, i in col.items()}
    return m


boundary_map = load_boundary_map(boundary_tsv)

prophage_seqs = defaultdict(lambda: defaultdict(list))  # prophage_seqs[mag_id][seq_name] = [{"kind", "window"}, ...]
with open(score_tsv) as f:
    header = f.readline().rstrip("\n").split("\t")
    col_length = header.index("length")
    for line in f:
        parts = line.rstrip("\n").split("\t")
        if not parts or not parts[0]:
            continue
        full_seqname = parts[0]
        base, _, suffix = full_seqname.partition("||")
        kind = suffix_kind(suffix)
        if kind is None or (mode != "both" and kind != mode) or "__" not in base:
            continue
        try:
            length = float(parts[col_length])
        except (IndexError, ValueError):
            continue
        if length < min_length_bp:
            continue
        mag_id, seq_name = base.split("__", 1)
        entry = {"kind": kind, "window": None}
        if kind == "partial":
            row = boundary_map.get(full_seqname)
            if row is not None:
                entry["window"] = bp_window(row)
        prophage_seqs[mag_id][seq_name].append(entry)

ffn_files = sorted(set(
    glob.glob(str(bact_ffn_dir / "*.ffn")) + glob.glob(str(bact_ffn_dir / "*.fna")) +
    glob.glob(str(bact_ffn_dir / "*.fa")) + glob.glob(str(bact_ffn_dir / "*.fasta"))
))

tagged = []
for fp in ffn_files:
    fp = Path(fp)
    mag_id = fp.stem
    seqs = prophage_seqs.get(mag_id)
    if not seqs:
        continue
    seq_name_candidates = sorted(seqs, key=len, reverse=True)  # longest first, avoids prefix clashes
    for line in open(fp, errors="ignore"):
        if not line.startswith(">"):
            continue
        tokens = line[1:].split()
        if not tokens:
            continue
        seqid = tokens[0]
        seq_name = next((s for s in seq_name_candidates if seqid == s or seqid.startswith(s + "_")), None)
        if seq_name is None:
            continue
        entries = seqs[seq_name]
        if any(e["kind"] == "full" for e in entries):
            tagged.append(seqid)
            continue
        windows = [e["window"] for e in entries if e["kind"] == "partial" and e["window"] is not None]
        if not windows:
            continue
        try:
            gene_start = int(float(tokens[2]))
            gene_end = int(float(tokens[4]))
        except (IndexError, ValueError):
            continue
        gene_lo, gene_hi = min(gene_start, gene_end), max(gene_start, gene_end)
        if any(lo <= gene_lo and gene_hi <= hi for lo, hi in windows):
            tagged.append(seqid)

with open(out_tsv, "w") as f:
    f.write("gene_id\tbiogroup\n")
    for gene_id in tagged:
        f.write(f"{gene_id}\t{new_biogroup}\n")
PYEOF
