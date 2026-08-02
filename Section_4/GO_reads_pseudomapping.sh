#!/bin/bash

set -euo pipefail

usage() {
    cat <<'USAGE'
Usage:
  kallisto_for_defense_v4_peatland_flags.sh -p PREFIX -r READS_DIR -i KALLISTO_INDEX -t THREADS -m SORT_MEMORY -s STRANDNESS

Required arguments:
  -p  Prefix used to locate the non-rRNA reads and name the output files.
  -r  Reads directory containing mapping_PREFIX_defense/sortmerna/.
  -i  Kallisto index file.
  -t  Number of threads used by Kallisto and samtools.
  -m  Memory available to each samtools sort thread (for example: 48G).
  -s  Library strand specificity: unstranded, fr, or rf.

Strandness values:
  unstranded  Do not pass a strand-specificity flag to Kallisto.
  fr          Pass --fr-stranded to Kallisto.
  rf          Pass --rf-stranded to Kallisto.

USAGE
    exit 1
}

########################################
# ARGUMENTS
########################################

PREFIX=""
READS_DIR=""
KALLISTO_INDEX=""
THREADS=""
SORT_MEMORY=""
STRANDNESS=""

while getopts ":p:r:i:t:m:s:h" opt; do
    case "${opt}" in
        p) PREFIX="${OPTARG}" ;;
        r) READS_DIR="${OPTARG}" ;;
        i) KALLISTO_INDEX="${OPTARG}" ;;
        t) THREADS="${OPTARG}" ;;
        m) SORT_MEMORY="${OPTARG}" ;;
        s) STRANDNESS="${OPTARG}" ;;
        h) usage ;;
        :) echo "ERROR: option -${OPTARG} requires an argument." >&2; usage ;;
        \?) echo "ERROR: invalid option: -${OPTARG}" >&2; usage ;;
    esac
done

[[ -n "${PREFIX}" && -n "${READS_DIR}" && -n "${KALLISTO_INDEX}" && \
   -n "${THREADS}" && -n "${SORT_MEMORY}" && -n "${STRANDNESS}" ]] || usage

[[ "${THREADS}" =~ ^[1-9][0-9]*$ ]] || {
    echo "ERROR: THREADS must be a positive integer: ${THREADS}" >&2
    exit 1
}

case "${STRANDNESS}" in
    unstranded)
        STRAND_FLAG=()
        ;;
    fr)
        STRAND_FLAG=(--fr-stranded)
        ;;
    rf)
        STRAND_FLAG=(--rf-stranded)
        ;;
    *)
        echo "ERROR: STRANDNESS must be unstranded, fr, or rf: ${STRANDNESS}" >&2
        exit 1
        ;;
esac

########################################
# CONFIGURATION
########################################

SORTMERNA_DIR="${READS_DIR}/mapping_${PREFIX}_defense/sortmerna"
NONRRNA_R1="${SORTMERNA_DIR}/${PREFIX}_nonrrna_fwd.fq"
NONRRNA_R2="${SORTMERNA_DIR}/${PREFIX}_nonrrna_rev.fq"

KALLISTO_OUTDIR="${READS_DIR}/mapping_${PREFIX}_defense/kallisto_output"

RAW_PSEUDOBAM="${KALLISTO_OUTDIR}/pseudoalignments.bam"
SORTED_PSEUDOBAM="${KALLISTO_OUTDIR}/${PREFIX}.pseudobam.sorted.bam"
DEPTH_TSV="${KALLISTO_OUTDIR}/${PREFIX}.pseudobam.depth.tsv"
ABUNDANCE_TSV="${KALLISTO_OUTDIR}/abundance.tsv"

FORCE=0

########################################
# INPUT CHECKS
########################################

[[ -s "${NONRRNA_R1}" ]] || {
    echo "ERROR: R1 not found or empty: ${NONRRNA_R1}" >&2
    exit 1
}

[[ -s "${NONRRNA_R2}" ]] || {
    echo "ERROR: R2 not found or empty: ${NONRRNA_R2}" >&2
    exit 1
}

[[ -s "${KALLISTO_INDEX}" ]] || {
    echo "ERROR: Kallisto index not found or empty: ${KALLISTO_INDEX}" >&2
    exit 1
}

mkdir -p "${KALLISTO_OUTDIR}"

########################################
# KALLISTO + PSEUDOBAM
########################################

if [[ "${FORCE}" -eq 1 || ! -s "${ABUNDANCE_TSV}" || ! -s "${SORTED_PSEUDOBAM}" ]]; then
    echo "[$(date)] Running Kallisto..."

    ml conda
    conda activate kallisto
    ml samtools

    kallisto quant \
        -i "${KALLISTO_INDEX}" \
        -o "${KALLISTO_OUTDIR}" \
        -t "${THREADS}" \
        "${STRAND_FLAG[@]}" \
        --pseudobam \
        "${NONRRNA_R1}" \
        "${NONRRNA_R2}"

    [[ -s "${ABUNDANCE_TSV}" ]] || {
        echo "ERROR: Kallisto did not create abundance.tsv" >&2
        exit 1
    }

    [[ -s "${RAW_PSEUDOBAM}" ]] || {
        echo "ERROR: Kallisto did not create pseudoalignments.bam" >&2
        exit 1
    }

    samtools quickcheck -v "${RAW_PSEUDOBAM}"

    echo "[$(date)] Sorting pseudobam..."
    samtools sort \
        -@ "${THREADS}" \
        -m "${SORT_MEMORY}" \
        -o "${SORTED_PSEUDOBAM}" \
        "${RAW_PSEUDOBAM}"

    samtools quickcheck -v "${SORTED_PSEUDOBAM}"
    samtools index -@ "${THREADS}" "${SORTED_PSEUDOBAM}"

    rm -f "${RAW_PSEUDOBAM}"

    conda deactivate || true
else
    echo "[$(date)] Existing abundance.tsv and sorted pseudobam found; skipping Kallisto."
fi

########################################
# BASE-LEVEL COVERAGE
########################################

if [[ "${FORCE}" -eq 1 || ! -s "${DEPTH_TSV}" ]]; then
    echo "[$(date)] Calculating per-base depth..."

    ml samtools

    [[ -s "${SORTED_PSEUDOBAM}" ]] || {
        echo "ERROR: sorted pseudobam not found: ${SORTED_PSEUDOBAM}" >&2
        exit 1
    }

    samtools quickcheck -v "${SORTED_PSEUDOBAM}"

    samtools depth -aa -@ "${THREADS}" "${SORTED_PSEUDOBAM}" > "${DEPTH_TSV}.tmp"
    mv "${DEPTH_TSV}.tmp" "${DEPTH_TSV}"
else
    echo "[$(date)] Existing depth table found; skipping samtools depth."
fi

########################################
# FINAL VALIDATION
########################################

[[ -s "${ABUNDANCE_TSV}" ]] || {
    echo "ERROR: final abundance.tsv missing or empty" >&2
    exit 1
}

[[ -s "${DEPTH_TSV}" ]] || {
    echo "ERROR: final depth table missing or empty" >&2
    exit 1
}

echo "[$(date)] Finished successfully."
echo "Inputs ready for defense_expression_HKmedian:"
echo "  KALLISTO_TSV=${ABUNDANCE_TSV}"
echo "  DEPTH_TSV=${DEPTH_TSV}"
