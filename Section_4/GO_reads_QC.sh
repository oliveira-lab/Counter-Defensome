#!/bin/bash
set -euo pipefail
shopt -s nullglob

usage() {
    cat <<'USAGE'
Usage:
  treat_reads-2.sh -i INPUT_DIR -t THREADS -m MEMORY_MB -d SORTMERNA_DB

Required arguments:
  -i  Directory containing the paired-end FASTQ files.
  -t  Number of threads to use.
  -m  Memory available to SortMeRNA, in MB.
  -d  SortMeRNA rRNA database path.

USAGE
    exit 1
}

while getopts "i:t:m:d:" opt; do
    case $opt in
        i) READS_DIR="$OPTARG" ;;
        t) THREADS="$OPTARG" ;;
        m) MEMORY="$OPTARG" ;;
        d) SMR_DB="$OPTARG" ;;
        *) usage ;;
    esac
done

if [[ -z "${READS_DIR:-}" || -z "${THREADS:-}" || -z "${MEMORY:-}" || -z "${SMR_DB:-}" ]]; then
    usage
fi

OUTDIR="${READS_DIR}/mapping_defense"
PAIR_DIR="${OUTDIR}/per_pair_processing"
FASTP_DIR="${PAIR_DIR}/fastp"
SORTMERNA_DIR="${PAIR_DIR}/sortmerna"
REPAIRED_DIR="${PAIR_DIR}/paired_nonrrna_repaired"
MANIFEST_DIR="${OUTDIR}/manifests"

FINAL_R1="${OUTDIR}/clean_nonrrna_R1.fastq.gz"
FINAL_R2="${OUTDIR}/clean_nonrrna_R2.fastq.gz"
PROCESSING_MANIFEST="${MANIFEST_DIR}/processing_manifest.tsv"

mkdir -p "${OUTDIR}" "${FASTP_DIR}" "${SORTMERNA_DIR}" "${REPAIRED_DIR}" "${MANIFEST_DIR}"

echo -e "sample_id\tinput_R1\tinput_R2\tfastp_R1\tfastp_R2\tnonrrna_R1\tnonrrna_R2" > "${PROCESSING_MANIFEST}"

cd "${READS_DIR}"
found=0

for r1 in *_1_*.fastq *_1_*.fq *_1_*.fastq.gz *_1_*.fq.gz *_R1_*.fastq *_R1_*.fq *_R1_*.fastq.gz *_R1_*.fq.gz; do
    [[ -e "$r1" ]] || continue
    [[ "$r1" == *"_single"* ]] && continue

    if [[ "$r1" == *"_1_"* ]]; then
        r2="${r1/_1_/_2_}"
        sample_id="${r1/_1_/_}"
    else
        r2="${r1/_R1_/_R2_}"
        sample_id="${r1/_R1_/_}"
    fi

    [[ -f "$r2" ]] || continue
    found=1

    sample_id="${sample_id%.fastq.gz}"
    sample_id="${sample_id%.fq.gz}"
    sample_id="${sample_id%.fastq}"
    sample_id="${sample_id%.fq}"
    safe_id=$(echo "${sample_id}" | sed 's#[/ ]#_#g')

    FASTP_R1="${FASTP_DIR}/${safe_id}_R1.fastp.fastq.gz"
    FASTP_R2="${FASTP_DIR}/${safe_id}_R2.fastp.fastq.gz"

    SMR_SAMPLE_DIR="${SORTMERNA_DIR}/${safe_id}"
    SMR_WORK="${SMR_SAMPLE_DIR}/work"
    mkdir -p "${SMR_SAMPLE_DIR}" "${SMR_WORK}"

    RAW_NONRRNA_R1="${SMR_SAMPLE_DIR}/${safe_id}_nonrrna_fwd.fq.gz"
    RAW_NONRRNA_R2="${SMR_SAMPLE_DIR}/${safe_id}_nonrrna_rev.fq.gz"

    REPAIRED_R1="${REPAIRED_DIR}/${safe_id}_nonrrna_paired_R1.fastq.gz"
    REPAIRED_R2="${REPAIRED_DIR}/${safe_id}_nonrrna_paired_R2.fastq.gz"
    SINGLETONS="${REPAIRED_DIR}/${safe_id}_nonrrna_singletons.fastq.gz"

    ml fastp
    fastp -i "${READS_DIR}/${r1}" -I "${READS_DIR}/${r2}" -o "${FASTP_R1}" -O "${FASTP_R2}" \
          -q 25 -e 25 -u 20 --length_required 50 --detect_adapter_for_pe --thread "${THREADS}" \
          --json "${FASTP_DIR}/${safe_id}.fastp.json" --html "${FASTP_DIR}/${safe_id}.fastp.html"
    module unload fastp

    ml sortmerna
    sortmerna --ref "${SMR_DB}" --reads "${FASTP_R1}" --reads "${FASTP_R2}" \
              --workdir "${SMR_WORK}" --threads "${THREADS}" -m "${MEMORY}" \
              --fastx --num_alignments 1 --max_pos 20 --paired_in --out2 \
              --other "${SMR_SAMPLE_DIR}/${safe_id}_nonrrna" \
              --aligned "${SMR_SAMPLE_DIR}/${safe_id}_rrna"
    module unload sortmerna/4.3.4

    ml bbtools
    repair.sh in1="${RAW_NONRRNA_R1}" in2="${RAW_NONRRNA_R2}" \
              out1="${REPAIRED_R1}" out2="${REPAIRED_R2}" outs="${SINGLETONS}"
    module unload bbtools

    echo -e "${safe_id}\t${READS_DIR}/${r1}\t${READS_DIR}/${r2}\t${FASTP_R1}\t${FASTP_R2}\t${REPAIRED_R1}\t${REPAIRED_R2}" >> "${PROCESSING_MANIFEST}"
done

if [[ "${found}" -ne 1 ]]; then
    echo "ERROR: No FASTQ pair detected in ${READS_DIR}" >&2
    exit 1
fi

rm -f "${FINAL_R1}" "${FINAL_R2}"
tail -n +2 "${PROCESSING_MANIFEST}" | while IFS=$'\t' read -r id in1 in2 f1 f2 rep1 rep2; do
    cat "${rep1}" >> "${FINAL_R1}"
    cat "${rep2}" >> "${FINAL_R2}"
done

n1=$(zcat "${FINAL_R1}" | awk 'END{print NR/4}')
echo -e "FINAL_R1\tFINAL_R2\tN_READ_PAIRS\n${FINAL_R1}\t${FINAL_R2}\t${n1}" > "${OUTDIR}/reads_treatment_summary.tsv"

echo "Workflow finished $(date)"
