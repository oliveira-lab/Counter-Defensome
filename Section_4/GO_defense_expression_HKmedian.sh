#!/bin/bash
set -euo pipefail
shopt -s nullglob

usage() {
    cat <<'USAGE'
Usage:
  defense_expression_HKmedian_v5_flags.sh -p PREFIX -r READS_DIR -d DEFENSE_DIR -v VIRAL_FFN_DIR -b BACT_FFN_DIR -f DF_REF -a ADF_REF -k HOUSEKEEPING_GENE_LIST -n MIN_READS -c MIN_COVERAGE

Required arguments:
  -p  Dataset prefix.
  -r  Directory containing the mapping outputs.
  -d  Directory containing DefenseFinder outputs.
  -v  Directory containing viral CDS FASTA files.
  -b  Directory containing bacterial CDS FASTA files.
  -f  Defense HMM reference list.
  -a  Counter-defense HMM reference list.
  -k  Housekeeping gene-name list.
  -n  Minimum mapped-read threshold.
  -c  Minimum covered-fraction threshold (eg. 0.5 for 50%).

USAGE
    exit 1
}

PREFIX=""
READS_DIR=""
DEFENSE_DIR=""
VIRAL_FFN_DIR=""
BACT_FFN_DIR=""
DF_REF=""
ADF_REF=""
HOUSEKEEPING_GENE_LIST=""
MIN_READS=""
MIN_COVERAGE=""

while getopts ":p:r:d:v:b:f:a:k:n:c:h" opt; do
    case "${opt}" in
        p) PREFIX="${OPTARG}" ;;
        r) READS_DIR="${OPTARG}" ;;
        d) DEFENSE_DIR="${OPTARG}" ;;
        v) VIRAL_FFN_DIR="${OPTARG}" ;;
        b) BACT_FFN_DIR="${OPTARG}" ;;
        f) DF_REF="${OPTARG}" ;;
        a) ADF_REF="${OPTARG}" ;;
        k) HOUSEKEEPING_GENE_LIST="${OPTARG}" ;;
        n) MIN_READS="${OPTARG}" ;;
        c) MIN_COVERAGE="${OPTARG}" ;;
        h) usage ;;
        *) usage ;;
    esac
done

[[ -n "${PREFIX}" ]] || usage
[[ -n "${READS_DIR}" ]] || usage
[[ -n "${DEFENSE_DIR}" ]] || usage
[[ -n "${VIRAL_FFN_DIR}" ]] || usage
[[ -n "${BACT_FFN_DIR}" ]] || usage
[[ -n "${DF_REF}" ]] || usage
[[ -n "${ADF_REF}" ]] || usage
[[ -n "${HOUSEKEEPING_GENE_LIST}" ]] || usage
[[ -n "${MIN_READS}" ]] || usage
[[ -n "${MIN_COVERAGE}" ]] || usage

########################################
# CONFIGURATION
########################################

# Kallisto/SAMtools outputs
KALLISTO_OUTDIR="${READS_DIR}/mapping_${PREFIX}_defense/kallisto_output"
KALLISTO_TSV="${KALLISTO_OUTDIR}/abundance.tsv"
DEPTH_TSV="${KALLISTO_OUTDIR}/${PREFIX}.pseudobam.depth.tsv"

# Prokka/Prodigal per-genome annotation outputs for the bacterial MAGs
BACT_PROKKA_ANNOT="${BACT_FFN_DIR%/ffn_cds}/prokka_prodigal_link/tables_per_genome/prodigal_to_prokka_annotation.tsv"

OUTDIR="${READS_DIR}/expression_${PREFIX}_defense_systems_hmmer_v7_HKmedian_ratio"

mkdir -p "${OUTDIR}"

########################################
# INPUT CHECKS
########################################

[[ -s "${KALLISTO_TSV}" ]] || { echo "ERROR: KALLISTO_TSV not found: ${KALLISTO_TSV}" >&2; exit 1; }
[[ -s "${DEPTH_TSV}" ]] || { echo "ERROR: DEPTH_TSV not found: ${DEPTH_TSV}" >&2; exit 1; }
[[ -d "${DEFENSE_DIR}" ]] || { echo "ERROR: DEFENSE_DIR not found: ${DEFENSE_DIR}" >&2; exit 1; }
[[ -d "${VIRAL_FFN_DIR}" ]] || { echo "ERROR: VIRAL_FFN_DIR not found: ${VIRAL_FFN_DIR}" >&2; exit 1; }
[[ -d "${BACT_FFN_DIR}" ]] || { echo "ERROR: BACT_FFN_DIR not found: ${BACT_FFN_DIR}" >&2; exit 1; }
[[ -s "${DF_REF}" ]] || { echo "ERROR: DF_REF not found: ${DF_REF}" >&2; exit 1; }
[[ -s "${ADF_REF}" ]] || { echo "ERROR: ADF_REF not found: ${ADF_REF}" >&2; exit 1; }
[[ -s "${BACT_PROKKA_ANNOT}" ]] || { echo "ERROR: BACT_PROKKA_ANNOT not found: ${BACT_PROKKA_ANNOT}" >&2; exit 1; }
[[ -s "${HOUSEKEEPING_GENE_LIST}" ]] || { echo "ERROR: HOUSEKEEPING_GENE_LIST not found: ${HOUSEKEEPING_GENE_LIST}" >&2; exit 1; }

PY="${OUTDIR}/build_HKmedian_ratio_reports_v2.py"

cat > "${PY}" <<'PYCODE'
#!/usr/bin/env python3

import argparse
import glob
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")

SYSTEM_COLS = [
    "gene_id", "raw_hit_id", "replicon", "gene_name", "sys_id",
    "defense_system", "subtype", "activity", "source_level"
]
HMM_COLS = SYSTEM_COLS + ["hit_score"]


def clean(x):
    return str(x).replace("/", "_").replace(" ", "_").replace(":", "_")


def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument("--prefix", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--kallisto", required=True)
    ap.add_argument("--depth", required=True)
    ap.add_argument("--defense-dir", required=True)
    ap.add_argument("--viral-ffn", required=True)
    ap.add_argument("--bact-ffn", required=True)
    ap.add_argument("--df-ref", required=True)
    ap.add_argument("--adf-ref", required=True)
    ap.add_argument("--bact-prokka-annot", required=True)
    ap.add_argument("--housekeeping-gene-list", required=True)
    ap.add_argument("--min-reads", type=float, default=5)
    ap.add_argument("--min-coverage", type=float, default=0.50)
    return ap.parse_args()


def read_hmm_ref(path, activity):
    rows = []
    with open(path) as f:
        for line in f:
            if not line.strip():
                continue
            p = line.strip().split()
            if len(p) < 3:
                continue
            rows.append({
                "defense_system": p[0],
                "subtype": p[1],
                "gene_name": p[2],
                "activity": activity,
            })
    return pd.DataFrame(rows, columns=["defense_system", "subtype", "gene_name", "activity"])


def make_gene_biogroup_map(viral_dir, bact_dir):
    rows = []
    for bg, d in [("viral", viral_dir), ("bacterial", bact_dir)]:
        for ext in ("*.ffn", "*.fna", "*.fa", "*.fasta"):
            for fp in sorted(glob.glob(str(Path(d) / "**" / ext), recursive=True)):
                base = Path(fp).name
                with open(fp, errors="ignore") as f:
                    for line in f:
                        if line.startswith(">"):
                            gid = line[1:].strip().split()[0]
                            rows.append((gid, bg, base))
    return pd.DataFrame(rows, columns=["gene_id", "biogroup", "ffn_file"])


def read_coverage(depth_tsv):
    cov = {}
    total = {}
    with open(depth_tsv) as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) < 3:
                continue
            gid = p[0]
            dp = int(float(p[2]))
            total[gid] = total.get(gid, 0) + 1
            if dp > 0:
                cov[gid] = cov.get(gid, 0) + 1

    rows = []
    for gid, npos in total.items():
        covered = cov.get(gid, 0)
        rows.append({
            "gene_id": gid,
            "covered_bases": covered,
            "positions_in_pseudobam": npos,
            "covered_fraction": covered / npos if npos > 0 else 0,
        })
    return pd.DataFrame(rows)


def add_rpk(df):
    """Calculate reads per kilobase using Kallisto effective length."""
    df = df.copy()
    denom = df["eff_length"].replace(0, np.nan) / 1000.0
    df["rpk"] = df["mapped_reads"] / denom
    df["rpk"] = df["rpk"].replace([np.inf, -np.inf], np.nan).fillna(0)
    return df


def read_housekeeping_gene_names(path):
    names = []
    with open(path) as f:
        for line in f:
            x = line.strip()
            if x and not x.startswith("#"):
                names.append(x)
    return set(names)


def read_housekeeping_annotation(path, hk_gene_names):
    annot = pd.read_csv(path, sep="\t")

    required = {"prodigal_gene_id", "gene"}
    missing = required - set(annot.columns)
    if missing:
        raise ValueError(f"Housekeeping annotation missing columns: {missing}")

    annot["gene_clean"] = (
        annot["gene"]
        .astype(str)
        .str.strip()
        .str.replace(r"_\d+$", "", regex=True)
    )

    hk = annot[annot["gene_clean"].isin(hk_gene_names)].copy()

    hk = hk.rename(columns={"prodigal_gene_id": "gene_id"})
    hk = hk.drop_duplicates("gene_id")

    return hk


def compute_hk_median_baseline(expr_support, hk_annot):
    """
    Calculate one community-level housekeeping baseline:

        median RPK of bacterial housekeeping gene records

    Eligible housekeeping records must already pass:
      - mapped_reads >= MIN_READS
      - covered_fraction >= MIN_COVERAGE

    Supported housekeeping records necessarily have RPK > 0.
    """
    hk = expr_support.merge(
        hk_annot[["gene_id", "gene_clean"]].drop_duplicates("gene_id"),
        on="gene_id",
        how="inner",
    )
    hk = hk[hk["biogroup"] == "bacterial"].copy()

    if hk.empty:
        raise ValueError(
            "No bacterial housekeeping genes passed the reads and coverage filters; "
            "cannot calculate the median housekeeping RPK baseline."
        )

    baseline = float(hk["rpk"].median())

    if not np.isfinite(baseline) or baseline <= 0:
        raise ValueError(f"Invalid median housekeeping RPK baseline: {baseline}")

    hk_records = hk.sort_values(["gene_clean", "rpk"], ascending=[True, False])

    hk_summary = pd.DataFrame([{
        "normalization": "median_RPK_ratio",
        "housekeeping_baseline_rpk": baseline,
        "formula": "RPK_gene / median_RPK_housekeeping",
        "scale_factor": 1,
        "n_housekeeping_gene_ids_in_annotation": hk_annot["gene_id"].nunique(),
        "n_housekeeping_gene_ids_passing_support_filters": hk["gene_id"].nunique(),
        "n_housekeeping_gene_names_passing_support_filters": hk["gene_clean"].nunique(),
    }])

    return baseline, hk_records, hk_summary

def make_master_expression(kallisto_tsv, gene_map, coverage_df, min_reads, min_coverage, hk_annot):
    ab = pd.read_csv(kallisto_tsv, sep="\t")
    ab = ab.rename(columns={
        "target_id": "gene_id",
        "est_counts": "mapped_reads",
        "length": "length_bp",
        "eff_length": "eff_length",
    })

    for c in ["mapped_reads", "length_bp", "eff_length"]:
        ab[c] = pd.to_numeric(ab[c], errors="coerce").fillna(0)
    gm = gene_map.drop_duplicates("gene_id").copy()
    m = ab.merge(gm, on="gene_id", how="left")
    m["biogroup"] = m["biogroup"].fillna("unknown")
    m["ffn_file"] = m["ffn_file"].fillna("unknown")

    m = m.merge(coverage_df, on="gene_id", how="left")
    m["covered_bases"] = m["covered_bases"].fillna(0).astype(int)
    m["positions_in_pseudobam"] = m["positions_in_pseudobam"].fillna(0).astype(int)
    m["covered_fraction"] = m["covered_fraction"].fillna(0)

    raw = add_rpk(m)

    support = raw[
        (raw["mapped_reads"] >= min_reads) &
        (raw["covered_fraction"] >= min_coverage)
    ].copy()

    baseline, hk_records, hk_summary = compute_hk_median_baseline(
        support, hk_annot
    )
    hk_summary["support_filter_min_reads"] = min_reads
    hk_summary["support_filter_min_coverage"] = min_coverage

    # housekeeping RPK
    raw["HK_norm_RPK_median"] = np.where(
        baseline > 0, raw["rpk"] / baseline, 0
    )
    support["HK_norm_RPK_median"] = np.where(
        baseline > 0, support["rpk"] / baseline, 0
    )

    cols = [
        "gene_id", "mapped_reads", "length_bp", "eff_length",
        "covered_bases", "positions_in_pseudobam", "covered_fraction",
        "rpk",
        "HK_norm_RPK_median",
        "biogroup", "ffn_file",
    ]

    # 'support' is the final analysis table after reads and coverage filters.
    return raw[cols], support[cols], baseline, hk_records, hk_summary


def parse_systems(defense_dir):
    system_files = glob.glob(str(Path(defense_dir) / "**" / "*_defense_finder_systems.tsv"), recursive=True)
    all_rows = []
    for sys_fp in system_files:
        sys_fp = Path(sys_fp)
        genes_fp = Path(str(sys_fp).replace("_defense_finder_systems.tsv", "_defense_finder_genes.tsv"))
        if not genes_fp.exists():
            raise FileNotFoundError(f"Missing genes file for {sys_fp}: {genes_fp}")

        genes = pd.read_csv(genes_fp, sep="\t")
        needed_genes = {"replicon", "hit_id", "gene_name", "sys_id", "type", "subtype", "activity"}
        if not needed_genes.issubset(genes.columns):
            raise ValueError(f"Missing columns in {genes_fp}: {needed_genes - set(genes.columns)}")

        for _, r in genes.iterrows():
            gid = str(r["hit_id"]).strip()
            if not gid or gid == "nan":
                continue
            all_rows.append({
                "gene_id": gid,
                "raw_hit_id": gid,
                "replicon": str(r["replicon"]),
                "gene_name": str(r["gene_name"]),
                "sys_id": str(r["sys_id"]),
                "defense_system": str(r["type"]),
                "subtype": str(r["subtype"]),
                "activity": str(r["activity"]),
                "source_level": "systems",
            })
    return pd.DataFrame(all_rows, columns=SYSTEM_COLS)


def build_validated_hmm_reference(df_ref_path, adf_ref_path):
    df_ref = read_hmm_ref(df_ref_path, "Defense")
    adf_ref = read_hmm_ref(adf_ref_path, "Antidefense")

    collisions = sorted(
        set(df_ref["gene_name"]) & set(adf_ref["gene_name"])
    )

    if collisions:
        raise ValueError(
            "The same gene_name occurs in both DF_REF and ADF_REF. "
            "These genes cannot be assigned unambiguously:\n"
            + "\n".join(collisions)
        )

    return pd.concat([df_ref, adf_ref], ignore_index=True)


def parse_hmmer(defense_dir, ref_df):
    files = glob.glob(str(Path(defense_dir) / "**" / "*_defense_finder_hmmer.tsv"), recursive=True)
    ref = ref_df.drop_duplicates("gene_name").copy()
    ref_map = ref.set_index("gene_name")[["defense_system", "subtype", "activity"]].to_dict("index")
    best = {}

    for fp in files:
        try:
            df = pd.read_csv(fp, sep="\t")
        except Exception:
            continue
        needed = {"hit_id", "replicon", "gene_name", "hit_score"}
        if not needed.issubset(df.columns):
            continue
        df["hit_score"] = pd.to_numeric(df["hit_score"], errors="coerce").fillna(-np.inf)
        for _, r in df.iterrows():
            gname = str(r["gene_name"])
            if gname not in ref_map:
                continue
            gid = str(r["hit_id"])
            score = float(r["hit_score"])
            if gid not in best or score > best[gid]["hit_score"]:
                info = ref_map[gname]
                best[gid] = {
                    "gene_id": gid,
                    "raw_hit_id": gid,
                    "replicon": str(r["replicon"]),
                    "gene_name": gname,
                    "sys_id": str(r["replicon"]) + "_" + info["defense_system"],
                    "defense_system": info["defense_system"],
                    "subtype": info["subtype"],
                    "activity": info["activity"],
                    "source_level": "hmmer",
                    "hit_score": score,
                }
    return pd.DataFrame(best.values(), columns=HMM_COLS)


def metric_col():
    return "HK_norm_RPK_median"

def summarize(df, group_cols):
    if df.empty:
        return pd.DataFrame()
    mcol = metric_col()
    return (
        df.groupby(group_cols, as_index=False)
        .agg(
            n_genes=("gene_id", "nunique"),
            supported_genes=("gene_id", "nunique"),
            mapped_reads=("mapped_reads", "sum"),
            sum_HK_norm_expression=(mcol, "sum"),
            mean_HK_norm_expression=(mcol, "mean"),
            median_HK_norm_expression=(mcol, "median"),
            max_HK_norm_expression=(mcol, "max"),
            mean_covered_fraction=("covered_fraction", "mean"),
        )
        .sort_values("sum_HK_norm_expression", ascending=False)
    )


def write_standard_outputs(df, out, prefix):
    out.mkdir(parents=True, exist_ok=True)
    df.to_csv(out / f"{prefix}_defense_gene_HKmedian_ratio.tsv", sep="\t", index=False)
    summarize(df, ["subtype"]).to_csv(out / f"{prefix}_defense_subtype_HKmedian_ratio.tsv", sep="\t", index=False)
    summarize(df, ["defense_system"]).to_csv(out / f"{prefix}_defense_system_HKmedian_ratio.tsv", sep="\t", index=False)
    summarize(df, ["biogroup"]).to_csv(out / f"{prefix}_defense_biogroup_HKmedian_ratio.tsv", sep="\t", index=False)
    summarize(df, ["biogroup", "defense_system"]).to_csv(out / f"{prefix}_defense_system_biogroup_HKmedian_ratio.tsv", sep="\t", index=False)
    summarize(df, ["biogroup", "subtype"]).to_csv(out / f"{prefix}_defense_subtype_biogroup_HKmedian_ratio.tsv", sep="\t", index=False)
    summarize(df, ["biogroup", "replicon"]).to_csv(out / f"{prefix}_defense_replicon_HKmedian_ratio.tsv", sep="\t", index=False)
    summarize(df, ["biogroup", "activity"]).to_csv(out / f"{prefix}_defense_activity_biogroup_HKmedian_ratio.tsv", sep="\t", index=False)
    if "gene_name" in df.columns:
        summarize(df, ["gene_name"]).to_csv(out / f"{prefix}_defense_hmmer_signature_HKmedian_ratio.tsv", sep="\t", index=False)
        summarize(df, ["biogroup", "gene_name"]).to_csv(out / f"{prefix}_defense_hmmer_signature_biogroup_HKmedian_ratio.tsv", sep="\t", index=False)


def make_advanced_reports(df, outdir, prefix, label):
    d = outdir / clean(label)
    d.mkdir(parents=True, exist_ok=True)
    if df.empty:
        (d / "NO_RECORDS.txt").write_text(f"No records for {label}.\n")
        return
    df.to_csv(d / f"{prefix}_{clean(label)}_gene_HKmedian_ratio.tsv", sep="\t", index=False)
    summarize(df, ["biogroup"]).to_csv(d / f"{prefix}_{clean(label)}_by_biogroup.tsv", sep="\t", index=False)
    summarize(df, ["activity"]).to_csv(
    d / f"{prefix}_{clean(label)}_by_activity.tsv",
    sep="\t",
    index=False)
    summarize(df, ["defense_system"]).to_csv(d / f"{prefix}_{clean(label)}_by_defense_system.tsv", sep="\t", index=False)
    summarize(df, ["subtype"]).to_csv(d / f"{prefix}_{clean(label)}_by_subtype.tsv", sep="\t", index=False)
    summarize(df, ["biogroup", "defense_system"]).to_csv(d / f"{prefix}_{clean(label)}_by_biogroup_system.tsv", sep="\t", index=False)
    summarize(df, ["biogroup", "subtype"]).to_csv(d / f"{prefix}_{clean(label)}_by_biogroup_subtype.tsv", sep="\t", index=False)
    summarize(df, ["biogroup", "replicon"]).to_csv(d / f"{prefix}_{clean(label)}_by_biogroup_replicon.tsv", sep="\t", index=False)


def write_split_outputs(df, base_out, prefix):
    split_root = base_out / "by_biogroup"
    split_root.mkdir(parents=True, exist_ok=True)
    for bg in ["viral", "bacterial"]:
        sub = df[df["biogroup"] == bg].copy()
        bg_out = split_root / bg
        bg_adv = bg_out / "advanced_reports"
        bg_out.mkdir(parents=True, exist_ok=True)
        bg_adv.mkdir(parents=True, exist_ok=True)
        if sub.empty:
            (bg_out / f"{prefix}_{bg}_NO_RECORDS.txt").write_text(f"No {bg} records found.\n")
            continue
        write_standard_outputs(sub, bg_out, f"{prefix}_{bg}")
        make_advanced_reports(sub, bg_adv, f"{prefix}_{bg}", "all_activity")
        make_advanced_reports(sub[sub["activity"] == "Defense"].copy(), bg_adv, f"{prefix}_{bg}", "Defense")
        make_advanced_reports(sub[sub["activity"] == "Antidefense"].copy(), bg_adv, f"{prefix}_{bg}", "Antidefense")


def write_reports(df, outdir, prefix):
    all_out = outdir / "all_biogroups"
    adv = all_out / "advanced_reports"
    all_out.mkdir(parents=True, exist_ok=True)
    adv.mkdir(parents=True, exist_ok=True)
    write_standard_outputs(df, all_out, prefix)
    make_advanced_reports(df, adv, prefix, "all_activity")
    make_advanced_reports(df[df["activity"] == "Defense"].copy(), adv, prefix, "Defense")
    make_advanced_reports(df[df["activity"] == "Antidefense"].copy(), adv, prefix, "Antidefense")
    write_split_outputs(df, outdir, prefix)


def technical_controls(raw, support, outdir, prefix, baseline, min_reads):
    d = outdir / "technical_controls"
    d.mkdir(parents=True, exist_ok=True)
    mcol = metric_col()

    for label, df in [
        ("raw_all_genes", raw),
        ("after_reads_and_coverage_filter", support),
    ]:
        (
            df.groupby("biogroup", as_index=False)
            .agg(
                n_genes=("gene_id", "nunique"),
                genes_passing_read_threshold=("mapped_reads",lambda x: int((x >= min_reads).sum())),
                mapped_reads=("mapped_reads", "sum"),
                    sum_HK_median_ratio=(mcol, "sum"),
                mean_HK_median_ratio=(mcol, "mean"),
                median_HK_median_ratio=(mcol, "median"),
                mean_covered_fraction=("covered_fraction", "mean"),
            )
            .to_csv(
                d / f"{prefix}_technical_{label}_biogroup_summary_HKmedian_ratio.tsv",
                sep="\t",
                index=False,
            )
        )

    pd.DataFrame([{
        "method": "median_RPK_ratio",
        "housekeeping_baseline_rpk": baseline,
        "formula": "RPK_gene / median_RPK_housekeeping",
        "scale_factor": 1,
    }]).to_csv(
        d / f"{prefix}_housekeeping_baseline_median_ratio.tsv",
        sep="\t",
        index=False,
    )

def main():
    args = parse_args()
    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)

    # Uses the housekeeping gene-name list
    hk_list_out = out / f"{args.prefix}_housekeeping_gene_names_used.txt"
    hk_gene_names = read_housekeeping_gene_names(args.housekeeping_gene_list)
    hk_list_out.write_text("\n".join(sorted(hk_gene_names)) + "\n")

    hk_annot = read_housekeeping_annotation(args.bact_prokka_annot, hk_gene_names)
    hk_annot.to_csv(out / f"{args.prefix}_housekeeping_annotation_used.tsv", sep="\t", index=False)

    gene_map = make_gene_biogroup_map(args.viral_ffn, args.bact_ffn)
    gene_map.to_csv(out / f"{args.prefix}_gene_biogroup_map.tsv", sep="\t", index=False)

    dup = gene_map["gene_id"][gene_map["gene_id"].duplicated()].drop_duplicates()
    dup_file = out / f"{args.prefix}_duplicated_gene_ids_in_ffn_map.txt"
    dup.to_csv(dup_file, index=False, header=False)

    biogroup_conflicts = (
        gene_map.groupby("gene_id")["biogroup"]
        .nunique()
    )
    biogroup_conflicts = biogroup_conflicts[biogroup_conflicts > 1]

    if not biogroup_conflicts.empty:
        raise ValueError(
            f"{len(biogroup_conflicts)} gene IDs map to more than one biogroup. "
            f"See duplicated IDs in: {dup_file}"
        )

    coverage_df = read_coverage(args.depth)
    raw, support, baseline, hk_records, hk_summary = make_master_expression(
        args.kallisto,
        gene_map,
        coverage_df,
        args.min_reads,
        args.min_coverage,
        hk_annot,
    )

    raw.to_csv(out / f"{args.prefix}_gene_expression_raw_unfiltered_HKmedian_ratio.tsv", sep="\t", index=False)
    support.to_csv(out / f"{args.prefix}_gene_expression_SUPPORTED_reads_coverage_HKmedian_ratio.tsv", sep="\t", index=False)
    hk_records.to_csv(out / f"{args.prefix}_housekeeping_expression_records_USED_FOR_BASELINE.tsv", sep="\t", index=False)
    hk_summary.to_csv(out / f"{args.prefix}_housekeeping_baseline_summary.tsv", sep="\t", index=False)

    sys_map = parse_systems(args.defense_dir)
    ref = build_validated_hmm_reference(
        args.df_ref,
        args.adf_ref,
    )
    hmm_map = parse_hmmer(args.defense_dir, ref)

    method_out = out / f"normalized_{args.prefix}_median_ratio"
    method_out.mkdir(parents=True, exist_ok=True)

    method_cols = [
        "gene_id", "mapped_reads", "length_bp", "eff_length",
        "covered_bases", "positions_in_pseudobam", "covered_fraction", "rpk",
        "HK_norm_RPK_median", "biogroup", "ffn_file",
    ]

    raw[method_cols].to_csv(
        method_out / f"{args.prefix}_raw_all_genes_median_ratio.tsv",
        sep="\t",
        index=False,
    )
    support[method_cols].to_csv(
        method_out / f"{args.prefix}_support_genes_median_ratio.tsv",
        sep="\t",
        index=False,
    )

    technical_controls(
        raw,
        support,
        method_out,
        args.prefix,
        baseline,
        args.min_reads,
    )

    system_out = method_out / "system_level"
    hmmer_out = method_out / "hmmer_level"
    system_out.mkdir(parents=True, exist_ok=True)
    hmmer_out.mkdir(parents=True, exist_ok=True)

    sys_map.to_csv(
        system_out / f"{args.prefix}_defense_gene_map.tsv",
        sep="\t",
        index=False,
    )
    # Merge defense annotations with the reads/coverage-supported table.
    sys_master = sys_map.merge(support, on="gene_id", how="inner")
    write_reports(sys_master, system_out, args.prefix)

    hmm_map.to_csv(
        hmmer_out / f"{args.prefix}_defense_gene_map.tsv",
        sep="\t",
        index=False,
    )
    # Merge defense annotations with the reads/coverage-supported table.
    hmm_master = hmm_map.merge(support, on="gene_id", how="inner")
    write_reports(hmm_master, hmmer_out, args.prefix)

    print("[DONE]")
    print(f"Main output: {out}")
    print(f"Normalization folder: {method_out}")
    print("Formula: HK_norm_RPK_median = RPK_gene / median(RPK_housekeeping)")
    print("Scale factor: 1 (no multiplication by 10^6)")
    print("Selection filters: mapped reads and covered fraction only")
    print(f"Median housekeeping RPK baseline: {baseline}")


if __name__ == "__main__":
    main()
PYCODE

chmod +x "${PY}"

ml python

python3 "${PY}" \
  --prefix "${PREFIX}" \
  --outdir "${OUTDIR}" \
  --kallisto "${KALLISTO_TSV}" \
  --depth "${DEPTH_TSV}" \
  --defense-dir "${DEFENSE_DIR}" \
  --viral-ffn "${VIRAL_FFN_DIR}" \
  --bact-ffn "${BACT_FFN_DIR}" \
  --df-ref "${DF_REF}" \
  --adf-ref "${ADF_REF}" \
  --bact-prokka-annot "${BACT_PROKKA_ANNOT}" \
  --housekeeping-gene-list "${HOUSEKEEPING_GENE_LIST}" \
  --min-reads "${MIN_READS}" \
  --min-coverage "${MIN_COVERAGE}"

echo "[$(date)] Finished successfully."
