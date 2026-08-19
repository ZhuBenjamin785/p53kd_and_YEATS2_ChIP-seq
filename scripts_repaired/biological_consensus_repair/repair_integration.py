#!/usr/bin/env python3
"""Regenerate matched RNA/H4K16ac integration, correlations, and Fisher tests.

The script never modifies upstream inputs.  Eligible RNA genes have finite
log2FoldChange and padj values.  Eligible ChIP genes are symbols assigned to a
tested DiffBind peak with finite Fold and FDR.  Gene-level ChIP effects are the
median Fold across all assigned peaks; mean and nearest-TSS summaries are
written as sensitivity analyses.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import false_discovery_control, fisher_exact, pearsonr, spearmanr


ROOT = Path(__file__).resolve().parents[3]


def arguments() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--rna", type=Path, default=ROOT / "shared/rna_seq_dea/shp53_vs_shLacZ_0hr/results.csv")
    p.add_argument("--chip", type=Path, default=ROOT / "shared/chipseq_summary_plots/p53KD/diffbind_peak_annotations.csv")
    p.add_argument("--outdir", type=Path, default=ROOT / "shared/biological_consensus_repaired/integration")
    p.add_argument("--rna-fdr", type=float, default=0.05)
    p.add_argument("--rna-min-abs-lfc", type=float, default=1.0)
    p.add_argument("--chip-fdr", type=float, default=0.05)
    p.add_argument("--expected-all-universe-size", type=int, default=2100,
                   help="fail unless this size is reproduced; use 0 for a labelled sensitivity input")
    return p.parse_args()


def clean_symbol(series: pd.Series) -> pd.Series:
    return series.astype("string").str.strip().str.upper().replace({"": pd.NA, "NA": pd.NA})


def collapse_rna(rna: pd.DataFrame) -> pd.DataFrame:
    """Collapse the five duplicated symbols without selecting on significance."""
    cols = ["log2FoldChange", "stat", "padj", "baseMean"]
    for col in cols:
        rna[col] = pd.to_numeric(rna[col], errors="coerce")
    rna["Gene"] = clean_symbol(rna["gene_name"])
    eligible = rna.dropna(subset=["Gene", "log2FoldChange", "padj"]).copy()
    # Median is deterministic and does not choose the most significant duplicate.
    gene = eligible.groupby("Gene", as_index=False).agg(
        RNA_log2FC=("log2FoldChange", "median"),
        RNA_stat=("stat", "median"),
        RNA_padj=("padj", "min"),
        RNA_baseMean=("baseMean", "sum"),
        RNA_ID_rows=("gene_id", "nunique"),
    )
    return gene


def collapse_chip(chip: pd.DataFrame, promoter_only: bool) -> pd.DataFrame:
    for col in ("Fold", "FDR", "distanceToTSS"):
        chip[col] = pd.to_numeric(chip[col], errors="coerce")
    chip["Gene"] = clean_symbol(chip["SYMBOL"])
    chip = chip.dropna(subset=["Gene", "Fold", "FDR"]).copy()
    if promoter_only:
        chip = chip[chip["annotation"].astype("string").str.contains("promoter", case=False, na=False)].copy()
    nearest_index = (
        chip.assign(abs_distance=chip["distanceToTSS"].abs().fillna(np.inf))
        .sort_values(["Gene", "abs_distance", "FDR", "Fold"], ascending=[True, True, True, False])
        .drop_duplicates("Gene").set_index("Gene")
    )
    grouped = chip.groupby("Gene", as_index=False).agg(
        H4K16ac_Fold_median=("Fold", "median"),
        H4K16ac_Fold_mean=("Fold", "mean"),
        H4K16ac_min_FDR=("FDR", "min"),
        H4K16ac_peak_count=("Fold", "size"),
    )
    grouped["H4K16ac_Fold_nearest_TSS"] = grouped["Gene"].map(nearest_index["Fold"])
    return grouped


def correlation_rows(frame: pd.DataFrame, scope: str) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for summary in ("median", "mean", "nearest_TSS"):
        x = frame[f"H4K16ac_Fold_{summary}"].to_numpy(float)
        y = frame["RNA_log2FC"].to_numpy(float)
        for method, result in (("Pearson", pearsonr(x, y)), ("Spearman", spearmanr(x, y))):
            rows.append({
                "scope": scope, "chip_gene_summary": summary, "method": method,
                "N": len(frame), "correlation": float(result.statistic),
                "p_value": float(result.pvalue),
            })
    return rows


def write_gene_set(path: Path, genes: set[str]) -> None:
    pd.DataFrame({"Gene": sorted(genes)}).to_csv(path, index=False)


def main() -> int:
    args = arguments()
    for path in (args.rna, args.chip):
        if not path.is_file():
            raise SystemExit(f"Missing input: {path}")
    args.outdir.mkdir(parents=True, exist_ok=True)
    tables = args.outdir / "tables"
    sets = args.outdir / "gene_sets"
    fea_sets = args.outdir.parent / "fea"
    tables.mkdir(exist_ok=True)
    sets.mkdir(exist_ok=True)
    fea_sets.mkdir(exist_ok=True)

    rna_raw = pd.read_csv(args.rna)
    chip_raw = pd.read_csv(args.chip)
    required_rna = {"gene_id", "gene_name", "baseMean", "log2FoldChange", "stat", "padj"}
    required_chip = {"SYMBOL", "Fold", "FDR", "annotation", "distanceToTSS"}
    if missing := required_rna - set(rna_raw.columns):
        raise SystemExit(f"RNA input missing: {sorted(missing)}")
    if missing := required_chip - set(chip_raw.columns):
        raise SystemExit(f"ChIP input missing: {sorted(missing)}")

    rna = collapse_rna(rna_raw)
    rna_sig = rna[(rna.RNA_padj <= args.rna_fdr) & (rna.RNA_log2FC.abs() >= args.rna_min_abs_lfc)].copy()
    rna_down = set(rna_sig.loc[rna_sig.RNA_log2FC < 0, "Gene"])
    rna_up = set(rna_sig.loc[rna_sig.RNA_log2FC > 0, "Gene"])
    write_gene_set(fea_sets / "rna_eligible_universe.csv", set(rna.Gene))
    write_gene_set(fea_sets / "rna_up.csv", rna_up)
    write_gene_set(fea_sets / "rna_down.csv", rna_down)

    chip_all = collapse_chip(chip_raw.copy(), False)
    chip_promoter = collapse_chip(chip_raw.copy(), True)

    raw_chip = chip_raw.copy()
    raw_chip["Gene"] = clean_symbol(raw_chip["SYMBOL"])
    raw_chip["Fold"] = pd.to_numeric(raw_chip["Fold"], errors="coerce")
    raw_chip["FDR"] = pd.to_numeric(raw_chip["FDR"], errors="coerce")
    chip_sig = raw_chip.dropna(subset=["Gene", "Fold", "FDR"])
    chip_sig = chip_sig[chip_sig.FDR <= args.chip_fdr]

    def directional_chip_sets(peaks: pd.DataFrame) -> tuple[set[str], set[str], set[str]]:
        losses = set(peaks.loc[peaks.Fold < 0, "Gene"])
        gains = set(peaks.loc[peaks.Fold > 0, "Gene"])
        mixed = losses & gains
        return losses - mixed, gains - mixed, mixed

    loss, gain, ambiguous = directional_chip_sets(chip_sig)
    promoter_chip_sig = chip_sig[
        chip_sig["annotation"].astype("string").str.contains("promoter", case=False, na=False)
    ]
    promoter_loss, promoter_gain, promoter_ambiguous = directional_chip_sets(promoter_chip_sig)

    correlations: list[dict[str, object]] = []
    universe_rows = []
    matched_frames: dict[str, pd.DataFrame] = {}
    for scope, chip in (("all_peaks", chip_all), ("promoter_peaks", chip_promoter)):
        matched = rna.merge(chip, on="Gene", how="inner", validate="one_to_one").sort_values("Gene")
        matched_frames[scope] = matched
        matched.to_csv(tables / f"{scope}__full_eligible_gene_integration.csv", index=False)
        write_gene_set(sets / f"{scope}__eligible_universe.csv", set(matched.Gene))
        correlations.extend(correlation_rows(matched, scope))
        universe_rows.append({
            "scope": scope, "RNA_eligible_symbols": len(rna),
            "ChIP_eligible_symbols": len(chip), "matched_universe_symbols": len(matched),
        })

    all_universe = set(matched_frames["all_peaks"].Gene)
    promoter_universe = set(matched_frames["promoter_peaks"].Gene)
    if args.expected_all_universe_size and len(all_universe) != args.expected_all_universe_size:
        raise SystemExit(
            f"Expected matched universe size {args.expected_all_universe_size:,}, observed {len(all_universe):,}"
        )

    # Significant integrated categories are restricted to the matching universe.
    all_categories = {
        "loss_down": loss & rna_down & all_universe,
        "gain_up": gain & rna_up & all_universe,
        "loss_up": loss & rna_up & all_universe,
        "gain_down": gain & rna_down & all_universe,
    }
    promoter_categories = {
        "loss_down": promoter_loss & rna_down & promoter_universe,
        "gain_up": promoter_gain & rna_up & promoter_universe,
        "loss_up": promoter_loss & rna_up & promoter_universe,
        "gain_down": promoter_gain & rna_down & promoter_universe,
    }
    for scope, categories in (("all_peaks", all_categories), ("promoter_peaks", promoter_categories)):
        for name, genes in categories.items():
            write_gene_set(sets / f"{scope}__{name}.csv", genes)

    tests = []
    for label, chip_set, rna_set in (
        ("H4K16ac loss vs RNA down", loss, rna_down),
        ("H4K16ac gain vs RNA up", gain, rna_up),
    ):
        a_set = chip_set & all_universe
        b_set = rna_set & all_universe
        a = len(a_set & b_set)
        b = len(a_set - b_set)
        c = len(b_set - a_set)
        d = len(all_universe - a_set - b_set)
        result = fisher_exact([[a, b], [c, d]], alternative="greater")
        tests.append({
            "test": label, "alternative": "greater", "universe_genes": len(all_universe),
            "chip_genes": len(a_set), "rna_genes": len(b_set), "overlap": a,
            "chip_and_rna": a, "chip_and_not_rna": b,
            "no_chip_and_rna": c, "no_chip_and_not_rna": d,
            "odds_ratio": float(result.statistic), "fisher_p_value": float(result.pvalue),
        })
    tests_df = pd.DataFrame(tests)
    tests_df["BH_FDR_two_directional_tests"] = false_discovery_control(tests_df.fisher_p_value.to_numpy(), method="bh")
    tests_df.to_csv(tables / "fisher_exact_results.csv", index=False)
    pd.DataFrame(correlations).to_csv(tables / "full_universe_correlation_statistics.csv", index=False)
    pd.DataFrame(universe_rows).to_csv(tables / "matched_universe_audit.csv", index=False)
    pd.DataFrame({"scope": ["all_peaks"] * len(ambiguous), "Gene": sorted(ambiguous)}).to_csv(
        tables / "ambiguous_significant_chip_direction_genes.csv", index=False)
    pd.DataFrame({"scope": ["promoter_peaks"] * len(promoter_ambiguous), "Gene": sorted(promoter_ambiguous)}).to_csv(
        tables / "ambiguous_significant_promoter_chip_direction_genes.csv", index=False)

    category_counts = []
    for scope, categories in (("all_peaks", all_categories), ("promoter_peaks", promoter_categories)):
        for name, genes in categories.items():
            category_counts.append({"scope": scope, "category": name, "count": len(genes)})
    pd.DataFrame(category_counts).to_csv(tables / "significant_category_counts.csv", index=False)

    parameters = {
        "rna_input": str(args.rna.resolve()), "chip_input": str(args.chip.resolve()),
        "rna_fdr": args.rna_fdr, "rna_min_abs_lfc": args.rna_min_abs_lfc,
        "chip_fdr": args.chip_fdr, "fisher_alternative": "greater",
        "primary_chip_gene_summary_for_correlation": "median Fold across all assigned tested peaks",
        "sensitivity_chip_gene_summaries": ["mean", "nearest absolute distanceToTSS"],
        "ambiguous_significant_chip_genes_excluded_from_directional_sets": sorted(ambiguous),
        "ambiguous_significant_promoter_chip_genes_excluded_from_directional_sets": sorted(promoter_ambiguous),
    }
    (args.outdir / "parameters.json").write_text(json.dumps(parameters, indent=2) + "\n")
    print(pd.DataFrame(universe_rows).to_string(index=False))
    print("\nFisher tests")
    print(tests_df.to_string(index=False))
    print("\nPrimary full-universe correlations")
    print(pd.DataFrame(correlations).query("chip_gene_summary == 'median'").to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
