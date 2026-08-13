#!/usr/bin/env python3
                      
"""Direct comparison of p53 KD and YEATS2 KD H4K16ac redistribution."""

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

root = Path.cwd()
plot_root = root / "chipseq_summary_plots"
outdir = plot_root / "p53KD_vs_YEATS2KD"
p53_dir = plot_root / "p53KD"
y2_dir = plot_root / "YEATS2KD"
p53_color, y2_color = "#E45756", "#72B7B2"
gain_color, loss_color = "#D55E00", "#0072B2"

sns.set_theme(style="white", context="notebook", font_scale=1.15)
plt.rcParams.update({"axes.titleweight": "bold", "axes.titlesize": 14,
                     "axes.labelsize": 11, "legend.fontsize": 10,
                     "savefig.bbox": "tight"})


def symmetric_bound(values, minimum=1e-9, pad=1.08):
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    return max(minimum, (np.max(np.abs(values)) if values.size else minimum) * pad)


def load_diff(path):
    df = pd.read_csv(path)
    df["Fold"] = pd.to_numeric(df["Fold"], errors="coerce")
    df["FDR"] = pd.to_numeric(df["FDR"], errors="coerce")
    return df[df.fdr < .05].copy()


def load_gene_effects(path):
    df = pd.read_csv(path)
    fold_col = next((c for c in ["Fold", "log2FC"] if c in df), None)
    fdr_col = next((c for c in ["FDR", "padj"] if c in df), None)
    if not fold_col or not fdr_col or "geneId" not in df:
        raise ValueError(f"Missing Fold/FDR/geneId columns in {path}")
    df = df[(pd.to_numeric(df[fdr_col], errors="coerce") < .05) & df.geneId.notna()].copy()
    df["geneId"] = df.geneId.astype(str).str.replace(r"\.0$", "", regex=True)
    df["Fold_value"] = pd.to_numeric(df[fold_col], errors="coerce")
    return df.groupby("geneId", as_index=False)["Fold_value"].median()


def interval_overlap_counts(left, right):
    overlap_any, concordant = 0, 0
    right_by_chr = {chrom: group.sort_values("start") for chrom, group in right.groupby("seqnames")}
    for row in left.itertuples():
        candidates = right_by_chr.get(row.seqnames)
        if candidates is None: continue
        hits = candidates[(candidates.start <= row.end) & (candidates.end >= row.start)]
        if len(hits):
            overlap_any += 1
            if np.any(np.sign(hits.Fold.to_numpy()) == np.sign(row.Fold)):
                concordant += 1
    return overlap_any, concordant


def main():
    required = [
        root / "diffbind_results/DiffBind_all_peaks.csv",
        root / "diffbind_results_yeats2/DiffBind_YEATS2_all_peaks.csv",
        p53_dir / "diffbind_peak_annotations.csv",
        y2_dir / "diffbind_peak_annotations.csv",
        p53_dir / "TSS_average_profiles.csv",
        y2_dir / "TSS_average_profiles.csv",
    ]
    missing = [str(path) for path in required if not path.is_file()]
    if missing: raise FileNotFoundError("Run both dataset summary jobs first. Missing:\n" + "\n".join(missing))
    outdir.mkdir(parents=True, exist_ok=True)

    p53 = load_diff(required[0]); y2 = load_diff(required[1])
    p53_genes = load_gene_effects(required[2]); y2_genes = load_gene_effects(required[3])
    shared = p53_genes.merge(y2_genes, on="geneId", suffixes=("_p53KD", "_YEATS2KD"))
    shared["Direction"] = np.where(
        np.sign(shared.Fold_value_p53KD) == np.sign(shared.Fold_value_YEATS2KD),
        "Concordant", "Discordant")
    shared.to_csv(outdir / "shared_nearest_gene_direction_concordance.csv", index=False)

    p53_gene_set, y2_gene_set = set(p53_genes.geneId), set(y2_genes.geneId)
    gene_counts = pd.DataFrame({
        "Category": ["p53 KD only", "Shared", "YEATS2 KD only"],
        "Genes": [len(p53_gene_set - y2_gene_set), len(p53_gene_set & y2_gene_set), len(y2_gene_set - p53_gene_set)]})
    gene_counts.to_csv(outdir / "nearest_gene_overlap_summary.csv", index=False)
    overlap_any, overlap_concordant = interval_overlap_counts(p53, y2)
    pd.DataFrame({"Metric": ["Significant p53 KD peaks", "Significant YEATS2 KD peaks",
                              "p53 KD peaks overlapping a YEATS2 KD peak", "Overlaps with concordant direction"],
                  "Count": [len(p53), len(y2), overlap_any, overlap_concordant]}).to_csv(
                      outdir / "exact_peak_overlap_summary.csv", index=False)

    fig, axes = plt.subplots(2, 2, figsize=(13, 10))
    effect_counts = pd.DataFrame({
        "Dataset": ["p53 KD", "p53 KD", "YEATS2 KD", "YEATS2 KD"],
        "Direction": ["Gain", "Loss", "Gain", "Loss"],
        "Peaks": [(p53.Fold > 0).sum(), (p53.Fold < 0).sum(), (y2.Fold > 0).sum(), (y2.Fold < 0).sum()]})
    sns.barplot(data=effect_counts, x="Dataset", y="Peaks", hue="Direction",
                palette={"Gain": gain_color, "Loss": loss_color}, ax=axes[0, 0])
    axes[0, 0].set_title("Significant H4K16ac redistribution"); axes[0, 0].legend(frameon=False)

    sns.barplot(data=gene_counts, x="Genes", y="Category",
                hue="Category", palette=[p53_color, "#8F63B8", y2_color], legend=False, ax=axes[0, 1])
    axes[0, 1].set_title("Nearest-gene overlap (FDR < 0.05)"); axes[0, 1].set_ylabel(None)

    palette = {"Concordant": "#4D9221", "Discordant": "#C51B7D"}
    sns.scatterplot(data=shared, x="Fold_value_p53KD", y="Fold_value_YEATS2KD",
                    hue="Direction", palette=palette, alpha=.65, s=28, ax=axes[1, 0])
    axes[1, 0].axhline(0, color="0.5", ls="--", lw=1); axes[1, 0].axvline(0, color="0.5", ls="--", lw=1)
    scatter_limit = symmetric_bound(np.concatenate([
        shared.Fold_value_p53KD.to_numpy(), shared.Fold_value_YEATS2KD.to_numpy()
    ]))
    axes[1, 0].set_xlim(-scatter_limit, scatter_limit)
    axes[1, 0].set_ylim(-scatter_limit, scatter_limit)
    axes[1, 0].set_aspect("equal", adjustable="box")
    axes[1, 0].set(title="Shared nearest-gene direction concordance",
                   xlabel="p53 KD DiffBind log2 fold-change", ylabel="YEATS2 KD DiffBind log2 fold-change")
    axes[1, 0].legend(frameon=False)

    p53_prof = pd.read_csv(required[4]); y2_prof = pd.read_csv(required[5])
    axes[1, 1].plot(p53_prof.distance_bp, p53_prof["p53 KD"] - p53_prof["Scramble"],
                    color=p53_color, lw=2.3, label="p53 KD − Scramble")
    axes[1, 1].plot(y2_prof.distance_bp, y2_prof["YEATS2 KD"] - y2_prof["Scramble"],
                    color=y2_color, lw=2.3, label="YEATS2 KD − Scramble")
    axes[1, 1].axhline(0, color="0.5", ls="--", lw=1); axes[1, 1].axvline(0, color="0.5", ls=":", lw=1)
    axes[1, 1].set_xlim(-2000, 2000)
    axes[1, 1].set(title="TSS-centered H4K16ac redistribution",
                   xlabel="Distance from TSS (bp)", ylabel="Δ mean log2(ChIP/Input)")
    axes[1, 1].legend(frameon=False)
    for ax in axes.flat: sns.despine(ax=ax)
    fig.suptitle("p53 KD versus YEATS2 KD: redistributed H4K16ac", fontsize=17, fontweight="bold", y=1.01)
    fig.tight_layout(); fig.savefig(outdir / "p53KD_vs_YEATS2KD_summary.pdf"); plt.close(fig)

    (outdir / "figure_notes.txt").write_text(
        "Differential peaks are DiffBind results at FDR < 0.05.\n"
        "Nearest genes are assigned by ChIPseeker; gene overlap does not imply exact peak overlap.\n"
        "TSS profiles show KD minus Scramble from log2(ChIP/Input) tracks generated from dm6-normalized BAMs.\n"
        "The comparison emphasizes H4K16ac redistribution rather than requiring exact peak overlap.\n")


if __name__ == "__main__":
    main()
