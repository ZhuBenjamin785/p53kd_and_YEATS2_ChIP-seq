#!/usr/bin/env python3
                      
"""Create a volcano plot from RNA-seq DESeq2 results."""

from pathlib import Path
import argparse

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


project_dir = Path(__file__).resolve().parent
default_results = project_dir / "rna_seq_dea" / "results.csv"
default_output = project_dir / "rna_seq_dea" / "volcano.png"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results", type=Path, default=default_results)
    parser.add_argument("--output", type=Path, default=default_output)
    parser.add_argument("--padj", type=float, default=0.05)
    parser.add_argument("--min-log2fc", type=float, default=1.0)
    parser.add_argument("--labels", type=int, default=15)
    args = parser.parse_args()

    if not args.results.exists():
        raise SystemExit("DESeq2 results not found: {}".format(args.results))
    df = pd.read_csv(args.results)
    gene_column = next((c for c in ("gene_name", "Geneid", "gene", "symbol") if c in df.columns), None)
    if gene_column is None or not {"padj", "log2FoldChange"}.issubset(df.columns):
        raise SystemExit("Results need a gene column, log2FoldChange, and padj.")

    plot_df = df[[gene_column, "log2FoldChange", "padj"]].copy()
    plot_df["log2FoldChange"] = pd.to_numeric(plot_df["log2FoldChange"], errors="coerce")
    plot_df["padj"] = pd.to_numeric(plot_df["padj"], errors="coerce")
    plot_df = plot_df.dropna(subset=["log2FoldChange", "padj"])
    plot_df["neg_log10_padj"] = -np.log10(plot_df["padj"].clip(lower=np.finfo(float).tiny))
    plot_df["category"] = "not significant"
    significant = plot_df["padj"] <= args.padj
    plot_df.loc[significant & (plot_df["log2FoldChange"] >= args.min_log2fc), "category"] = "upregulated"
    plot_df.loc[significant & (plot_df["log2FoldChange"] <= -args.min_log2fc), "category"] = "downregulated"

    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(9, 7))
    sns.scatterplot(
        data=plot_df, x="log2FoldChange", y="neg_log10_padj", hue="category",
        hue_order=["not significant", "upregulated", "downregulated"],
        palette={"not significant": "lightgray", "upregulated": "firebrick", "downregulated": "royalblue"},
        alpha=0.75, s=22, linewidth=0, ax=ax,
    )
    ax.axvline(args.min_log2fc, color="black", linestyle="--", linewidth=1)
    ax.axvline(-args.min_log2fc, color="black", linestyle="--", linewidth=1)
    ax.axhline(-np.log10(args.padj), color="black", linestyle="--", linewidth=1)
    label_df = plot_df[plot_df["category"] != "not significant"].sort_values("padj").head(args.labels)
    for _, row in label_df.iterrows():
        ax.annotate(str(row[gene_column]), (row["log2FoldChange"], row["neg_log10_padj"]), xytext=(4, 4), textcoords="offset points", fontsize=8)
    ax.set_xlabel("log2 fold change")
    ax.set_ylabel("-log10 adjusted p-value")
    ax.set_title("RNA-seq differential expression")
    fig.tight_layout()
    fig.savefig(args.output, dpi=200)
    plt.close(fig)
    print("Volcano plot written to {}".format(args.output))


if __name__ == "__main__":
    main()
